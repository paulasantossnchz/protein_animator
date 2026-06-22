"""
pdb_fetcher.py
Automates retrieval of homologous PDB codes via NCBI BLAST and the RCSB Search
API, and downloads the corresponding structure files to disk.

Public API
----------
blast_fasta_vs_pdb   : BLASTP a FASTA sequence against the PDB database.
search_homologs_rcsb : Find homologs for a known PDB code via RCSB Search API.
download_pdbs        : Download .pdb files for a list of PDB codes.
"""

import os
import time
import requests
from Bio.Blast import NCBIWWW, NCBIXML


# ── Module-level constants ────────────────────────────────────────────────────

_DEFAULT_OUTPUT_DIR  = os.path.join("data", "input_pdbs")
_RCSB_FASTA_URL      = "https://www.rcsb.org/fasta/entry/{code}/download"
_RCSB_SEARCH_URL     = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_DOWNLOAD_URL   = "https://files.rcsb.org/download/{code}.pdb"
_MAX_RETRIES         = 3
_RETRY_BACKOFF       = 2.0   # segundos base entre reintentos (escala con el intento)


# ── Auxiliary functions ───────────────────────────────────────────────────────

def _http_with_retries(method: str, url: str, timeout: int, retries: int = _MAX_RETRIES, **kwargs):
    """
    Perform an HTTP request retrying on transient failures.

    The RCSB endpoints occasionally answer a successful query with a sporadic
    HTTP 4xx/5xx or drop the connection; retrying the identical request a few
    times with a short, increasing back-off recovers from those glitches and
    prevents the pipeline from wrongly reporting "no homologs".

    Returns the last Response obtained (its status_code is 200 on success), or
    None if every attempt raised a connection error.
    """
    response  = None
    attempt   = 0
    succeeded = False

    while attempt < retries and not succeeded:
        try:
            response = requests.request(method, url, timeout=timeout, **kwargs)
            if response.status_code == 200:
                succeeded = True
            else:
                print(
                    f"[_http_with_retries] HTTP {response.status_code} en "
                    f"{url} (intento {attempt + 1}/{retries})."
                )
        except Exception as exc:
            print(
                f"[_http_with_retries] Error de conexión en {url} "
                f"(intento {attempt + 1}/{retries}): {exc}"
            )
            response = None

        is_last = (attempt >= retries - 1)
        if not succeeded and not is_last:
            time.sleep(_RETRY_BACKOFF * (attempt + 1))

        attempt = attempt + 1

    return response

def _extract_pdb_entry(hit_id: str) -> str:
    """
    Return a CODE_CHAIN string from a BLAST hit identifier.

    BLAST returns PDB hits with IDs like 'pdb|1ABC|A'; this function
    extracts both the accession and the chain letter.
    Chain defaults to 'A' when absent from the identifier.
    """
    parts = hit_id.split("|")
    code  = parts[1][:4].upper() if len(parts) >= 2 else hit_id[:4].upper()
    chain = parts[2].strip().upper() if (len(parts) >= 3 and parts[2].strip()) else "A"
    return f"{code}_{chain}"


def _parse_first_fasta_sequence(fasta_text: str) -> str:
    """
    Extract the amino-acid sequence of the first entry in a FASTA block.

    Stops collecting residues when a second header line ('>') is encountered,
    using boolean flags instead of break/continue to control the loop.
    """
    lines = fasta_text.strip().split("\n")
    sequence_parts = []
    first_header_found  = False
    second_header_found = False

    for line in lines:
        is_header = line.startswith(">")

        if is_header and not first_header_found:
            first_header_found = True
        elif is_header and first_header_found:
            second_header_found = True

        if first_header_found and not second_header_found and not is_header:
            sequence_parts.append(line.strip())

    return "".join(sequence_parts)


def _fetch_pdb_sequence(pdb_code: str) -> str:
    """
    Download the FASTA of the first chain for *pdb_code* from RCSB.

    Returns an empty string and prints the error if the request fails.
    """
    sequence = ""
    url = _RCSB_FASTA_URL.format(code=pdb_code.upper())

    response = _http_with_retries("GET", url, 15)
    if response is not None and response.status_code == 200:
        sequence = _parse_first_fasta_sequence(response.text)
    else:
        print(f"[_fetch_pdb_sequence] No se pudo obtener la FASTA de {pdb_code}.")

    return sequence


def _build_rcsb_homolog_query(
    sequence: str,
    identity_cutoff: float,
    evalue_cutoff: float,
    resolution_cutoff: float,
    max_results: int,
) -> dict:
    """
    Build the RCSB Search API v2 JSON query that combines sequence similarity
    with a crystallographic resolution filter.

    The return type is "polymer_instance" (not "entry") so the API reports the
    exact chain that matched the sequence search, e.g. "1G4Y.B".  Assuming
    chain A would silently load the wrong polymer in complexes where the
    homolog sits on chain B/C/…, which breaks the downstream atom matching.

    identity_cutoff  : value in [0, 1] (e.g. 0.30 for 30 %).
    evalue_cutoff    : maximum E-value accepted by the sequence service.
    resolution_cutoff: maximum resolution in Å for the attribute filter.
    max_results      : maximum number of chain hits to request (top by score).
    """
    sequence_node = {
        "type": "terminal",
        "service": "sequence",
        "parameters": {
            "evalue_cutoff":    evalue_cutoff,
            "identity_cutoff":  identity_cutoff,
            "sequence_type":    "protein",
            "value":            sequence,
        },
    }

    resolution_node = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator":  "less_or_equal",
            "value":     resolution_cutoff,
            "negation":  False,
        },
    }

    query = {
        "query": {
            "type":             "group",
            "logical_operator": "and",
            "nodes":            [sequence_node, resolution_node],
        },
        "return_type": "polymer_instance",
        "request_options": {
            "paginate": {"start": 0, "rows": max_results},
            "sort":     [{"sort_by": "score", "direction": "desc"}],
        },
    }

    return query


# ── Public API ────────────────────────────────────────────────────────────────

def blast_fasta_vs_pdb(
    sequence: str,
    identity_threshold: float = 30.0,
    coverage_threshold: float  = 50.0,
) -> list:
    """
    Run BLASTP against the NCBI PDB database and return filtered PDB codes.

    Parameters
    ----------
    sequence           : Raw amino-acid sequence or FASTA-formatted string.
    identity_threshold : Minimum % sequence identity to retain a hit (default 30 %).
    coverage_threshold : Minimum query-coverage percentage to retain a hit (default 50 %).

    Returns
    -------
    Deduplicated list of 4-letter PDB codes passing both thresholds, or an
    empty list if BLAST fails or no hits survive filtering.
    """
    filtered_pdbs = []

    try:
        result_handle  = NCBIWWW.qblast("blastp", "pdb", sequence)
        blast_records  = list(NCBIXML.parse(result_handle))
        result_handle.close()

        for record in blast_records:
            query_length = max(record.query_length, 1)

            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    identity_pct = (hsp.identities / hsp.align_length) * 100
                    coverage_pct = (
                        (hsp.query_end - hsp.query_start) / query_length
                    ) * 100

                    passes_filter = (
                        identity_pct >= identity_threshold
                        and coverage_pct >= coverage_threshold
                    )

                    if passes_filter:
                        pdb_entry = _extract_pdb_entry(alignment.hit_id)
                        if pdb_entry not in filtered_pdbs:
                            filtered_pdbs.append(pdb_entry)

    except Exception as exc:
        print(f"[blast_fasta_vs_pdb] BLAST error: {exc}")
        filtered_pdbs = []

    return filtered_pdbs


def search_homologs_rcsb(
    pdb_code: str,
    identity_cutoff: float  = 0.30,
    evalue_cutoff: float    = 1.0,
    resolution_cutoff: float = 3.0,
    max_results: int        = 15,
) -> list:
    """
    Query the RCSB PDB Search API for structures homologous to *pdb_code*.

    The function first retrieves the sequence of the first chain of *pdb_code*
    from RCSB, then issues a combined query that enforces both sequence
    similarity and a crystallographic resolution limit.

    Hits are returned as 'CODE_CHAIN' strings using the chain that actually
    matched the search (e.g. '1G4Y_B'), the query structure is excluded, and
    only the first matching chain of each PDB entry is kept (so a structure
    with several homologous copies does not flood the animation with near
    duplicates).

    Parameters
    ----------
    pdb_code          : 4-letter PDB identifier of the reference structure.
    identity_cutoff   : Minimum sequence identity in [0, 1] (default 0.30 = 30 %).
    evalue_cutoff     : Maximum E-value for the sequence search (default 1.0).
    resolution_cutoff : Maximum resolution in Å for returned structures (default 3.0 Å).
    max_results       : Maximum number of chain hits to request (default 15).

    Returns
    -------
    List of 'CODE_CHAIN' identifiers (excluding *pdb_code* itself), or an empty
    list if the sequence cannot be fetched or the API call fails.
    """
    homologs   = []
    seen_codes = set()

    sequence = _fetch_pdb_sequence(pdb_code)

    if sequence:
        try:
            query    = _build_rcsb_homolog_query(
                sequence, identity_cutoff, evalue_cutoff,
                resolution_cutoff, max_results,
            )
            response = _http_with_retries("POST", _RCSB_SEARCH_URL, 30, json=query)

            if response is not None and response.status_code == 200:
                data = response.json()
                for result in data.get("result_set", []):
                    # polymer_instance identifiers look like '1G4Y.B'
                    identifier = result.get("identifier", "").upper().replace(".", "_")
                    code       = identifier.split("_")[0]
                    is_query   = (code == pdb_code.upper()[:4])
                    is_new     = (code not in seen_codes)
                    if identifier and not is_query and is_new:
                        seen_codes.add(code)
                        homologs.append(identifier)
            else:
                status_code = response.status_code if response is not None else "sin respuesta"
                print(
                    f"[search_homologs_rcsb] La búsqueda RCSB falló tras "
                    f"{_MAX_RETRIES} intentos (HTTP {status_code}) para {pdb_code}."
                )

        except Exception as exc:
            print(f"[search_homologs_rcsb] Error querying RCSB for {pdb_code}: {exc}")
            homologs = []

    else:
        print(
            f"[search_homologs_rcsb] Empty sequence for {pdb_code}; "
            "search aborted."
        )

    return homologs


def download_pdbs(
    pdb_entries: list,
    output_dir: str = _DEFAULT_OUTPUT_DIR,
) -> list:
    """
    Download a .pdb file for each entry in *pdb_entries* and save to *output_dir*.

    Each entry may be a plain 4-letter code ('1ABC') or include a chain
    ('1ABC_A').  The chain is preserved and paired with the downloaded file
    path so the caller can later filter atoms to the correct chain.

    Files are opened with the classic open()/close() pattern (no context
    manager). A try/finally block guarantees the file handle is always closed
    even if the write raises an unexpected exception.

    Parameters
    ----------
    pdb_entries : Iterable of identifiers in '1ABC' or '1ABC_A' format.
    output_dir  : Destination directory; created automatically if absent.

    Returns
    -------
    List of (filepath, chain) tuples for entries successfully written to disk.
    """
    os.makedirs(output_dir, exist_ok=True)
    downloaded = []

    for entry in pdb_entries:
        parts      = entry.upper().split("_")
        code_upper = parts[0][:4]
        chain      = parts[1] if len(parts) >= 2 else "A"
        url        = _RCSB_DOWNLOAD_URL.format(code=code_upper)
        filepath   = os.path.join(output_dir, f"{code_upper}.pdb")

        try:
            response = _http_with_retries("GET", url, 30)

            if response is not None and response.status_code == 200:
                archivo = open(filepath, "w", encoding="utf-8")
                try:
                    archivo.write(response.text)
                finally:
                    archivo.close()

                downloaded.append((filepath, chain))
                print(f"[download_pdbs] {code_upper} chain {chain} → {filepath}")

            else:
                status_code = response.status_code if response is not None else "sin respuesta"
                print(
                    f"[download_pdbs] HTTP {status_code} "
                    f"for {code_upper}; skipping."
                )

        except Exception as exc:
            print(f"[download_pdbs] Error downloading {code_upper}: {exc}")

    return downloaded
