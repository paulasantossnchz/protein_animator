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
import requests
from Bio.Blast import NCBIWWW, NCBIXML


# ── Module-level constants ────────────────────────────────────────────────────

_DEFAULT_OUTPUT_DIR  = os.path.join("data", "input_pdbs")
_RCSB_FASTA_URL      = "https://www.rcsb.org/fasta/entry/{code}/download"
_RCSB_SEARCH_URL     = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_DOWNLOAD_URL   = "https://files.rcsb.org/download/{code}.pdb"


# ── Auxiliary functions ───────────────────────────────────────────────────────

def _extract_pdb_code(hit_id: str) -> str:
    """
    Return the 4-letter PDB accession from a BLAST hit identifier.

    BLAST returns PDB hits with IDs like 'pdb|1ABC|A'; this function
    isolates the accession portion regardless of the number of '|' fields.
    """
    parts = hit_id.split("|")
    if len(parts) >= 2:
        code = parts[1][:4].upper()
    else:
        code = hit_id[:4].upper()
    return code


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

    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            sequence = _parse_first_fasta_sequence(response.text)
        else:
            print(
                f"[_fetch_pdb_sequence] HTTP {response.status_code} "
                f"while fetching FASTA for {pdb_code}."
            )
    except Exception as exc:
        print(f"[_fetch_pdb_sequence] Connection error for {pdb_code}: {exc}")

    return sequence


def _build_rcsb_homolog_query(
    sequence: str,
    identity_cutoff: float,
    evalue_cutoff: float,
    resolution_cutoff: float,
) -> dict:
    """
    Build the RCSB Search API v2 JSON query that combines sequence similarity
    with a crystallographic resolution filter.

    identity_cutoff  : value in [0, 1] (e.g. 0.30 for 30 %).
    evalue_cutoff    : maximum E-value accepted by the sequence service.
    resolution_cutoff: maximum resolution in Å for the attribute filter.
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
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
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
                        pdb_code = _extract_pdb_code(alignment.hit_id)
                        if pdb_code not in filtered_pdbs:
                            filtered_pdbs.append(pdb_code)

    except Exception as exc:
        print(f"[blast_fasta_vs_pdb] BLAST error: {exc}")
        filtered_pdbs = []

    return filtered_pdbs


def search_homologs_rcsb(
    pdb_code: str,
    identity_cutoff: float  = 0.30,
    evalue_cutoff: float    = 1.0,
    resolution_cutoff: float = 3.0,
) -> list:
    """
    Query the RCSB PDB Search API for structures homologous to *pdb_code*.

    The function first retrieves the sequence of the first chain of *pdb_code*
    from RCSB, then issues a combined query that enforces both sequence
    similarity and a crystallographic resolution limit.

    Parameters
    ----------
    pdb_code          : 4-letter PDB identifier of the reference structure.
    identity_cutoff   : Minimum sequence identity in [0, 1] (default 0.30 = 30 %).
    evalue_cutoff     : Maximum E-value for the sequence search (default 1.0).
    resolution_cutoff : Maximum resolution in Å for returned structures (default 3.0 Å).

    Returns
    -------
    List of PDB entry identifiers (excluding *pdb_code* itself), or an empty
    list if the sequence cannot be fetched or the API call fails.
    """
    homologs = []

    sequence = _fetch_pdb_sequence(pdb_code)

    if sequence:
        try:
            query    = _build_rcsb_homolog_query(
                sequence, identity_cutoff, evalue_cutoff, resolution_cutoff
            )
            response = requests.post(_RCSB_SEARCH_URL, json=query, timeout=30)

            if response.status_code == 200:
                data = response.json()
                for result in data.get("result_set", []):
                    identifier = result.get("identifier", "").upper()
                    is_query   = (identifier == pdb_code.upper())
                    if identifier and not is_query:
                        homologs.append(identifier)
            else:
                print(
                    f"[search_homologs_rcsb] RCSB Search API returned "
                    f"HTTP {response.status_code} for query on {pdb_code}."
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
    pdb_codes: list,
    output_dir: str = _DEFAULT_OUTPUT_DIR,
) -> list:
    """
    Download a .pdb file for each code in *pdb_codes* and save to *output_dir*.

    Files are opened with the classic open()/close() pattern (no context
    manager). A try/finally block guarantees the file handle is always closed
    even if the write raises an unexpected exception.

    Parameters
    ----------
    pdb_codes  : Iterable of 4-letter PDB identifiers.
    output_dir : Destination directory; created automatically if absent.

    Returns
    -------
    List of PDB codes that were successfully written to disk.
    """
    os.makedirs(output_dir, exist_ok=True)
    downloaded = []

    for code in pdb_codes:
        code_upper = code.upper()
        url        = _RCSB_DOWNLOAD_URL.format(code=code_upper)
        filepath   = os.path.join(output_dir, f"{code_upper}.pdb")

        try:
            response = requests.get(url, timeout=30)

            if response.status_code == 200:
                archivo = open(filepath, "w", encoding="utf-8")
                try:
                    archivo.write(response.text)
                finally:
                    archivo.close()

                downloaded.append(code_upper)
                print(f"[download_pdbs] {code_upper} → {filepath}")

            else:
                print(
                    f"[download_pdbs] HTTP {response.status_code} "
                    f"for {code_upper}; skipping."
                )

        except Exception as exc:
            print(f"[download_pdbs] Error downloading {code_upper}: {exc}")

    return downloaded
