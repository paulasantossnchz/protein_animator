"""
app.py
Gradio-based interactive interface for the Protein Structure Animator pipeline.

Connects the three backend modules in sequence:
  pdb_fetcher  → retrieve / download PDB structures
  interpolator → align structures and build the morph trajectory in PyMOL
  animator     → render each frame and export an animated GIF

Run with:
    python app.py
"""

import os
import gradio as gr

from src.pdb_fetcher  import blast_fasta_vs_pdb, search_homologs_rcsb, download_pdbs
from src.interpolator import interpolate_structures
from src.animator     import animate_morph


# ── Constants ─────────────────────────────────────────────────────────────────

_INPUT_DIR  = os.path.join("data", "input_pdbs")
_OUTPUT_DIR = os.path.join("data", "output_animations")


# ── Pipeline helpers ──────────────────────────────────────────────────────────

def _codes_from_fasta(
    fasta_seq: str,
    identity_threshold: float,
    coverage_threshold: float,
) -> tuple:
    """
    Run BLASTP against the PDB database and return (pdb_codes, status_msg).

    Returns an empty list with an explanatory message if no hits pass the
    thresholds or if BLAST itself fails.
    """
    pdb_codes = []
    status    = ""

    result = blast_fasta_vs_pdb(
        fasta_seq.strip(),
        identity_threshold=identity_threshold,
        coverage_threshold=coverage_threshold,
    )

    if len(result) > 0:
        pdb_codes = result
        status    = f"BLAST: {len(pdb_codes)} hit(s) encontrados."
    else:
        status = (
            "BLAST no encontró hits por encima de los umbrales seleccionados. "
            "Prueba a reducir la identidad o la cobertura mínima."
        )

    return pdb_codes, status


def _codes_from_pdb(
    pdb_code_entry: str,
    identity_threshold: float,
    resolution_cutoff: float,
) -> tuple:
    """
    Search RCSB for structures homologous to *pdb_code_entry* and return
    (pdb_entries, status_msg).

    pdb_code_entry may be a plain 4-letter code ('1CLL') or include a chain
    ('1CLL_A').  The chain is preserved for the reference entry; homologs
    returned by RCSB default to chain A.
    identity_threshold is converted from percentage to fraction [0, 1].
    """
    pdb_entries       = []
    status            = ""
    parts             = pdb_code_entry.upper().split("_")
    pdb_code          = parts[0][:4]
    ref_chain         = parts[1] if len(parts) >= 2 else "A"
    ref_entry         = f"{pdb_code}_{ref_chain}"
    identity_fraction = identity_threshold / 100.0

    homologs = search_homologs_rcsb(
        pdb_code,
        identity_cutoff=identity_fraction,
        resolution_cutoff=resolution_cutoff,
    )

    pdb_entries = [ref_entry] + homologs

    if len(homologs) > 0:
        status = (
            f"RCSB Search: {len(homologs)} homólogo(s) encontrado(s) para '{ref_entry}'. "
            f"Total de estructuras: {len(pdb_entries)}."
        )
    else:
        status = (
            f"No se encontraron homólogos para '{ref_entry}' con los parámetros actuales. "
            "Solo hay una estructura; se necesitan al menos 2 para animar."
        )

    return pdb_entries, status


def run_pipeline(
    fasta_seq: str,
    pdb_code: str,
    identity_threshold: float,
    coverage_threshold: float,
    resolution_cutoff: float,
    morph_steps: int,
    fps: int,
) -> tuple:
    """
    Orchestrate the full pipeline and return (pdb_list_text, status_msg, gif_path).

    Priority: FASTA sequence > PDB code.  If both fields are empty an error
    message is returned immediately without running any computation.

    All failures are communicated through the returned status string rather
    than exceptions, so Gradio can display them in the UI.
    """
    pdb_list_text = ""
    status        = ""
    gif_path      = None

    has_fasta    = bool(fasta_seq and fasta_seq.strip())
    has_pdb_code = bool(pdb_code  and pdb_code.strip())
    input_valid  = has_fasta or has_pdb_code

    pdb_codes    = []
    codes_status = ""

    if not input_valid:
        status = (
            "Error: introduce una secuencia FASTA (pestaña 'Secuencia FASTA') "
            "o un código PDB (pestaña 'Código PDB') antes de continuar."
        )

    elif has_fasta:
        pdb_codes, codes_status = _codes_from_fasta(
            fasta_seq, identity_threshold, coverage_threshold
        )

    else:
        pdb_codes, codes_status = _codes_from_pdb(
            pdb_code.strip(), identity_threshold, resolution_cutoff
        )

    enough_codes = input_valid and len(pdb_codes) >= 2

    if input_valid and not enough_codes:
        status = (
            codes_status
            + "  Se necesitan al menos 2 estructuras para generar la animación."
        )

    if enough_codes:
        pdb_list_text = "\n".join(pdb_codes)
        status        = codes_status + "  Descargando estructuras…"

        pdb_entries = download_pdbs(pdb_codes, output_dir=_INPUT_DIR)

        if len(pdb_entries) < 2:
            status = (
                f"Error: solo se descargaron {len(pdb_entries)} estructura(s) de "
                f"{len(pdb_codes)} solicitada(s). Se necesitan al menos 2."
            )

        else:
            status = (
                codes_status
                + f"  {len(pdb_entries)} estructura(s) descargadas. "
                + "Generando trayectoria de morphing…"
            )
            morph_name = interpolate_structures(pdb_entries, morph_steps=int(morph_steps))

            if morph_name:
                status = (
                    codes_status
                    + f"  Trayectoria '{morph_name}' lista. Exportando GIF…"
                )
                gif_path = animate_morph(
                    morph_name, output_dir=_OUTPUT_DIR, fps=int(fps)
                )

                if gif_path:
                    status = (
                        f"Animación generada correctamente.\n"
                        f"Archivo: {gif_path}\n"
                        f"Estructuras usadas: {', '.join(pdb_codes)}"
                    )
                else:
                    status = (
                        "Error: la exportación del GIF falló. "
                        "Consulta los mensajes de la consola para más detalles."
                    )

            else:
                status = (
                    "Error: la interpolación con PyMOL falló. "
                    "Comprueba que las estructuras descargadas comparten suficientes "
                    "residuos equivalentes para alinearse."
                )

    return pdb_list_text, status, gif_path


# ── Gradio interface ──────────────────────────────────────────────────────────

def build_interface() -> gr.Blocks:
    """
    Assemble and return the Gradio Blocks application.

    All widget definitions, layout, and event wiring are contained here.
    No side-effects on module import.
    """
    with gr.Blocks(title="Protein Structure Animator") as demo:

        gr.Markdown(
            """
            # 🧬 Protein Structure Animator
            Genera animaciones de *morphing* entre estructuras de proteínas homólogas.

            **Modo FASTA** — ejecuta BLAST contra la base de datos PDB y filtra por identidad y cobertura.
            **Modo PDB** — busca homólogos directamente en RCSB y filtra por identidad y resolución.
            """
        )

        # ── Input block ───────────────────────────────────────────────────────
        gr.Markdown("### 1. Fuente de datos")

        with gr.Tabs():
            with gr.TabItem("Secuencia FASTA"):
                fasta_input = gr.Textbox(
                    label="Secuencia FASTA",
                    placeholder=(
                        ">Mi_proteina\n"
                        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTT..."
                    ),
                    lines=8,
                    info=(
                        "Pega aquí tu secuencia en formato FASTA o como cadena cruda "
                        "de aminoácidos. Se lanzará BLASTP contra la base de datos PDB."
                    ),
                )

            with gr.TabItem("Código PDB"):
                pdb_input = gr.Textbox(
                    label="Código PDB de referencia",
                    placeholder="Ej: 1CRN o 1CRN_A",
                    max_lines=1,
                    info=(
                        "Introduce un código PDB de 4 letras, opcionalmente con cadena "
                        "(p.ej. 1CLL_A). Se buscarán homólogos mediante la API de RCSB."
                    ),
                )

        # ── Parameter block ───────────────────────────────────────────────────
        gr.Markdown("### 2. Parámetros de búsqueda e interpolación")

        with gr.Row():
            identity_slider = gr.Slider(
                minimum=10, maximum=100, value=30, step=5,
                label="Identidad mínima (%)",
                info="Umbral de identidad de secuencia para conservar un hit.",
            )
            coverage_slider = gr.Slider(
                minimum=10, maximum=100, value=50, step=5,
                label="Cobertura mínima (%) — solo BLAST",
                info="Porcentaje de la query cubierto por el alineamiento.",
            )

        with gr.Row():
            resolution_slider = gr.Slider(
                minimum=1.0, maximum=5.0, value=3.0, step=0.5,
                label="Resolución máxima (Å) — solo RCSB",
                info="Filtra estructuras con resolución peor que este valor.",
            )
            morph_steps_slider = gr.Slider(
                minimum=5, maximum=100, value=30, step=5,
                label="Pasos de interpolación",
                info="Fotogramas interpolados entre cada par de estructuras consecutivas.",
            )

        with gr.Row():
            fps_slider = gr.Slider(
                minimum=5, maximum=60, value=25, step=5,
                label="FPS de la animación",
                info="Velocidad de reproducción del GIF de salida.",
            )

        # ── Action ────────────────────────────────────────────────────────────
        gr.Markdown("### 3. Ejecutar")

        generate_btn = gr.Button(
            value="▶  Generar Animación",
            variant="primary",
            size="lg",
        )

        # ── Output block ──────────────────────────────────────────────────────
        gr.Markdown("### 4. Resultados")

        with gr.Row():
            pdb_list_output = gr.Textbox(
                label="Códigos PDB recuperados",
                lines=5,
                interactive=False,
                placeholder="Los códigos PDB seleccionados aparecerán aquí una vez completada la búsqueda.",
            )
            status_output = gr.Textbox(
                label="Estado del pipeline",
                lines=5,
                interactive=False,
                placeholder="El progreso y los mensajes de error del pipeline aparecerán aquí.",
            )

        gif_output = gr.Image(
            label="Animación generada",
            type="filepath",
        )

        # ── Event wiring ──────────────────────────────────────────────────────
        generate_btn.click(
            fn=run_pipeline,
            inputs=[
                fasta_input,
                pdb_input,
                identity_slider,
                coverage_slider,
                resolution_slider,
                morph_steps_slider,
                fps_slider,
            ],
            outputs=[pdb_list_output, status_output, gif_output],
        )

    return demo


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    app = build_interface()
    app.launch(share=False, server_name="0.0.0.0", server_port=7860, theme=gr.themes.Soft())
