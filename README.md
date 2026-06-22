# 🧬 Protein Structure Animator

Herramienta que genera una **animación** (GIF) del posible cambio conformacional de una
proteína a partir de varias de sus estructuras experimentales. Cada estructura homóloga
es un fotograma; ordenadas por similitud e interpoladas, ofrecen una aproximación visual
a los movimientos que la proteína podría experimentar *in vivo*.

La búsqueda de estructuras es **automática**: basta con una secuencia (vía BLAST) o un
código PDB (vía la API de búsqueda de RCSB).

---

## Características

- **Dos modos de entrada**: secuencia FASTA (BLASTP contra el PDB) o código PDB de
  referencia (búsqueda de homólogos en RCSB).
- **Selección automática de la cadena correcta** (no asume la cadena A).
- **Filtros** de identidad de secuencia, cobertura y resolución máxima.
- **Ordenación por RMSD** para que la transición entre fotogramas sea continua.
- **Interpolación lineal** de coordenadas para suavizar los saltos entre estructuras muy
  diferentes.
- **Reintentos automáticos** ante fallos transitorios de la API de RCSB.
- Interfaz web con **Gradio**; salida en GIF animado.

---

## Requisitos e instalación

```bash
pip install -r requirements.txt
```

Dependencias: `biopython`, `requests`, `numpy`, `pillow`, `gradio`, `pymol-open-source`.

> **PyMOL**: si la instalación de `pymol-open-source` por pip da problemas, suele ser más
> fácil mediante conda:
> ```bash
> conda install -c conda-forge pymol-open-source
> ```

---

## Uso

```bash
python app.py
```

Abre la interfaz en <http://localhost:7860> y elige uno de los dos modos:

- **Pestaña «Código PDB»** — introduce un código de 4 letras (`1CLL`) o con cadena
  (`1CLL_A`). La búsqueda en RCSB es casi instantánea.
- **Pestaña «Secuencia FASTA»** — pega una secuencia. Lanza BLASTP contra el PDB
  (puede tardar 1–3 min por la cola del servidor de NCBI).

Ajusta los parámetros (identidad, cobertura, resolución, pasos de interpolación, FPS),
pulsa **«Generar Animación»** y el GIF aparecerá en el panel de resultados (y se guarda
en `data/output_animations/`).

> El control de **identidad mínima** arranca en el 30 %. Subirlo restringe la búsqueda a
> estructuras más parecidas (menos, pero de la misma proteína); bajarlo incluye homólogos
> más lejanos. Para pruebas rápidas, baja los *pasos de interpolación* a 10–15.

---

## Ejemplos verificados

### Modo PDB

| Código | Identidad | Comentario |
|---|:---:|---|
| `1CLL` | 70 % | Calmodulina: cambio conformacional **muy marcado** (extendida → compacta). Mejor ejemplo. |
| `2HHB` | 90 % | Hemoglobina (cambio cuaternario sutil). |
| `1UBQ` | 90 % | Ubiquitina (valida el manejo de cadenas distintas de A). |
| `6LU7` | 70 % | Proteasa principal de SARS-CoV-2. |

> A identidad alta, proteínas rígidas (lisozima, CDK2…) producen animaciones casi
> estáticas: el movimiento visible depende de lo diferentes que sean las estructuras.

### Modo FASTA

Calmodulina (sube la identidad a ~60 % para centrar la búsqueda):

```
>calmodulina_1CLL
ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNG
TIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVD
EMIREADIDGDGQVNYEEFVQMMTAK
```

---

## Cómo funciona

```
Entrada (FASTA o código PDB)
        │
        ▼
src/pdb_fetcher.py   →  busca homólogos (BLAST / RCSB) y descarga los .pdb
        │
        ▼
src/interpolator.py  →  filtra cadena · alinea · ordena por RMSD ·
                        reduce a átomos comunes · interpola (PyMOL headless)
        │
        ▼
src/animator.py      →  renderiza cada fotograma y monta el GIF (Pillow)
        │
        ▼
   GIF animado  (data/output_animations/)
```

`app.py` orquesta las tres etapas y expone la interfaz Gradio.

---

## Estructura del proyecto

```
protein_animator/
├── app.py                     # Interfaz Gradio y orquestación del pipeline
├── requirements.txt
├── README.md
├── SantosPaula_informe.qmd    # Informe del proyecto (Quarto)
├── src/
│   ├── pdb_fetcher.py         # Búsqueda (BLAST/RCSB) y descarga de estructuras
│   ├── interpolator.py        # Alineamiento, ordenación e interpolación (PyMOL)
│   └── animator.py            # Render de fotogramas y montaje del GIF
└── data/
    ├── input_pdbs/            # PDB descargados (ignorado por git)
    └── output_animations/     # GIF generados (ignorado por git)
```

---

## Limitaciones

- La interpolación lineal es una aproximación geométrica simple; en transiciones muy
  grandes puede generar conformaciones intermedias no físicas.
- La animación no es una trayectoria energéticamente realista, sino una secuencia
  ordenada de estados experimentales.

---

Más detalles sobre el diseño, las decisiones y los resultados en el informe:
`SantosPaula_informe.qmd` (o el HTML renderizado).
