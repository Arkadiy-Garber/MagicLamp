# MagicLamp
========
### A software package for annotation of genomic datasets using discreet HMM sets.

### Citing MagicLamp
Please cite our [pre-print](https://www.biorxiv.org/content/10.64898/2026.07.19.739457v1.full.pdf)

Special thanks to AstrobioMike and his [bit](https://github.com/AstrobioMike/bit) software package for enabling easy access to NCBI's RefSeq and GenBank assemblies.

## Installation (Conda is required for this software)

    git clone https://github.com/Arkadiy-Garber/MagicLamp.git
    cd MagicLamp
    bash setup.sh
    conda activate magiclamp

`setup.sh` builds the `magiclamp` conda environment in two tiers:

1. **Full install (attempted first):** the core tools (`hmmer`, `blast`, `prodigal`, `diamond`, `metabat2`, `bit`) **plus** the R stack used by the optional `--makeplots` figures.
2. **Minimal install (automatic fallback):** if the full solve fails — most often because of an R-package conflict — `setup.sh` automatically retries with just the core tools. In minimal mode **every genie runs normally**; only the optional `--makeplots` figures are unavailable.

After the environment is created, `setup.sh` sources **`magiclamp.paths`**, which writes all of the per-genie HMM directory variables (`iron_hmms`, `litho_hmms`, `atp_hmms`, `abx_hmms`, `motility_hmms`, …) and the `rscripts` path into the conda environment's `activate.d` hook. `magiclamp.paths` is the single source of truth for these variables — if you add a genie, add its variable there. It also warns (without failing) if any expected `hmms/` subdirectory is missing.

## Usage

Run any genie through the `MagicLamp.py` dispatcher:

    MagicLamp.py <Genie> -bin_dir <input_directory> -bin_ext <extension> -out <output_directory> -t <threads>

List all available genies and a one-line description of each:

    MagicLamp.py help

### Available genies

| Genie | Target | HMM variable |
|-------|--------|--------------|
| **FeGenie** | Iron acquisition, oxidation, reduction, and storage genes and operons | `iron_hmms` |
| **LithoGenie** | Chemolithoautotrophic metabolisms (S, H₂, N, methane, and more) | `litho_hmms` |
| **Lucifer** | Light-sensing and light-producing (bioluminescence/rhodopsin) genes | `lux_hmms` |
| **ATPGenie** | ATP synthases (F-, V-, and A-type) | `atp_hmms` |
| **PortGenie** | Sodium antiporters and symporters | `portna_hmms` |
| **RnfGenie** | Rnf complex (`rnfABCDGE`) | `rnf_hmms` |
| **AbxGenie** | Antibiotic-biosynthesis genes (subcategorized) | `abx_hmms` |
| **MotiliGenie** | Motility, chemotaxis, thermotaxis, and pili genes | `motility_hmms` |
| **RiboGenie** | Ribosomal-protein and translation-machinery genes | `ribo_hmms` |
| **ResistiGenie** | Resistance genes (heavy metals, UV, ROS, radioactivity, antibiotics) | `resist_hmms` |
| **SporeGenie** | Sporulation genes, subcategorized by sporulation stage | `spore_hmms` |
| **OmniGenie** | Generic runner for any genie's HMM set (see below) | `<genie>_hmms` |
| **HmmGenie** | Annotation with a user-provided set of HMMs | — |

### Basic usage

Point a genie at a directory of genome assemblies (nucleotide FASTA), giving the file extension so it knows which files to process:

    MagicLamp.py ATPGenie -bin_dir genomes/ -bin_ext fna -out atpgenie_out -t 16

- `-bin_dir` — directory containing your genomes/bins.
- `-bin_ext` — file extension of those files, **without** the leading period (e.g. `fna`, `fa`, `fasta`).
- `-out` — output directory (each genie has its own default, e.g. `atpgenie_out`, `ribogenie_out`).
- `-t` — number of threads for HMMER (and DIAMOND, where applicable).

By default, MagicLamp calls ORFs from nucleotide contigs (via Prodigal) before running HMMER.

### Common input-format flags

    # GenBank-format bins
    MagicLamp.py RiboGenie -bin_dir genomes/ -bin_ext gbk --gbk -out ribogenie_out -t 16

    # Metagenomic / metatranscriptomic assemblies
    MagicLamp.py LithoGenie -bin_dir mags/ -bin_ext fa --meta -out litho_out -t 16

### Plots, normalization, and other options

- `--makeplots` — generate heatmaps/figures via R (**requires the full install**; unavailable in minimal mode).
- `--norm` — normalize gene counts to the number of predicted ORFs per genome instead of raw counts.
- `--all_results` — report all hits regardless of clustering/operon structure (where supported).
- `-ref <db.faa>` — reference protein database for a DIAMOND-based confirmation step.
- `-d <int>` — maximum distance (in genes) between genes to be considered part of the same cluster/operon.
- `-inflation <int>` — inflation factor for final gene-category counts (default 1000).

### LithoGenie: highlighting a metabolism in the heatmap

    MagicLamp.py LithoGenie -bin_dir genomes/ -bin_ext fna -cat sulfur --makeplots -out litho_out -t 16

Valid `-cat` values: `sulfur, hydrogen, methane, nitrogen, oxygen, carbon-monoxide, C1compounds, carbon, urea, halogenetated-compounds, arsenic, selenium, nitriles, iron, ALL` (default: `ALL`).

### OmniGenie: run any genie's HMM set generically

`OmniGenie` is a genie-agnostic runner. Pass `-genie <name>`, where `<name>` is the **prefix of that genie's HMM environment variable** (i.e. the part before `_hmms`). OmniGenie resolves the library from `${<name>_hmms}`, falling back to `hmms/<name>/`.

    MagicLamp.py OmniGenie -genie atp    -bin_dir genomes/ -bin_ext fna -out genie_out -t 16
    MagicLamp.py OmniGenie -genie abx    -bin_dir genomes/ -bin_ext fna -out genie_out -t 16
    MagicLamp.py OmniGenie -genie portna -bin_dir genomes/ -bin_ext fna -out genie_out -t 16

Valid `-genie` values match the HMM-variable prefixes in the table above: `iron, litho, lux, atp, portna, rnf, abx, motility, ribo, resist, spore`.

### Using your own HMM set (HmmGenie)

To annotate with a custom set of HMMs, place the raw `.hmm` files in a directory and point HmmGenie at it:

    MagicLamp.py HmmGenie -bin_dir genomes/ -bin_ext fna -hmm_dir my_hmms/ -out hmmgenie_out -t 16
