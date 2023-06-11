# MagicLamp - YfGenie
### A software package for annotation of genomic datasets using discreet HMM sets.


### Citing MagicLamp
There is no official publication for MagicLamp. If it was useful for your work, please cite as follows:

Garber, AI., Ramirez, GA., Merino, N., Pavia MJ., McAllister, SM. (2020) MagicLamp: toolkit for annotation of genomic data using discreet and curated HMM sets. 2023: MagicLamp, GitHub repository: https://github.com/Arkadiy-Garber/MagicLamp.

## Installation (Conda is required for this software)
    git clone https://github.com/Arkadiy-Garber/MagicLamp.git
    cd MagicLamp
    bash setup.sh
    conda activate magiclamp

## Usage

### Using your own HMM set with YfGenie
    YfGenie.py --hmm -d HMMs_dir -a GCF_023585845.1 -o GCF_023585845.1 -t 16

_- In the above command _GCF_023585845.1_ represents the RefSeq assembly accession._

_- HMMs_dir is the folder containing raw HMM files. See the subfolders inside the hmms directory to see what these look like._

_- You can also provide a meta-data file via the -m argument with gene and pathway names for each provided HMM (formatted after the _hmm_summary.csv_ file in this repo)._

### Relying on annotation within GFF file (recommended if the gene/pathway is well-annotated)
    YfGenie.py --gff -y genes.tsv -a GCF_023585845.1 -o GCF_023585845.1 -t 16

_- genes.tsv is a single-column file listing gene names of interest (example file of the same names can be found in this repo)._

### Simply extracting amino acid usage frequencies and GC content from the provided genome assembly
    YfGenie.py --gc -a GCF_023585845.1 -o GCF_023585845.1 -t 16

_- this will generate a single-line TSV file that lists usage frequences for each amino acid residue._

### The works (YfGenie can be run in multiple modes at once)
    YfGenie.py --hmm --gff --gc -d HMMs_dir -y genes.tsv -a GCF_023585845.1 -o GCF_023585845.1 -t 16

### You can also provide locally annotated files via the -c CONTIGS, -g GFF, and -p PROTs arguments
    YfGenie.py --hmm --gff --gc -d HMMs_dir -y genes.tsv -c genome.fa -g genome.gff -p genome.faa -o genome_out -t 16




