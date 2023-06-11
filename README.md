# MagicLamp
### A software package for annotation of genomic datasets using discreet HMM sets.


### Citing MagicLamp
There is no official publication for MagicLamp. If it was useful for your work, please cite as follows:

Garber, AI., Ramirez, GA., Merino, N., Pavia MJ., McAllister, SM. (2020) MagicLamp: toolkit for annotation of genomic data using discreet and curated HMM sets. 2023: MagicLamp, GitHub repository: https://github.com/Arkadiy-Garber/MagicLamp.

## Installation (Conda is required for this software)
    git clone https://github.com/Arkadiy-Garber/MagicLamp.git
    cd MagicLamp
    bash setup.sh
    conda activate magiclamp

### Using your own HMM set with YfGenie
    YfGenie.py --hmm -d HMMs_dir -a GCF_023585845.1 -o GCF_023585845.1 -t 16

- In the above command _GCF_023585845.1_ represents the RefSeq assembly accession.
- HMMs_dir is the folder containing raw HMM files. See the subfolders inside the hmms directory to see what these look like.
- You can also provide a meta-data file via the -m argument with gene and pathway names for each provided HMM (formatted after the _hmm_summary.csv_ file in this repo).

### Relying on pre-existing annotations rom the GFF file (recommended only if the gene/pathway is well-annotated)
    YfGenie.py --gff -y genes.tsv -a GCF_023585845.1 -o GCF_023585845.1 -t 16

- genes.tsv is a single-column file listing gene names of interest (example file of the same names can be found in this repo).

### Simply extracting amino acid usgae frequencies and GC content from the provided genome assembly
    YfGenie.py --gc -a GCF_023585845.1 -o GCF_023585845.1 -t 16

- this will generate a single-line TSV file that lists usage frequences for each amino acid residue.

### The works (YfGenie can be run in multiple modes at once)
    YfGenie.py --hmm --gff --gc -d HMMs_dir -y genes.tsv -a GCF_023585845.1 -o GCF_023585845.1 -t 16


<sup>*Phobius is not available through Anaconda, so the executables (phobius.pl, phobius.options, phobius.model, decodeanhmm, and decodeanhmm.64bit) are included in this repository. Standalone copy of Phobius was obtained from the following website: https://phobius.sbc.su.se/data.html, and users should cite appropriately.*</sup>

