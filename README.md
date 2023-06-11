# MagicLamp
## A software package for annotation of genomic datasets using discreet HMM sets.


## Citing MagicLamp
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
- HMMs_dir is the folder containing HMM files.
- You can also provide a meta-data file via the -m argument with gene and pathway names for each provided HMM (formatted after the _hmm_summary.csv_ file in this repo).

### Relying on pre-existing annotations
    YfGenie.py --gff -y genes.tsv -a GCF_023585845.1 -o GCF_023585845.1 -t 16

- In the above command _GCF_023585845.1_ represents the RefSeq assembly accession.
- genes.tsv is a single-column file listing gene names of interest (example file of the same names can be found in this repo).

### Simply extracting amino acid usgae frequencies and GC content from the provided genome assembly
    YfGenie.py --gc -a GCF_023585845.1 -o GCF_023585845.1 -t 16

### The works (YfGenie can be run in multiple modes at once)
    YfGenie.py --hmm --gff --gc -d HMMs_dir -y genes.tsv -a GCF_023585845.1 -o GCF_023585845.1 -t 16


### Usage
*(for a list of all available genies)*

    MagicLamp.py help

 *(for a detailed help menu for each genie)*
 
      MagicLamp.py FeGenie -h
      MagicLamp.py LithoGenie -h
      MagicLamp.py WspGenie -h
      MagicLamp.py GasGenie -h
      MagicLamp.py MagnetoGenie -h
      MagicLamp.py RosGenie -h
      MagicLamp.py Lucifer -h
      MagicLamp.py HmmGenie -h
 
 
### Quick-start
#### simplest form of the command
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output

#### if you are providing gene-calls add the --orfs flag
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --orfs
 
#### you can also provide GenBank files using the --gbk flag
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --gbk
 
#### if you want to normalize the number of identified genes to the total number of ORFs in the dataset, use the --norm flag
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --norm
  
#### with LithoGenie, you can also specficy a specific category for heatmap analysis using the -cat flag
    MagicLamp.py LithoGenie -bin_dir genomes/ -bin_ext fna -out lithogenie_output -cat sulfur

#### if you already ran LithoGenie once, but want to re-do the heatmap file with a different category, use the --skip flag
    MagicLamp.py LithoGenie -bin_dir genomes/ -bin_ext fna -out lithogenie_output --cat iron --skip

#### you can also include prediction of signal peptides and transmembrane domains with [Phobius](https://phobius.sbc.su.se/), using the --phobius flag
    MagicLamp.py LithoGenie -bin_dir genomes/ -bin_ext fna -out lithogenie_output --cat iron --skip

#### full list of LithoGenie categories (default action of the program is to create a broad-category heatmap)
* iron
* sulfur
* nitrogen
* oxygen
* methane
* hydrogen
* arsenic
* nitriles
* manganese
* carbon-monoxide
* halogenetated-compounds
* C1compounds
* urea
* selenium

*many of the HMMs used by LithoGenie were developed and compiled by Karthik Anantharaman: https://github.com/kanantharaman/metabolic-hmms*

*Phobius is not available through Anaconda, so the executables (phobius.pl, phobius.options, phobius.model, decodeanhmm, and decodeanhmm.64bit) are included in this repository. Standalone copy of Phobius was obtained from the following website: https://phobius.sbc.su.se/data.html, and users should cite the following article:*

Käll, L., Krogh, A., & Sonnhammer, E. L. L. (2004). A combined transmembrane topology and signal peptide prediction method. Journal of Molecular Biology, 338(5), 1027–1036.

