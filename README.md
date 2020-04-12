## A more targetted platform for annotation of datasets using HMM sets.

### quick install (if you have Conda installed)

    bash setup.sh
    source activate magiccave

#### Required Dependencies
* Python (version 3.6 or higher)
* BLAST (version 2.7.1+)
* Prodigal (version 2.6.3)
* HMMER (version 3.2.1)

#### Optional Dependencies
* Diamond (version 0.9.22.123) -- if you performing cross-validation against a reference database
* R (version 3.5.1) -- if generating plots
* Rscript -- if generating plots

### installation if you don't have Conda

Put MagicLamp.py script $PATH into your bash profile

    export PATH=$PATH:$(pwd)
                                      

### Usage
    MagicLamp.py help
*(this will provide you with a list of available genies)*

* MagicLamp.py FeGenie: HMM-based identification and categorization of iron genes and iron gene operons in genomes and metagenomes.
* MagicLamp.py LithoGenie: HMM-based identification and categorization of genes and operons relevant to chemolithoautotrophic metabolisms.
* MagicLamp.py RosGenie: HMM-based identifications of all genes responsible for neutralization of reactive-oxygen species.
* MagicLamp.py MagnetoGenie: HMM-based identification of genes responsible for magnetosome formation.
* MagicLamp.py WspGenie: HMM-based identification of the Wsp operon.
* MagicLamp.py Lucifer: HMM-based identification of light-sensing and light-producing genes
* MagicLamp.py GasGenie: HMM-based identification of genes responsible for gas vesicle formation.
* MagicLamp.py HmmGenie: Identification of a user-provided set of HMMs


      MagicLamp.py FeGenie -h
 *(this will provide you with a detailed help menu for each genie)*
 
 
### Quick-start
#### simplest form of the command
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output

#### if you are providing gene-calls
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --orfs
 
#### you can also provide GenBank files
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --gbk
 
#### if you want to normalize the number of identified genes to the total number of ORFs in the dataset
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --norm
  
#### with LithoGenie, you can also specficy a specific category for heatmap analysis
    MagicLamp.py LithoGenie -bin_dir genomes/ -bin_ext fna -out lithogenie_output --cat sulfur

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

*(many of the HMMs used by LithoGenie were developed and compiled by Karthik Anantharaman: https://github.com/kanantharaman/metabolic-hmms)*


