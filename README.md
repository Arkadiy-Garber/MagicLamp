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
 *(this will provide you with a detailed help menu)*
 
 
 ### Quick-start
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output
    



