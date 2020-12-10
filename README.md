### A software package for annotation of genomic datasets using discreet HMM sets.

#### Citing MagicLamp
There is no official publication for MagicLamp. If it was useful for your work, you can cite it as follows:

Garber, AI and Miller, SR (2020) MagicLamp: toolkit for annotation of 'omics datasets using curated HMM sets. 2020: ParaHunter, GitHub repository: https://github.com/Arkadiy-Garber/MagicLamp.


### installation
#### quick install (if you have Conda installed)
    git clone 
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

#### installation if you don't have Conda

Put MagicLamp.py script $PATH into your bash profile

    export PATH=$PATH:$(pwd)
                                      

### Usage
*(for a list of available genies)*

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

#### if you are providing gene-calls
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --orfs
 
#### you can also provide GenBank files
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --gbk
 
#### if you want to normalize the number of identified genes to the total number of ORFs in the dataset
    MagicLamp.py FeGenie -bin_dir genomes/ -bin_ext fna -out fegenie_output --norm
  
#### with LithoGenie, you can also specficy a specific category for heatmap analysis
    MagicLamp.py LithoGenie -bin_dir genomes/ -bin_ext fna -out lithogenie_output --cat sulfur

#### if you already ran LithoGenie once, but want to re-do the heatmap file with a different category
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

*(many of the HMMs used by LithoGenie were developed and compiled by Karthik Anantharaman: https://github.com/kanantharaman/metabolic-hmms)*


