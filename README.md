# MagicLamp
### A software package for annotation of genomic datasets using discreet HMM sets.


## Citing MagicLamp
There is no official publication for MagicLamp. If it was useful for your work, please cite as follows:

Garber, AI., Ramirez, GA., Merino, N., Pavia MJ., McAllister, SM. (2020) MagicLamp: toolkit for annotation of 'omics datasets using curated HMM sets. 2021: MagicLamp, GitHub repository: https://github.com/Arkadiy-Garber/MagicLamp.

## Installation
### Quick install (if you have Conda installed)
    git clone https://github.com/Arkadiy-Garber/MagicLamp.git
    cd MagicLamp
    bash setup.sh
    source activate magiclamp
(if "source activate magiclamp" does not work, you can use "conda activate magiclamp")


### Installation if you don't have Conda

Put MagicLamp.py script $PATH into your bash profile

    export PATH=$PATH:$(pwd)

#### Dependencies (only if you did not use the setup.sh script to configure a Conda environment)
##### Required
* Python (version 3.6 or higher)
* BLAST (version 2.7.1+)
* Prodigal (version 2.6.3)
* HMMER (version 3.2.1)

##### Optional
* Diamond (version 0.9.22.123) -- if you performing cross-validation against a reference database
* R (version 3.5.1) -- if generating plots
* Rscript -- if generating plots


### Using HmmGenie
An improtant feature of MagicLamp is the ability of users to use their own HMM sets with HmmGenie, like so:

    MagicLamp.py HmmGenie -bin_dir bins/ -bin_ext fna -hmm_dir HMMs/ -hmm_ext hmm -rules rules-template.csv

In the above example, the HMMs/ directory contains the user-compiled or created HMMs, where each HMM file ends with a .hmm filename extension. The rules-template.csv file is provided in the main MagicLamp repository (should be in the same folder as this readme). Users are to fill out this file, providing information on each HMM that they wish to use with this program.

### Submitting custom HMM sets
We encourage users to submit their custom-built or compiled HMM sets for specific pathways or processes. Please email the following material to agarber4@asu.edu: 1) zipped folder containing the HMM files, 2) filled-out rules-template.csv file. Optionally, users can also provide 1) genome(s) known to encode the genetic pathway(s), and 2) genome(s) known not to encode the genetic pathway.
                               

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

