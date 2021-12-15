## MagicLamp
### A software package for annotation of genomic datasets using discreet HMM sets.


#### Citing MagicLamp
There is no official publication for MagicLamp. If it was useful for your work, please cite as follows:

Garber, AI., Ramirez, GA., Merino, N., Pavia MJ., McAllister, SM. (2020) MagicLamp: toolkit for annotation of 'omics datasets using curated HMM sets. 2021: MagicLamp, GitHub repository: https://github.com/Arkadiy-Garber/MagicLamp.


### Submitting custom HMM sets
We encourage users to submit their custom-built HMM sets for specific pathways or processes. Please email the following material to agarber4@asu.edu: 1) zipped folder containing the HMM files, 2) filled-out rules-template.csv file. Optionally, users can also provide 1) genome(s) known to encode the genetic pathway(s), and 2) genome(s) known not to encode the genetic pathway.


### Installation
#### quick install (if you have Conda installed)
    git clone https://github.com/Arkadiy-Garber/MagicLamp.git
    cd MagicLamp
    bash setup.sh
    source activate magiclamp
(if "source activate magiclamp" does not work, you can use "conda activate magiclamp")


#### installation if you don't have Conda

Put MagicLamp.py script $PATH into your bash profile

    export PATH=$PATH:$(pwd)

#### Required Dependencies
* Python (version 3.6 or higher)
* BLAST (version 2.7.1+)
* Prodigal (version 2.6.3)
* HMMER (version 3.2.1)

#### Optional Dependencies
* Diamond (version 0.9.22.123) -- if you performing cross-validation against a reference database
* R (version 3.5.1) -- if generating plots
* Rscript -- if generating plots
                                      

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



## Tutorial (Binder)

If you would like to learn more about what is going on under the hood of this program, and see examples of how to use it, please see the following slideshow and associated video:

[Slideshow](https://github.com/biovcnet/topic-functional-annotation/blob/master/Lesson-5/Lesson-5-Annotation-with-HMMs.pdf) | [Video presentation](https://www.youtube.com/watch?v=TUjxP_8d-MY)


## Walkthrough

If you'd like to try out this program in a virtual terminal session with all dependencies and test dataset pre-loaded, please follow the tutorials below.

These videos were made as part of the [Bioinformatics Virtual Coordination Network](https://biovcnet.github.io/) :)


### Tutorial 1: Using pre-existing HMM sets included with MagicLamp (e.g. lithotrophy, respiration, iron, ROS, etc.)

[![Binder](https://mybinder.org/badge_logo.svg)](https://gesis.mybinder.org/binder/v2/gh/biovcnet/bvcn-binder-magiclamp/master?urlpath=lab) <-- click here to begin

(Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the awesome [binder](https://mybinder.org/) team!)

You can also follow along in this linked video: 
[Video presentation 1](https://www.youtube.com/watch?v=zhvHgVkURwg)


Enter the MagicLamp repository

    cd MagicCave/

print the MagicLamp help menu

    MagicLamp.py help

print WspGenie help menu

    MagicLamp.py WspGenie -h

run WspGenie on test dataset

    MagicLamp.py WspGenie -bin_dir test_dataset/ -bin_ext fna -out wspgenie_out


go into the wspgenie output directory and check out the output file

    cd wspgenie_out/
    less -S wspgenie-summary.csv

check out the gene predictions

    cd ORF_calls/
    cd ../../

mv ORF calls to the main directory

    mv wspgenie_out/ORF_calls/ ./

print LithoGenie help menu

    MagicLamp.py LithoGenie -h

run LithoGenie on ORF calls

    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out

check out the output

    cd lithogenie_out/
    less -S lithogenie-summary.csv
    less lithogenie.ALL.heatmap.csv
    cd ../

re-run LithoGenie to create a .heatmap.csv for an element-of-interest

    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out --skip -cat sulfur
    # answer 'y' to the question
    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out --skip -cat iron

check out the new results

    cd lithogenie_out/
    less lithogenie.sulfur.heatmap.csv
    less lithogenie.iron.heatmap.csv

print the HmmGenie help menu

    MagicLamp.py HmmGenie -h

run HmmGenie with a set of HMMs for gas vesicle formation

    MagicLamp.py HmmGenie -hmm_dir MagicCave/hmms/gas/ -hmm_ext hmm -bin_dir test_dataset/ -bin_ext fna -out gas_out

check out the results and re-run HmmGenie with more stringent parameters

    MagicLamp.py HmmGenie -hmm_dir MagicCave/hmms/gas/ -hmm_ext hmm -bin_dir test_dataset/ -bin_ext fna -out gas_out -clu 5

check out the results

    cd gas_out/
    less -S genie-summary.csv



### Tutorial 2: Making your own HMMs to use with HmmGenie

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/biovcnet/bvcn-binder-hmmgenie/master?urlpath=lab) <-- click here to begin

You can also follow along in this linked video: 
[Video presentation 2](https://www.youtube.com/watch?v=LKX678JzRfU)

Run phmmer on the Cyc1 fasta file

    phmmer -A Cyc1.refseq.msa --tblout Cyc1.refseq.tblout -E 1E-20 Cyc1.faa ../refseq_db/refseq_nr.sample.faa

Build HMM file from MSA (multiple sequence alignment) file, using hmmbuild

    hmmbuild Cyc1.hmm Cyc1.refseq.msa

Query the Cyc1 HMM file against refseq database sample

    hmmsearch --tblout Cyc1.hmm.refseq.tblout Cyc1.hmm ../refseq_db/refseq_nr.sample.faa

Examine the output file. What do the bit scores look like for likely false positives

    less Cyc1.hmm.refseq.tblout

Move into directory containing MtrA FASTA file, and create an alignment using Muscle.

    muscle -in MtrA.faa -out MtrA.fa

Build HMM file from MSA (multiple sequence alignment) file, using hmmbuild

    hmmbuild MtrA.hmm MtrA.fa

Query the MtrA HMM file against refseq database sample

    hmmsearch --tblout MtrA.hmm.nr.tblout MtrA.hmm ../refseq_db/refseq_nr.sample.faa

Examine the output file. What do the bit scores look like for likely false positives

    less MtrA.hmm.nr.tblout

Move the HMM files into a single directory

    mv MtrA.hmm ../HMMs/
    mv Cyc1.hmm ../HMMs/

Check out the Pfam-derived HMM and bitscores.txt file

    less Catalase.hmm

Run HmmGenie (MagicLamp) on test dataset using the new HMM collection

    MagicLamp.py HmmGenie -hmm_dir HMMs/ -hmm_ext hmm -bin_dir test_data/ -bin_ext txt -out hmmgenie_out -eval 1E-1
    MagicLamp.py HmmGenie -hmm_dir HMMs/ -hmm_ext hmm -bin_dir test_data/ -bin_ext txt -out hmmgenie_out -bit HMMs/bitscores.txt
    MagicLamp.py HmmGenie -hmm_dir HMMs/ -hmm_ext hmm -bin_dir test_data/ -bin_ext txt -out hmmgenie_out -bit HMMs/bitscores.txt -clu 2

