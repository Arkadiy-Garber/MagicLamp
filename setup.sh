#!/usr/bin/env bash

printf "\n    ${GREEN}Setting up conda environment...${NC}\n\n"

## adding conda channels
conda config --add channels defaults 2> /dev/null
conda config --add channels bioconda 2> /dev/null
conda config --add channels conda-forge 2> /dev/null
conda config --add channels au-eoed 2> /dev/null

conda create -n magiclamp -c r r-dplyr, r-tibble, r-stringr r-fuzzyjoin r-RColorBrewer r-forcats r-plotly r-ggplot2 r-stringi r-reshape r-reshape2 r-tidyverse r-argparse r-ggdendro r-pvclust python=3.7 hmmer diamond prodigal blast metabat2 --yes

#Rscript -e 'install.packages("grid", repos = "http://cran.us.r-project.org")'
Rscript -e 'install.packages("broom", repos = "http://cran.us.r-project.org”)'
Rscript -e 'install.packages("ggpubr", repos = "http://cran.us.r-project.org”)'

## activating environment
conda activate magiclamp

## creating directory for conda-env-specific source files
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d

## adding paths:
echo '#!/bin/sh'" \

export rscripts=\"$(pwd)/rscripts\"

export PATH=\"$(pwd):"'$PATH'\"" \

export gas_hmms=\"$(pwd)/hmms/gas\"

export wsp_hmms=\"$(pwd)/hmms/wsp\"

export ros_hmms=\"$(pwd)/hmms/ros\"

export lux_hmms=\"$(pwd)/hmms/lux\"

export magneto_hmms=\"$(pwd)/hmms/magneto\"

export litho_hmms=\"$(pwd)/hmms/litho\"

export mn_hmms=\"$(pwd)/hmms/manganese\"

export circ_hmms=\"$(pwd)/hmms/circ\"

export iron_hmms=\"$(pwd)/hmms/iron\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

# re-activating environment so variable and PATH changes take effect
conda activate magiclamp


printf "\n        ${GREEN}DONE!${NC}\n\n"
