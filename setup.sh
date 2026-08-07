#!/usr/bin/env bash
#
# MagicLamp v3 setup.
#
# Tries to build the FULL environment (core bioinformatics tools + the R stack
# used by --makeplots). If the full solve fails (R packages are the usual
# culprit), it falls back to a MINIMAL environment with just the core tools:
#   hmmer, blast, prodigal, diamond, bit
# In minimal mode all genies run normally; only the optional --makeplots
# figures are unavailable.

set -uo pipefail

GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m'

ENV_NAME="magiclamp"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Core tools every genie needs (no R, no plotting).
CORE_PKGS=(hmmer blast prodigal diamond metabat2 bit python=3.7)

# R stack used only for --makeplots.
R_PKGS=(r-base r-dplyr r-tibble r-stringr r-fuzzyjoin r-RColorBrewer r-forcats
        r-plotly r-ggplot2 r-stringi r-reshape r-reshape2 r-tidyverse r-argparse
        r-ggdendro r-pvclust r-ggpubr r-broom)

CHANNELS=(-c conda-forge -c bioconda -c defaults -c astrobiomike -c r)

printf "\n    ${GREEN}Setting up the MagicLamp conda environment...${NC}\n\n"

## make sure conda channels are configured
conda config --add channels defaults    2> /dev/null || true
conda config --add channels bioconda     2> /dev/null || true
conda config --add channels conda-forge  2> /dev/null || true

## remove any stale environment of the same name so retries are clean
conda env remove -n "${ENV_NAME}" -y 2> /dev/null || true

FULL_INSTALL_OK=0

## ---- Attempt 1: FULL environment (core tools + R plotting stack) ----
printf "    ${GREEN}Attempting full install (core tools + R plotting stack)...${NC}\n"
if conda create -n "${ENV_NAME}" "${CHANNELS[@]}" \
        "${CORE_PKGS[@]}" "${R_PKGS[@]}" --yes; then
    FULL_INSTALL_OK=1
    printf "\n    ${GREEN}Full install succeeded (plotting via --makeplots enabled).${NC}\n"
else
    printf "\n    ${YELLOW}Full install failed (usually an R-package solve).${NC}\n"
    printf "    ${YELLOW}Falling back to a minimal install with core tools only...${NC}\n\n"

    ## ---- Attempt 2: MINIMAL environment (core tools only) ----
    conda env remove -n "${ENV_NAME}" -y 2> /dev/null || true
    if conda create -n "${ENV_NAME}" -c conda-forge -c bioconda -c defaults -c astrobiomike \
            "${CORE_PKGS[@]}" --yes; then
        printf "\n    ${YELLOW}Minimal install succeeded. All genies will run, but the optional\n"
        printf "    --makeplots figures are NOT available (no R stack).${NC}\n"
    else
        printf "\n    ${RED}Minimal install also failed. Please check your conda setup and channels.${NC}\n" >&2
        exit 1
    fi
fi

## activate the environment
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

## export all HMM library paths + rscripts + PATH into the env's activate hook.
## magiclamp.paths is the single source of truth for these variables.
source "${SCRIPT_DIR}/magiclamp.paths"

## re-activate so the exported variables and PATH take effect
conda deactivate
conda activate "${ENV_NAME}"

if [ "${FULL_INSTALL_OK}" -eq 1 ]; then
    printf "\n        ${GREEN}DONE! (full install)${NC}\n\n"
else
    printf "\n        ${GREEN}DONE! (minimal install — no --makeplots)${NC}\n\n"
fi
