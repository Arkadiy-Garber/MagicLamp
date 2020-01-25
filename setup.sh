
conda create -n magiccave hmmer diamond prodigal blast --yes

source activate magiccave

mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d


## adding paths:
echo '#!/bin/sh'" \

export PATH=\"$(pwd):"'$PATH'\"" \

export gas_hmms=\"$(pwd)/hmms/gas\"

export wsp_hmms=\"$(pwd)/hmms/wsp\"

export ros_hmms=\"$(pwd)/hmms/ros\"

export lux_hmms=\"$(pwd)/hmms/lux\"

export magneto_hmms=\"$(pwd)/hmms/magneto\"

export litho_hmms=\"$(pwd)/hmms/litho\"

export iron_hmms=\"$(pwd)/hmms/iron\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
