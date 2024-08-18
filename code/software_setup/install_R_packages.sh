#!/bin/bash

source project_env.sh

mkdir -p ${my_softwares}/R_Packages/
cd ${my_softwares}/R_Packages/
mkdir -p ${my_softwares}/R_Packages/4.4
mkdir -p ${my_softwares}/R_Packages/4.3
echo "R_LIBS_USER="${my_softwares}/R_Packages/%v""  > ${my_softwares}/R_Packages/.Renviron
ml r/4.4.0
Rscript ${my_code}/software_setup/R.setup.R