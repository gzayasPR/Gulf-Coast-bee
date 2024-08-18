#!/bin/bash
source project_env.sh

cd $my_code/software_setup/

cp ../project_env.sh .

echo "Starting Dockers"
bash install_dockers.sh
echo "Finished"

echo "Starting R"
bash install_R_packages.sh
echo "Finished"

echo "Starting Others"
bash install_softwares.sh
echo "Finished"

echo "Starting Conda"
bash install_conda_env.sh
echo "Finished"