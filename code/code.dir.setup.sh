#!/bin/bash
# Source the environment variables
source project_env.sh

cd $my_code/Population_Genomics

bash 0.Setup.sh


cd $my_code
cd $my_code/Variant_Calling/

bash 0.Setup.sh



cd $my_code
cd $my_code/ANGSD_code/
bash 0.Setup.sh