#!/bin/bash
source project_env.sh
# Ensure my_softwares is set
if [ -z "${my_softwares}" ]; then
  echo "my_softwares is not set. Exiting."
  exit 1
fi

# Create directory for custom installations
mkdir -p ${my_softwares}
cd ${my_softwares}

# Set environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"
export PKG_CONFIG_PATH=${my_softwares}/lib/pkgconfig


cd ${my_softwares}
BIN_VERSION="1.6.1"
mkdir -p ${my_softwares}/DeepVariant
cd ${my_softwares}/DeepVariant
ml apptainer
apptainer pull docker://google/deepvariant:"${BIN_VERSION}"

mkdir -p ${my_softwares}/GLNexus
cd ${my_softwares}/GLNexus
apptainer pull docker://quay.io/mlin/glnexus:v1.2.7
cd ${my_softwares}