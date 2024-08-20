#!/bin/bash
source project_env.sh

# Ensure my_softwares is set
if [ -z "${my_softwares}" ]; then
  echo "my_softwares is not set. Exiting."
  exit 1
fi

# Create directory for custom installations
mkdir -p "${my_softwares}"
cd "${my_softwares}"

# Set environment variables
export PATH="${my_softwares}/bin:$PATH"
export LD_LIBRARY_PATH="${my_softwares}/lib:$LD_LIBRARY_PATH"
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"
export PKG_CONFIG_PATH="${my_softwares}/lib/pkgconfig"

# Load Miniconda module
ml miniconda3

# Initialize Conda if necessary
conda init bash
source ~/.bashrc

# Create Conda environment in the same directory
if [ -d "${my_softwares}/pcangsd/env" ]; then
  echo "PCAngsd environment already exists."
  conda activate "${my_softwares}/pcangsd/env"
else
 # Clone and install PCAngsd
 git clone https://github.com/Rosemeis/pcangsd.git
 cd pcangsd
  conda env create --prefix "${my_softwares}/pcangsd/env" -f environment.yml
  conda activate "${my_softwares}/pcangsd/env"
  pip3 install -e .
fi

# Save PCAngsd directory to .bashrc
dir_save=$(pwd)
echo "export pcangsd=${dir_save}" >> ~/.bashrc
source ~/.bashrc

cd "${my_softwares}"

# Install Trim Galore in a new Conda environment
if [ -d "${my_softwares}/trim_galore/env" ]; then
  echo "Trim Galore environment already exists."
  conda activate "${my_softwares}/trim_galore/env"
else
  conda create --prefix "${my_softwares}/trim_galore/env" -c bioconda trim-galore -y
  conda activate "${my_softwares}/trim_galore/env"
fi

# Save Trim Galore directory to .bashrc
dir_save=$(pwd)
echo "export trim_galore=${dir_save}" >> ~/.bashrc
source ~/.bashrc

cd "${my_softwares}"

# Install Qualimap in a new Conda environment
if [ -d "${my_softwares}/qualimap/env" ]; then
  echo "Qualimap environment already exists."
  conda activate "${my_softwares}/qualimap/env"
else
  conda create --prefix "${my_softwares}/qualimap/env" -c bioconda qualimap -y
  conda activate "${my_softwares}/qualimap/env"
fi

# Install Qualimap in a new Conda environment
if [ -d "${my_softwares}/angsd/env" ]; then
  echo "Qualimap environment already exists."
  conda activate "${my_softwares}/angsd/env"
else
  conda create --prefix "${my_softwares}/angsd/env" -c bioconda angsd -y
  conda activate "${my_softwares}/angsd/env"
fi

# Install multiqc in a new Conda environment
if [ -d "${my_softwares}/multiqc/env" ]; then
  echo "Qualimap environment already exists."
  conda activate "${my_softwares}/multiqc/env"
else
  conda create --prefix "${my_softwares}/multiqc/env" -c bioconda angsd -y
  conda activate "${my_softwares}/multiqc/env"
fi
