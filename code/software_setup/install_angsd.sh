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

# Install zlib
wget http://www.zlib.net/zlib-1.3.1.tar.gz
tar -xzvf zlib-1.3.1.tar.gz
cd zlib-1.3.1
./configure --prefix=${my_softwares}
make
make install
cd ..

# Install bzip2 with -fPIC
wget https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz
tar -xzvf bzip2-1.0.8.tar.gz
cd bzip2-1.0.8
make CFLAGS="-fPIC"
make install PREFIX=${my_softwares}
cd ..

# Install xz
wget https://tukaani.org/xz/xz-5.2.5.tar.gz
tar -xzvf xz-5.2.5.tar.gz
cd xz-5.2.5
./configure --prefix=${my_softwares}
make
make install
cd ..

# Install OpenSSL (required for curl with TLS support)
wget https://www.openssl.org/source/openssl-1.1.1l.tar.gz
tar -xzvf openssl-1.1.1l.tar.gz
cd openssl-1.1.1l
./config --prefix=${my_softwares} --openssldir=${my_softwares}/ssl
make
make install
cd ..

# Install curl with OpenSSL
wget https://curl.se/download/curl-7.79.1.tar.gz
tar -xzvf curl-7.79.1.tar.gz
cd curl-7.79.1
./configure --prefix=${my_softwares} --with-zlib=${my_softwares} --with-ssl=${my_softwares}
make
make install
cd ..

# Set environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"
export PKG_CONFIG_PATH=${my_softwares}/lib/pkgconfig

# Clone and build HTSlib with submodules
git clone https://github.com/samtools/htslib.git
cd htslib
git submodule update --init --recursive
autoreconf -i
./configure --prefix=${my_softwares} --with-zlib=${my_softwares} --with-bzip2=${my_softwares} --with-lzma=${my_softwares} --with-curl=${my_softwares}
make
make install
cd ..

# Clone and compile ANGSD
git clone https://github.com/ANGSD/angsd.git
cd angsd

# Modify the Makefile for ANGSD to use local installations
sed -i "1iCXXFLAGS += -I${my_softwares}/include\nLDFLAGS += -L${my_softwares}/lib\nLDLIBS += -L${my_softwares}/lib -lbz2 -lcurl" Makefile
g++ -o angsd *.o${my_softwares}/htslib/libhts.a -L${my_softwares}//lib -lbz2 -lcurl -lz -lm -llzma -lpthread -lcrypto

# Compile ANGSD with explicit library paths
make HTSSRC=${my_softwares}/htslib/ CXXFLAGS="-I${my_softwares}/include" LDFLAGS="-L${my_softwares}/lib" LDLIBS="-L${my_softwares}/lib -lbz2 -lcurl"
echo 'export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc

cd ${my_softwares}

# Set environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"
export PKG_CONFIG_PATH=${my_softwares}/lib/pkgconfig

# Clone ngsRelate
git clone https://github.com/ANGSD/ngsRelate
cd ngsRelate

# Modify the Makefile for ngsRelate to use local installations and include libcrypto
sed -i "1iCXXFLAGS += -I${my_softwares}/include\nLDFLAGS += -L${my_softwares}/lib\nLDLIBS += -L${my_softwares}/lib -lbz2 -lcurl -lcrypto" Makefile

# Ensure correct library paths and compile ngsRelate
make HTSSRC=${my_softwares}/htslib/ CXXFLAGS="-I${my_softwares}/include" LDFLAGS="-L${my_softwares}/lib" LDLIBS="-L${my_softwares}/lib -lbz2 -lcurl -lcrypto"

# Compile ngsRelate with explicit library paths
g++ -O3 -o ngsRelate *.o ${my_softwares}/htslib/libhts.a -L${my_softwares}/lib -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -D__WITH_BCF__

# Verify the installation
./ngsRelate

cd ${my_softwares}
git clone --recursive https://github.com/clwgg/nQuire
cd nQuire
make submodules
make


cd ${my_softwares}
ml miniconda3
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd
conda env create -n "PCAngsd" -f environment.yml
conda activate PCAngsd
pip3 install -e .

cd ~
dir_save=$(pwd)/.local/lib/python3.12/site-packages/pcangsd 


echo "pcangsd=$dir_save"pcangsd > ~/.bashrc
source ~/.bashrc

cd ${my_softwares}
ml miniconda3
# Step 1: Create and activate a new conda environment
conda create -n ploidyNGS_env python=3.8 -y
conda activate ploidyNGS_env

# Step 2: Install required dependencies
conda install -c bioconda samtools -y
conda install -c conda-forge numpy pandas scipy matplotlib seaborn -y

# Step 3: Clone the ploidyNGS repository
git clone https://github.com/diriano/ploidyNGS.git
cd ploidyNGS

# Step 4: Install ploidyNGS
pip install .

# Install KING
cd ${my_softwares}
mkdir king
cd king
wget https://www.kingrelatedness.com/KINGcode.tar.gz
tar -xzvf KINGcode.tar.gz
c++ -lm -lz -O2 -fopenmp -o king *.cpp


# ldne
cd ${my_softwares}
mkdir ldne
cd ldne
cp ../ldne.zip .
unzip ldne.zip
cd LDNe/
unzip LDNe.zip
chmod +x LDNe


cd ${my_softwares}
mkdir plink1.9
cd plink1.9
cp ../plink_linux_x86_64_20240804.zip .
unzip plink_linux_x86_64_20240804.zip
chmod +x plink


plink_linux_x86_64_20240804.zip

