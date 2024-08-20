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
if [ -d "zlib-1.3.1" ]; then
  rm -rf zlib-1.3.1
fi
wget http://www.zlib.net/zlib-1.3.1.tar.gz
tar -xzvf zlib-1.3.1.tar.gz
cd zlib-1.3.1
./configure --prefix=${my_softwares}
make && make install || { echo "zlib installation failed"; exit 1; }
cd ..

# Install bzip2 with -fPIC
if [ -d "bzip2-1.0.8" ]; then
  rm -rf bzip2-1.0.8
fi
wget https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz
tar -xzvf bzip2-1.0.8.tar.gz
cd bzip2-1.0.8
make CFLAGS="-fPIC"
make install PREFIX=${my_softwares} || { echo "bzip2 installation failed"; exit 1; }
cd ..

# Install xz
if [ -d "xz-5.2.5" ]; then
  rm -rf xz-5.2.5
fi
wget https://tukaani.org/xz/xz-5.2.5.tar.gz
tar -xzvf xz-5.2.5.tar.gz
cd xz-5.2.5
./configure --prefix=${my_softwares}
make && make install || { echo "xz installation failed"; exit 1; }
cd ..

# Install OpenSSL (required for curl with TLS support)
if [ -d "openssl-1.1.1l" ]; then
  rm -rf openssl-1.1.1l
fi
wget https://www.openssl.org/source/openssl-1.1.1l.tar.gz
tar -xzvf openssl-1.1.1l.tar.gz
cd openssl-1.1.1l
./config --prefix=${my_softwares} --openssldir=${my_softwares}/ssl
make && make install || { echo "OpenSSL installation failed"; exit 1; }
cd ..

# Install curl with OpenSSL
if [ -d "curl-7.79.1" ]; then
  rm -rf curl-7.79.1
fi
wget https://curl.se/download/curl-7.79.1.tar.gz
tar -xzvf curl-7.79.1.tar.gz
cd curl-7.79.1
./configure --prefix=${my_softwares} --with-zlib=${my_softwares} --with-ssl=${my_softwares}
make && make install || { echo "curl installation failed"; exit 1; }
cd ..

# Install htslib
if [ -d "htslib" ]; then
  rm -rf htslib
fi
git clone https://github.com/samtools/htslib.git
cd htslib
git submodule update --init --recursive
autoreconf -i
./configure --prefix=${my_softwares} \
            --with-zlib=${my_softwares} \
            --with-bzip2=${my_softwares} \
            --with-lzma=${my_softwares} \
            --with-curl=${my_softwares}
make && make install || { echo "htslib installation failed"; exit 1; }
cd ..

wget https://raw.githubusercontent.com/ANGSD/angsd/master/misc/ngsadmix32.cpp
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
cd ..

# Clone ngsRelate
cd ${my_softwares}
git clone https://github.com/ANGSD/ngsRelate
cd ngsRelate

# Modify the Makefile for ngsRelate to use local installations and include libcrypto
sed -i "1iCXXFLAGS += -I${my_softwares}/include\nLDFLAGS += -L${my_softwares}/lib\nLDLIBS += -L${my_softwares}/lib -lbz2 -lcurl -lcrypto" Makefile

# Ensure correct library paths and compile ngsRelate
make HTSSRC=${my_softwares}/htslib/ CXXFLAGS="-I${my_softwares}/include" LDFLAGS="-L${my_softwares}/lib" LDLIBS="-L${my_softwares}/lib -lbz2 -lcurl -lcrypto"

# Compile ngsRelate with explicit library paths
g++ -O3 -o ngsRelate *.o ${my_softwares}/htslib/libhts.a -L${my_softwares}/lib -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -D__WITH_BCF__

