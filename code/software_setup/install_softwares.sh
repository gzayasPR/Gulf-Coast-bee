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
# Install zlib
wget http://www.zlib.net/zlib-1.3.1.tar.gz
tar -xzvf zlib-1.3.1.tar.gz
cd zlib-1.3.1
./configure --prefix=${my_softwares}
make
make install
cd ..

cd ${my_softwares}
# Install bzip2 with -fPIC
wget https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz
tar -xzvf bzip2-1.0.8.tar.gz
cd bzip2-1.0.8
make CFLAGS="-fPIC"
make install PREFIX=${my_softwares}
cd ..

cd ${my_softwares}
# Install xz
wget https://tukaani.org/xz/xz-5.2.5.tar.gz
tar -xzvf xz-5.2.5.tar.gz
cd xz-5.2.5
./configure --prefix=${my_softwares}
make
make install
cd ..

cd ${my_softwares}
# Install OpenSSL (required for curl with TLS support)
wget https://www.openssl.org/source/openssl-1.1.1l.tar.gz
tar -xzvf openssl-1.1.1l.tar.gz
cd openssl-1.1.1l
./config --prefix=${my_softwares} --openssldir=${my_softwares}/ssl
make
make install
cd ..

cd ${my_softwares}
# Install curl with OpenSSL
wget https://curl.se/download/curl-7.79.1.tar.gz
tar -xzvf curl-7.79.1.tar.gz
cd curl-7.79.1
./configure --prefix=${my_softwares} --with-zlib=${my_softwares} --with-ssl=${my_softwares}
make
make install
cd ..

cd ${my_softwares}
git clone https://github.com/samtools/htslib.git
cd htslib
git submodule update --init --recursive
autoreconf -i
./configure --prefix=${my_softwares} --with-zlib=${my_softwares} --with-bzip2=${my_softwares} --with-lzma=${my_softwares} --with-curl=${my_softwares}
make
make install
cd ..

cd ${my_softwares}
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

# Install KING
cd ${my_softwares}
mkdir king
cd king
wget https://www.kingrelatedness.com/KINGcode.tar.gz
tar -xzvf KINGcode.tar.gz
c++ -lm -lz -O2 -fopenmp -o king *.cpp

