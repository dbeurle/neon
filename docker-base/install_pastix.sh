#!/bin/bash

# Download and extract the stable PaStiX version
wget https://gforge.inria.fr/frs/download.php/latestfile/218/pastix_5.2.3.tar.bz2
tar -xvf pastix_5.2.3.tar.bz2
ls
cd pastix_5.2.3

# Build (only can do this once due to some makefile bug)
patch -p1 < ../cmake_patch.diff
mkdir build; cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPASTIX_WITH_MPI=OFF -DPASTIX_INT64=OFF -DBUILD_SHARED_LIBS=ON ../src
make all -j4 && make install

ln -s /usr/local/lib/libpastix.so /usr/lib64/libpastix.so
