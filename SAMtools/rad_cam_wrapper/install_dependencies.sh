#!/bin/bash

# Install zlib
tar -xvzf zlib-1.2.11.tar.gz
mkdir zlib-1.2.11-install
cd zlib-1.2.11-install
ZDIR=`pwd`
cd ../zlib-1.2.11
./configure
make test
make install prefix=${ZDIR}
cd ..
mv zlib-1.2.11-install zlib-1.2.11
cd zlib-1.2.11
ZDIR=`pwd`
cd ..

# Install HDF5
tar -xvzf hdf5-1.10.2.tar.gz
mkdir hdf5-1.10.2-install
cd hdf5-1.10.2-install
H5DIR=`pwd`
cd ../hdf5-1.10.2
./configure --with-zlib=${ZDIR} --prefix=${H5DIR}
make check
make install
cd ..
mv hdf5-1.10.2-install hdf5-1.10.2
cd hdf5-1.10.2
H5DIR=`pwd`
cd ..

# Install NetCDF C libraries
tar -xvzf netcdf-c-4.6.1.tar.gz
mkdir netcdf-c-4.6.1-install
cd netcdf-c-4.6.1-install
NCDIR=`pwd`
cd ../netcdf-c-4.6.1
CPPFLAGS=-I${H5DIR}/include LDFLAGS=-L${H5DIR}/lib ./configure --prefix=${NCDIR}
make check
make install
cd ..
mv netcdf-c-4.6.1-install netcdf-c-4.6.1
cd netcdf-c-4.6.1
NCDIR=`pwd`
cd ..

# Install NetCDF Fortran libraries
tar -xvzf netcdf-fortran-4.4.4.tar.gz
mkdir netcdf-fortran-4.4.4-install
cd netcdf-fortran-4.4.4-install
NFDIR=`pwd`
CPPFLAGS=-I${NCDIR}/include -LDFLAGS=-L${NCDIR}/lib ./configure --prefix=${NFDIR}
make check
make install
cd ..
mv netcdf-fortran-4.4.4-install netcdf-fortran-4.4.4
cd netcdf-fortran-4.4.4
NFDIR=`pwd`
cd ..
