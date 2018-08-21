# #!/bin/bash
# 
# # Install zlib
# tar -xvzf zlib-1.2.11.tar.gz
# mkdir zlib
# cd zlib; ZDIR=`pwd`; cd ..
# cd zlib-1.2.11
# ./configure
# make test
# make install prefix=${ZDIR}
# cd ..
# rm -r zlib-1.2.11
# 
# # Install HDF5
# tar -xvzf hdf5-1.10.2.tar.gz
# mkdir hdf5
# cd hdf5; H5DIR=`pwd`; cd ..
# cd hdf5-1.10.2
# ./configure --with-zlib=${ZDIR} --prefix=${H5DIR}
# make check
# make install
# cd ..
# rm -r hdf5-1.10.2
# 
# # Install NetCDF C libraries
# tar -xvzf netcdf-c-4.6.1.tar.gz
# mkdir netcdf-c
# cd netcdf-c; NCDIR=`pwd`; cd ..
# cd netcdf-c-4.6.1
# CPPFLAGS=-I${H5DIR}/include LDFLAGS=-L${H5DIR}/lib ./configure --prefix=${NCDIR}
# make check
# make install
# cd ..
# rm -r netcdf-c-4.6.1

cd netcdf-c
NCDIR=`pwd`
cd ..

# Install NetCDF Fortran libraries
tar -xvzf netcdf-fortran-4.4.4.tar.gz
mkdir netcdf-fortran
cd netcdf-fortran; NFDIR=`pwd`; cd ..
cd netcdf-fortran-4.4.4
CPPFLAGS=-I${NCDIR}/include LDFLAGS=-L${NCDIR}/lib ./configure --prefix=${NFDIR}
make check
make install
cd ..
rm -r netcdf-fortran-4.4.4
