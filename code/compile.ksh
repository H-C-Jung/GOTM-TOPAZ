#!/bin/ksh
# Change FC, CC, FLAGS : ./compilers/compiler.${FORTRAN_COMPILER}

order=${1}

export GOTMDIR=/home/hcjung/wdata/Model/GOTM-TOPAZ_v1.0/code

export NETCDFHOME=/usr/local/netcdf/363_intel14
export NETCDFINC=${NETCDFHOME}/include
export NETCDFLIBDIR=${NETCDFHOME}/lib

export FORTRAN_COMPILER=IFORT

#order :: exec, distclean 
make ${order} -f Makefile
