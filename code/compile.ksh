#!/bin/ksh
# Change FC, CC, FLAGS : ./compilers/compiler.${FORTRAN_COMPILER}

order=${1}

export SCENARIO=eastSea

export GOTMDIR=/home/hcjung/wdata/Model/GOTM-TOPAZ_v1.0/code

export FORTRAN_COMPILER=IFORT

export NETCDFHOME=/usr/local/netcdf/363_intel14
export NETCDFINC=${NETCDFHOME}/include
export NETCDFLIBDIR=${NETCDFHOME}/lib

export MPILIB=/usr/local/mpi/intel/mvapich2-1.4_14

#order :: exec, distclean 
make ${order} -f Makefile

ln -sf ${GOTMDIR}/src/gotm_prod_IFORT ../test_case/${SCENARIO}/gotm
