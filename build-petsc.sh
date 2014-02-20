export LANG=C
export LC_ALL=C

# How you select a petsc compiler
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

BLAS_DIR=location_of_BLAS_on_your_system
MPI_DIR=location_of_mpi_on_your_system
PREFIX=where_you_want_petsc_installed

BLAS_DIR=/opt/lapack/3.4.2/intel-13.1.1.163/
MPI_DIR=/opt/intel/impi/4.1.0.024/intel64/
PREFIX=/home/wpg/TeaLeaf/tealeaf_petsc

./configure --prefix=${PREFIX} --with-c++-support=1 --with-c-support=1 \
  --with-fortran=1 --with-x11=no --with-mpi=1 --with-hypre=0 --with=spooles=0 \
  --with-ml=0 --with-mpi-dir=${MPI_DIR} --with-blas-lapack-dir=${BLAS_DIR}/lib \
  --with-debugging=no PETSC_ARCH=linux-static --with-shared-libraries=0

export PETSC_ARCH=linux-static
export PETSC_DIR=${PWD}

make all
make install
