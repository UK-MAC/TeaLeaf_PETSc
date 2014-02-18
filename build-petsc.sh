export LANG=C
export LC_ALL=C

BLAS_DIR=
MPI_DIR=
PREFIX=.

./configure --with-prefix=${PREFIX} --with-c++-support=1 --with-c-support=1 \
  --with-fortran=1 --with-x11=no --with-mpi=1 --with-hypre=0 --with=spooles=0 \
  --with-ml=0 --with-mpi-dir=${MPI_DIR} --with-blas-lapack-dir=${BLAS_DIR}/lib \
  --with-debugging=no PETSC_ARCH=linux-static

export PETSC_ARCH=linux-static
export PETSC_DIR=${PWD}

make all
make install
