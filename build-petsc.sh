export LANG=C
export LC_ALL=C

# How you select a petsc compiler
export I_MPI_CC="icc -lmpi"
export I_MPI_CXX="icpc -lmpi"
export I_MPI_F77="ifort -lmpi"
export I_MPI_F90="ifort -lmpi"

BLAS_FLAGS="-mkl=sequential"
MPI_DIR=/sw/sdev/mpt-x86_64/2.11-p11155
PREFIX=/store/jsouthern/packages/petsc/dev

# Get PETSc from Bitbucket
cd ${HOME}
git clone https://bitbucket.org/petsc/petsc petsc-dev

cd petsc-dev

./configure --prefix=${PREFIX} --with-c++-support=1 --with-c-support=1 \
  --with-cc=${I_MPI_CC} --with-cxx=${I_MPI_CXX} --with-fc=${I_MPI_F77} \
  --with-fortran=1 --with-x11=no --with-mpi=1 --with-hypre=0 --with=spooles=0 \
  --with-ml=0 --with-mpiexec=mpirun --with-blas-lapack-lib=${BLAS_FLAGS} \
  --with-debugging=no PETSC_ARCH=linux-static --with-shared-libraries=0 \
  --with-x=0 --with-threadcomm=1 --with-openmp=1 --with-pthreadclasses=1

export PETSC_ARCH=linux-static
export PETSC_DIR=${PWD}

make all
make install
