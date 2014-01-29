## PETSc

- Download a copy of PETSc from mcs.anl.gov/petsc/, we use version 3.4.3
- Extract the downloaded copy of PETSc, and copy in the `build-petsc.sh` found in the TeaLeaf repository.
- Edit the `build-petsc.sh` script and change the BLAS_DIR, MPI_DIR and INSTALL_DIR environment variables.
  - If you are only building PETSc for use with TeaLeaf, the consider installing to `libs/petsc` in the TeaLeaf directory.
- Run the script to build and install PETSc.
