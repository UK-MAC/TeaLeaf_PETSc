# Example of how do you download, install and compile PETSc 3.14.3, assuming installing in /home/usid, with MPICH
# using the GNU compiler. It could work with others but might need maths libs in lib paths.
# The variables OTHER_LIBS and HYPRE_LIB can be used to set these up

cd /home/usid
git clone -b maint https://bitbucket.org/petsc/petsc petsc-3.14.3
cd petsc-3.14.3
./configure --download-fblaslapack --download-hypre --download-mpich --with-debugging=0 --with-c2html=0
make PETSC_DIR=/home/usid/petsc-3.14.3 PETSC_ARCH=arch-linux2-c-opt all
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/usid/petsc-3.14.3/arch-linux2-c-opt/lib/
cd /home/usid/TeaLeaf_PETSc
COM_PATH_P=home/usid/petsc-3.14.3 make COMPILER=GNU

