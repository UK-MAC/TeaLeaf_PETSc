# Example of how do you download, install and compile PETSc 3.5.2, assuming installing in /home/usid, with MPICH
# using the GNU compiler. It could work with others but might need maths libs in lib paths.
# The variables OTHER_LIBS and HYPRE_LIB can be used to set these up

cd /home/sdargavi/projects/dependencies
git clone -b maint https://bitbucket.org/petsc/petsc petsc-3.5.2
cd petsc-3.5.2
./configure --download-fblaslapack --download-hypre --download-mpich --with-debugging=0 --with-c2html=0
make PETSC_DIR=/home/sdargavi/projects/dependencies/petsc-3.5.2 PETSC_ARCH=arch-linux2-c-opt all
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/usid/petsc-3.5.2/arch-linux2-c-opt/lib/
cd /home/sdargavi/projects/TeaLeaf_PETSc
COM_PATH_P=/home/sdargavi/projects/dependencies/petsc-3.5.2 make COMPILER=GNU

