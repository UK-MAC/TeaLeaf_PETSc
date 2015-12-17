#define EXTERNAL_FACE -1

MODULE PETScTeaLeaf

  USE definitions_module

  IMPLICIT NONE

#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscdm.h"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"

  INTEGER :: perr
  INTEGER :: mpisize

  ! The maximum number of multigrid levels
  ! It will automatically stop when it coarsens
  ! down to 1 global cell
  PetscInt,parameter :: max_nlevels=10
  !PetscInt,parameter :: max_nlevels=2

  KSP :: kspObj
  PC  :: pcObj
  Vec :: Sol
  Vec :: X
  Vec :: B
  Mat :: A
  Vec :: XLoc
  Vec :: RHSLoc
  
  ! The DMs representing each coarse grid
  DM  :: petscDA(max_nlevels)
  ! The restrictors are stored in this structure
  Mat :: ZT(max_nlevels)
  ! The prolongators are stored in this structure
  Mat :: ZTA(max_nlevels)
  ! We're not currently using this
  ! PETSc is creating the prolongators
  Mat :: Z(max_nlevels)
  ! We're not currently using this
  ! PETSc is creating the coarse space matrices
  Mat :: E(max_nlevels)

CONTAINS

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE setup_petsc(eps,max_iters)

  USE data_module

#include "finclude/petscsys.h"

   ! ~~~~~~~~~~~~~

  INTEGER :: max_iters
  REAL(kind=8) :: eps

  INTEGER :: nlevels, our_level, petsc_level, levels

  PetscInt :: refine_x=2, refine_y=2, refine_z=1
  PetscInt :: actual_refine_x, actual_refine_y, actual_refine_z
  PetscInt :: x_cells, y_cells, z_cells
  PetscInt :: x_cells_old, y_cells_old, z_cells_old
  PetscInt,allocatable,dimension(:) :: lx_coarse,ly_coarse
  PetscViewer :: viewer
  PetscReal :: ztafill=5.0,efill=1.0

  Vec :: rowsum
  PetscInt :: location_min, location_max
  PetscReal :: value_min, value_max

  ! ~~~~~~~~~~~~~

  ! Initialise PETSc
  CALL PetscInitialize(PETSC_NULL_CHARACTER,perr)
  
  ! ~~~~~~~~~~~~~
  ! Create a PETSc DM object to represent the structured grid
  ! ~~~~~~~~~~~~~

  CALL DMDACreate2D(PETSC_COMM_WORLD,                      &
#if PETSC_VER == 342
                    DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, &
#else
                    DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,     &
#endif
                    DMDA_STENCIL_STAR,                     &
                    grid%x_cells,                          &
                    grid%y_cells,                          &
                    px,py,                                 &
                    1,1,                                   &
                    lx(1:px),                              &
                    ly(1:py),                              &
                    petscDA(1),perr)

  ! ~~~~~~~~~~~~~   
  ! Use the DM coarsen to create the coarse grids 
  ! ~~~~~~~~~~~~~
   
  ! ~~~~~~~~~~~~~
  ! Note we are going to introduce two different ways of 
  ! numbering the multigrid levels
  ! our_level and petsc_level
  ! This may seem ridiculous but is done deliberately 
  ! in order to abstract away from the internal petsc level
  ! numbering for when we want to start doing our own nested cycles
  ! that may not fit into standard PETSc cycles
  ! 
  ! PETSc MG levels work in the opposite order to those in our code. If there are N levels:
  !           PETSc   our_level
  ! Fine     N - 1      1
  !          N - 2      2
  !            .        .
  !            1      N - 1
  ! Coarse     0        N
    
  ! Therefore  ---  petsc_level = N - our_level  
  ! ~~~~~~~~~~~~~   

  ! Loop over the levels and coarsen the DM 
  x_cells_old = grid%x_cells
  y_cells_old = grid%y_cells
  !z_cells_old = grid%z_cells

  actual_refine_x=refine_x;
  actual_refine_y=refine_y;
  actual_refine_z=refine_z;

  level_loop_1: DO our_level = 2, max_nlevels
   
    CALL DMDASetRefinementFactor(petscDA(our_level-1),       &
                                 refine_x,                   &
                                 refine_y,                   &
                                 refine_z,                   &
                                 perr)
                                 
    ! need this to force the coarsening factors to be set from the refinement factors                                 
    CALL DMSetFromOptions(petscDA(our_level-1),perr) 

    ! get the refinement factors actually being used, set from the options file or command line
    !CALL DMDAGetRefinementFactor(petscDA(our_level-1),              &
    !                             actual_refine_x,                   &
    !                             actual_refine_y,                   &
    !                             actual_refine_z,                   &
    !                             perr)

    if (x_cells_old > px*actual_refine_x .and. y_cells_old > py*actual_refine_y) then
      CALL DMCoarsen(petscDA(our_level-1),PETSC_COMM_WORLD,petscDA(our_level),perr)
    
      ! Get the number of cells in x, y and z on the coarse level
      call DMDAGetInfo(petscDA(our_level), &
                       PETSC_NULL_INTEGER, &
                       x_cells, &
                       y_cells, &
                       z_cells, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, &
                       perr)

      ! We have to resort to getting the information from successive levels
      ! rather than from the DMDAGetRefinementFactor due to problems with a
      ! missing fortran interface
      actual_refine_x=(x_cells_old+x_cells-1)/x_cells;
      actual_refine_y=(y_cells_old+y_cells-1)/y_cells;
      !actual_refine_z=z_cells_old/z_cells;
      x_cells_old=x_cells
      y_cells_old=y_cells
      !z_cells_old=z_cells
    else
      allocate(lx_coarse(px),ly_coarse(py))
      lx_coarse=1; ly_coarse=1;
      CALL DMDACreate2D(PETSC_COMM_WORLD,                      &
#if PETSC_VER == 342
                        DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, &
#else
                        DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,     &
#endif
                        DMDA_STENCIL_STAR,                     &
                        px,                                    &
                        py,                                    &
                        px,py,                                 &
                        1,1,                                   &
                        lx_coarse(1:px),                       &
                        ly_coarse(1:py),                       &
                        petscDA(our_level),perr)
      deallocate(lx_coarse,ly_coarse)
      nlevels = our_level
      exit level_loop_1
    endif
        
    ! Otherwise use all levels
    if (our_level == max_nlevels) then
      nlevels = max_nlevels        
    end if
       
  ENDDO level_loop_1 
  
  ! ~~~~~~~~~~~~~   
  ! Build the restrictor on each level
  ! ~~~~~~~~~~~~~
  
  ! Loop over the levels
  DO our_level=2, nlevels
    
    ! Create the restrictor going from our_level to our_level-1
    CALL DMCreateAggregates(petscDA(our_level), petscDA(our_level - 1), ZT(our_level - 1), perr)

    CALL DMCreateGlobalVector(petscDA(our_level), rowsum, perr)
    CALL VecSetFromOptions(rowsum, perr)
    CALL MatGetRowSum(ZT(our_level - 1), rowsum, perr)
    CALL VecMin(rowsum, location_min, value_min, perr)
    CALL VecMax(rowsum, location_max, value_max, perr)
    CALL VecDestroy(rowsum, perr)
    IF (parallel%boss) THEN
      WRITE(6,*) "Level ",our_level-1," rowsum ZT min=",location_min,value_min, &
                                      " rowsum ZT max=",location_max,value_max
    ENDIF

    ! Explicltly go and create the prolongator by transposing
    ! For some reason PETSc isn't letting us just call PCMGSetInterpolation and give
    ! it the restrictor, the doc says it is clever enough to work out whether you are 
    ! passing in the restrictor or prolongator
    CALL MatTranspose(ZT(our_level - 1), MAT_INITIAL_MATRIX, Z(our_level - 1), perr)
      
    CALL DMCreateGlobalVector(petscDA(our_level-1), rowsum, perr)
    CALL VecSetFromOptions(rowsum, perr)
    CALL MatGetRowSum(Z(our_level - 1), rowsum, perr)
    CALL VecMin(rowsum, location_min, value_min, perr)
    CALL VecMax(rowsum, location_max, value_max, perr)
    CALL VecDestroy(rowsum, perr)
    IF (parallel%boss) THEN
      WRITE(6,*) "Level ",our_level-1," rowsum Z  min=",location_min,value_min, &
                                      " rowsum Z  max=",location_max,value_max
    ENDIF
  ENDDO
  
  ! ~~~~~~~~~~~~~   
  ! Setup the KSP Solver
  ! ~~~~~~~~~~~~~  

  CALL KSPCreate(MPI_COMM_WORLD,kspObj,perr)
#if PETSC_VER == 342
  CALL KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,max_iters,perr)
#else
  CALL KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,max_iters,perr)
#endif

  CALL KSPSetCheckNormIteration(kspObj,-1,perr)
  CALL KSPSetFromOptions(kspObj,perr)

  IF(parallel%boss) THEN
    WRITE(g_out,*)
    WRITE(g_out,*)'PETSc Setup:'
    WRITE(g_out,*)'Absolute Tolerance set to', eps
    WRITE(g_out,*)'max_iters set to', max_iters
  ENDIF  
  
  ! ~~~~~~~~~~~~~   
  ! Setup the PC for the ksp
  ! ~~~~~~~~~~~~~    
  call PCCreate(MPI_COMM_WORLD, pcObj, perr)
  
  ! Set the PC type to MG
  call PCSetType(pcObj, PCMG, perr)
  
  ! Enable the user to set the number of levels for the MG PC from the command line
  CALL PCSetFromOptions(pcObj, perr)
  CALL PCMGGetLevels(pcObj, levels, perr)
  if (levels > 0) nlevels=min(levels,nlevels)

  ! Set the number of levels for the MG
  ! You must do this before calling any other PETSc 
  ! multigrid routines
  call PCMGSetLevels(pcObj, nlevels, PETSC_NULL_OBJECT, perr)
  
  ! Tell PETSc to use a Galerkin projection to form the 
  ! coarse space matrices
  call PCMGSetGalerkin(pcObj, .TRUE., perr)

  ! Loop over the levels
  ! from the second grid to the coarsest grid 
  do our_level = 1, nlevels-1
    
    ! Get our_level numbering
    petsc_level = nlevels - our_level
            
    ! Set the prolongator
    ! We don't have to bother setting the restrictor as PETSc
    ! will just take the transpose
    call PCMGSetInterpolation(pcObj, petsc_level, Z(our_level), perr)
        
  end do

  ! Enable the user to set options for the PC from the command line
  CALL PCSetFromOptions(pcObj,perr) 

  ! Set the PC to the KSP
  call KSPSetPC(kspObj, pcObj, perr)
  
  ! ~~~~~~~~~~~~~   
  ! Create the sparsity of the A matrix
  ! and preallocate
  ! setupMatA_petsc actually sets the values
  ! ~~~~~~~~~~~~~    

  CALL MPI_Comm_Size(MPI_COMM_WORLD,mpisize,perr)

  IF(mpisize .EQ. 1) THEN
    CALL DMSetMatType(petscDA(1),'seqaij',perr)
  ELSE
    CALL DMSetMatType(petscDA(1),'mpiaij',perr)
  ENDIF

  CALL DMCreateMatrix(petscDA(1),A,perr)
  
  ! ~~~~~~~~~~~~~   
  ! Create the vectors we need (sol and RHS) from 
  ! the DM on the top grid
  ! ~~~~~~~~~~~~~    

  ! Setup the initial solution vector
  CALL DMCreateGlobalVector(petscDA(1),X,perr)

  ! Duplicate Vector to setup B
  CALL DMCreateGlobalVector(petscDA(1),B,perr)

  ! Local Vector for RHS Vector
  CALL DMCreateLocalVector(petscDA(1),RHSLoc,perr)

  ! Local Vector for X Vector
  CALL DMCreateLocalVector(petscDA(1),XLoc,perr)
  
  ! ~~~~~~~~~~~~~  

  total_cg_iter = 0
  total_cheby_iter = 0
  total_petsc_iter = 0

END SUBROUTINE setup_petsc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE cleanup_petsc()

#include "finclude/petscsys.h"

  CALL PetscFinalize(perr)

END SUBROUTINE cleanup_petsc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE setupSol_petsc(c,rx,ry)

    PetscErrorCode :: ierr
    INTEGER       :: c                                ! What chunk are we solving
    INTEGER       :: x_min,x_max,y_min,y_max
    INTEGER       :: g_xmin, g_xmax, g_ymin, g_ymax
    INTEGER       :: left,right,top,bottom
    INTEGER :: i,j,ilen,count
    REAL(KIND=8),pointer,dimension(:) :: rowdata
    INTEGER,pointer,dimension(:) :: rowloc
    REAL(KIND=8) :: val
    REAL(KIND=8),INTENT(IN)   :: rx,ry
    PetscScalar,pointer :: xv(:,:)

    count = 0
    x_min = chunks(c)%field%x_min
    x_max = chunks(c)%field%x_max
    y_min = chunks(c)%field%y_min
    y_max = chunks(c)%field%y_max
    g_xmin = 1
    g_xmax = grid%x_cells
    g_ymin = 1
    g_ymax = grid%y_cells
    left = chunks(c)%field%left
    right = chunks(c)%field%right
    top = chunks(c)%field%top
    bottom = chunks(c)%field%bottom

    ! u0 = chunks(c)%field%u

    CALL VecZeroEntries(X,perr)

    CALL DMDAVecGetArrayF90(petscDA(1),X,xv,perr)

    DO j = bottom, top
      DO i = left, right
          xv(i-1,j-1) = chunks(c)%field%u((i-left)+1,(j-bottom)+1)
      ENDDO
    ENDDO

    CALL DMDAVecRestoreArrayF90(petscDA(1),X,xv,perr)

END SUBROUTINE setupSol_petsc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE setupRHS_petsc(c,rx,ry)

    PetscErrorCode :: ierr
    INTEGER       :: c                                ! What chunk are we solving
    INTEGER       :: left,right,top,bottom
    INTEGER       :: x_min,x_max,y_min,y_max
    INTEGER       :: g_xmin, g_xmax, g_ymin, g_ymax
    INTEGER :: i,j,ilen,count,i_local,j_local
    REAL(KIND=8),pointer,dimension(:) :: rowdata
    INTEGER,pointer,dimension(:) :: rowloc
    REAL(KIND=8) :: val
    REAL(KIND=8),INTENT(IN)   :: rx,ry
    PetscScalar,pointer :: bv(:,:)

    count = 0
    left = chunks(c)%field%left
    right = chunks(c)%field%right
    top = chunks(c)%field%top
    bottom = chunks(c)%field%bottom
    x_min = chunks(c)%field%x_min
    x_max = chunks(c)%field%x_max
    y_min = chunks(c)%field%y_min
    y_max = chunks(c)%field%y_max
    g_xmin = 1
    g_xmax = grid%x_cells
    g_ymin = 1
    g_ymax = grid%y_cells

    CALL VecZeroEntries(B,perr)

    CALL DMDAVecGetArrayF90(petscDA(1),B,bv,perr)

    DO j = bottom, top
      DO i = left, right

        i_local = (i - left)+1
        j_local = (j - bottom)+1

        val = chunks(c)%field%u(i_local,j_local)

        bv(i-1,j-1) = val

      ENDDO
    ENDDO

    CALL DMDAVecRestoreArrayF90(petscDA(1),B,bv,perr)

END SUBROUTINE setupRHS_petsc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE getSolution_petsc(c)

    INTEGER :: i,j,ilen,rowloc
    PetscScalar :: sol_v(1)
    PetscOffset iss
    INTEGER :: count
    INTEGER       :: x_min,x_max,y_min,y_max
    INTEGER       :: c                                ! What chunk are we solving
    PetscScalar,pointer :: xv(:,:)
    INTEGER       :: left,right,top,bottom

    left = chunks(c)%field%left
    right = chunks(c)%field%right
    top = chunks(c)%field%top
    bottom = chunks(c)%field%bottom

    CALL DMDAVecGetArrayF90(petscDA(1),X,xv,perr)

    DO j = bottom, top
      DO i = left, right
        chunks(c)%field%u((i-left)+1,(j-bottom)+1) = xv(i-1,j-1)
      ENDDO
    ENDDO

    CALL DMDAVecRestoreArrayF90(petscDA(1),X,xv,perr)

END SUBROUTINE getSolution_petsc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE setupMatA_petsc(c,rx,ry)

#include "finclude/petscsys.h"

  INTEGER       :: c                                ! What chunk are we solving
  INTEGER       :: left,right,top,bottom
  INTEGER       :: x_min,x_max,y_min,y_max
  INTEGER       :: g_xmin, g_xmax, g_ymin, g_ymax
  REAL(KIND=8),INTENT(IN)   :: rx,ry
  INTEGER       :: i,j,i_g,j_g
  INTEGER       :: count
  REAL(KIND=8)  :: c1,c2,c3,c4,c5

  MatStencil  :: row(4,1)       ! 4 x Number of stencils in entry (adding 1 at a time for now, so 1)
  MatStencil  :: column(4,5)    ! 4 x Stencil Size: (i,j,k,m) * stencil Size (5 point stencil)
  PetscScalar :: stencil(5)     ! 4 x Stencil Size: (i,j,k,m) * stencil Size (5 point stencil)
  PetscReal   :: ztafill=5.0,efill=1.0

  PetscViewer :: viewer

  left = chunks(c)%field%left
  right = chunks(c)%field%right
  top = chunks(c)%field%top
  bottom = chunks(c)%field%bottom

  x_min = chunks(c)%field%x_min
  x_max = chunks(c)%field%x_max
  y_min = chunks(c)%field%y_min
  y_max = chunks(c)%field%y_max

  g_xmin = 1
  g_xmax = grid%x_cells
  g_ymin = 1
  g_ymax = grid%y_cells

  ! Zero Matrix
  CALL MatZeroEntries(A,perr)

  DO j_g = bottom,top
    DO i_g = left,right

        ! -1 is applied to all indexes to shift them to zero-index, presuming the fortran code maintains
        ! an index from 1 to max

        i = (i_g - left) + 1
        j = (j_g - bottom) + 1

        count = 1
        row(MatStencil_i,1) = i_g - 1
        row(MatStencil_j,1) = j_g - 1

        c2 = (-1*rx) * chunks(c)%field%Vector_Kx(i,j)
        c3 = (-1*rx) * chunks(c)%field%Vector_Kx(i+1,j)
        c4 = (-1*ry) * chunks(c)%field%Vector_Ky(i,j)
        c5 = (-1*ry) * chunks(c)%field%Vector_Ky(i,j+1)

        c1 = (1.0_8-c2-c3-c4-c5)

        IF(.NOT.(j_g == 1)) THEN     ! Not Bottom External Boundary
          stencil(count) = c4
          column(MatStencil_i,count) = i_g-1
          column(MatStencil_j,count) = (j_g-1) - 1
          count = count + 1
        ENDIF

        IF(.NOT.(i_g == 1)) THEN     ! Not Left External Boundary
          stencil(count) = c2
          column(MatStencil_i,count) = (i_g-1) - 1
          column(MatStencil_j,count) = j_g - 1
          count = count + 1
        ENDIF

        ! Cell Centered Value
        stencil(count) = c1
        column(MatStencil_i,count) = i_g - 1
        column(MatStencil_j,count) = j_g - 1
        count = count + 1

        IF(.NOT.(i_g == grid%x_cells)) THEN     ! Not Right External Boundary
          stencil(count) = c3
          column(MatStencil_i,count) = (i_g+1) - 1
          column(MatStencil_j,count) = j_g - 1
          count = count + 1
        ENDIF

        IF(.NOT.(j_g == grid%y_cells)) THEN     ! Not Top External Boundary
          stencil(count) = c5
          column(MatStencil_i,count) = i_g - 1
          column(MatStencil_j,count) = (j_g+1) - 1
          count = count + 1
        ENDIF

        CALL MatSetValuesStencil(A,1,row,count-1,column,stencil,INSERT_VALUES,perr)

    ENDDO
  ENDDO

  CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,perr)
  CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,perr)

END SUBROUTINE setupMatA_petsc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE solve_petsc(numit,error)

    !USE MPI
    USE data_module

    include "mpif.h"

    INTEGER,INTENT(INOUT) :: numit
    REAL(KIND=8),INTENT(INOUT) :: error
    INTEGER :: errcode, mpierr
    KSPConvergedReason :: reason
    PetscReal :: residual  ! Final residual

    CALL KSPSetOperators(kspObj,A,A,perr)
    CALL KSPSolve(kspObj,B,X,perr)
    CALL KSPGetIterationNumber(kspObj, numit, perr)

    CALL KSPGetConvergedReason(kspObj,reason,perr)
    CALL KSPGetResidualNorm(kspObj,residual,perr)
    error=residual

    IF(reason < 0) THEN
      WRITE(g_out,*) ' Error: Did not Converge. Calling MPI_Abort'
      WRITE(g_out,*) ' Divergence Reason:'

      IF(reason == -2) THEN
        WRITE(g_out,*) ' Diverged Null'
      ELSE IF(reason == -3) THEN
        WRITE(g_out,*) ' Diverged ITS, took ', numit, ' iterations'
      ELSE IF(reason == -4) THEN
        WRITE(g_out,*) ' Diverged DTOL'
      ELSE IF(reason == -5) THEN
        WRITE(g_out,*) ' Diverged Breakdown'
      ELSE IF(reason == -6) THEN
        WRITE(g_out,*) ' Diverged Breakdown BICG'
      ELSE IF(reason == -7) THEN
        WRITE(g_out,*) ' Diverged NonSymmetric'
      ELSE IF(reason == -8) THEN
        WRITE(g_out,*) ' Diverged Indefinite PC'
      ELSE IF(reason == -9) THEN
        WRITE(g_out,*) ' Diverged nanorinf'
      ELSE IF(reason == -10) THEN
        WRITE(g_out,*) ' Diverged indefinite mat'
      ENDIF

      CALL MPI_Abort(MPI_COMM_WORLD,errcode,mpierr)
    ENDIF 

    total_petsc_iter = total_petsc_iter + numit

END SUBROUTINE solve_petsc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Apply Paul Garrett Approach
! (1) Use CG to Execute 10 iteration, retrieve eigenvalues
! (2) Set Solver to Chebyshev type, set eigenvalues from CG execution
! (3) Set Residual for Chebyshev to reduce iteration count
! (4) Execute Chebyshev Solve
SUBROUTINE solve_petsc_pgcg(eps,max_iters,numit_cg,numit_cheby,error)

!    USE MPI
    USE data_module

#include "finclude/petscsys.h"

    INTEGER,INTENT(INOUT) :: numit_cg, numit_cheby
    INTEGER :: errcode, mpierr
    KSPConvergedReason :: reason
    PC :: tPC
    REAL(KIND=8) :: eps,error
    INTEGER :: max_iters
    PetscReal :: r(0:pgcg_cg_iter-1)    ! Real Component of EigenValue Array
    PetscReal :: c(0:pgcg_cg_iter-1)    ! Complex Component of EigenValue Array
    PetscInt  :: neig      ! Number of Eigenvalues computed
    PetscReal :: emax      ! Max EigenValue
    PetscReal :: emin      ! Min EigenValue
    PetscReal :: residual  ! Final residual

    IF(parallel%boss) WRITE(g_out,*) 'pgcg_cg_iter set to ', pgcg_cg_iter

    CALL KSPSetComputeEigenValues(kspObj,PETSC_TRUE,perr) ! Must be called before KSPSetup Routines to enable eigenvalue calculation

    CALL KSPSetType(kspObj,KSPCG,perr)
    CALL KSPGetPC(kspObj,tPC,perr)
    CALL PCSetType(tPC,PCBJACOBI,perr)

    CALL KSPSetOperators(kspObj,A,A,perr)
#if PETSC_VER == 342
    CALL KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,pgcg_cg_iter,perr)
#else
    CALL KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,pgcg_cg_iter,perr)
#endif
    CALL KSPSetInitialGuessNonzero(kspObj,PETSC_FALSE,perr)

    CALL KSPSolve(kspObj,B,X,perr)

    ! Don't check convergence reasons or iteration count here - unlikely to have converged yet

    ! Calculate the EigenValues from CG Solve after 10 iterations. Array should be returned sorted.
    CALL KSPComputeEigenValues(kspObj,pgcg_cg_iter,r,c,neig,perr)
    emax = r(neig-1)
    emin = r(0)

    IF(parallel%boss) WRITE(g_out,*) ' Emax:', emax, ', Emin:', emin

    CALL KSPGetIterationNumber(kspObj, numit_cg, perr)
    total_cg_iter = total_cg_iter + numit_cg

    CALL KSPGetConvergedReason(kspObj,reason,perr)
    CALL KSPGetResidualNorm(kspObj,residual,perr)
    error=residual

    IF(reason < 0) THEN   ! Did not converge in CG, Progress onto Cheby

      ! Main Solve is in Next section within Chebyshev

      CALL KSPSetComputeEigenValues(kspObj,PETSC_FALSE,perr) ! Disable EigenValue Calculation

      CALL KSPSetType(kspObj,KSPCHEBYSHEV,perr)
      CALL KSPGetPC(kspObj,tPC,perr)
      CALL PCSetType(tPC,PCBJACOBI,perr)

      CALL KSPSetInitialGuessNonzero(kspObj,PETSC_TRUE,perr)    ! Disable zeroing of results vector (reuse for residual)
      CALL KSPSetOperators(kspObj,A,A,perr)
#if PETSC_VER == 342
      CALL KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,max_iters,perr)
#else
      CALL KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,max_iters,perr)
#endif
      CALL KSPChebyshevSetEigenValues(kspObj,emax,emin,perr)

      CALL KSPSolve(kspObj,B,X,perr)

      CALL KSPGetIterationNumber(kspObj, numit_cheby, perr)

      CALL KSPGetConvergedReason(kspObj,reason,perr)
      CALL KSPGetResidualNorm(kspObj,residual,perr)
      error=residual
      IF(reason < 0) THEN
        WRITE(g_out,*) ' Error: Did not Converge. Calling MPI_Abort'
        WRITE(g_out,*) ' Divergence Reason:'

        IF(reason == -2) THEN
          WRITE(g_out,*) ' Diverged Null'
        ELSE IF(reason == -3) THEN
          WRITE(g_out,*) ' Cheby - Diverged ITS, took ', numit_cheby, ' iterations'
        ELSE IF(reason == -4) THEN
          WRITE(g_out,*) ' Diverged DTOL'
        ELSE IF(reason == -5) THEN
          WRITE(g_out,*) ' Diverged Breakdown'
        ELSE IF(reason == -6) THEN
          WRITE(g_out,*) ' Diverged Breakdown BICG'
        ELSE IF(reason == -7) THEN
          WRITE(g_out,*) ' Diverged NonSymmetric'
        ELSE IF(reason == -8) THEN
          WRITE(g_out,*) ' Diverged Indefinite PC'
        ELSE IF(reason == -9) THEN
          WRITE(g_out,*) ' Diverged nanorinf'
        ELSE IF(reason == -10) THEN
          WRITE(g_out,*) ' Diverged indefinite mat'
        ENDIF

        CALL MPI_Abort(MPI_COMM_WORLD,errcode,mpierr)
      ENDIF 

      total_cheby_iter = total_cheby_iter + numit_cheby

    ENDIF

END SUBROUTINE solve_petsc_pgcg

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE printXVec(fileName)

    !USE MPI

    IMPLICIT NONE

    include "mpif.h"

#   include "finclude/petscviewer.h"

    CHARACTER(LEN=*) :: fileName
    CHARACTER(LEN=8) :: id
    PetscViewer :: viewer

    ! Write Global Solution Vector
    CALL PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    CALL PetscViewerSetType(viewer,"ascii",perr)
    CALL PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX,perr)
    CALL PetscViewerFileSetName(viewer,fileName,perr)
    CALL VecView(X,viewer,perr)
    CALL PetscViewerDestroy(viewer,perr)

END SUBROUTINE printXVec

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE printBVec(fileName)

    !USE MPI

    IMPLICIT NONE

    include "mpif.h"

#   include "finclude/petscviewer.h"

    CHARACTER(LEN=*) :: fileName
    CHARACTER(LEN=8) :: id
    PetscViewer :: viewer

    CALL PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    CALL PetscViewerSetType(viewer,"ascii",perr)
    CALL PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX,perr)
    CALL PetscViewerFileSetName(viewer,fileName,perr)
    CALL VecView(B,viewer,perr)
    CALL PetscViewerDestroy(viewer,perr)

END SUBROUTINE printBVec

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE printMatA(fileName)

    !USE MPI

    IMPLICIT NONE

    include "mpif.h"

#   include "finclude/petscviewer.h"

    CHARACTER(LEN=*) :: fileName
    CHARACTER(LEN=8) :: id
    PetscViewer :: viewer

    CALL PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    CALL PetscViewerSetType(viewer,"ascii",perr)
    CALL PetscViewerSetFormat(viewer,PETSC_VIEWER_DEFAULT,perr)
    CALL PetscViewerFileSetName(viewer,fileName,perr)
    CALL MatView(A,viewer,perr)
    CALL PetscViewerDestroy(viewer,perr)

END SUBROUTINE printMatA

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE PETScTeaLeaf
