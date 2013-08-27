#define EXTERNAL_FACE -1

module PETScTeaLeaf

  use definitions_module

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

  KSP :: kspObj
  Vec :: Sol
  Vec :: X
  Vec :: B
  Mat :: A
  Vec :: XLoc
  Vec :: RHSLoc
  DM  :: petscDA

contains

subroutine setup_petsc(eps,max_iters)

#include "finclude/petscsys.h"

  INTEGER :: c,cx,cy
  INTEGER :: lx(1:px)
  INTEGER :: ly(1:py)
  REAL(kind=8) :: eps
  INTEGER :: max_iters

  call PetscInitialize(PETSC_NULL_CHARACTER,perr)

  ! px, py set in clover_decompose to be chunks_x and chunks_y
  ! clover_decompose MUST be called first

  DO cx=1,px
    c = cx
    lx(cx) = (chunks(c)%field%right - chunks(c)%field%left)+1
    !write(6,*) 'lx(',cx,')=',lx(cx)
  ENDDO

  DO cy=1,py
    c = ((cy-1)*px)+1
    ly(cy) = (chunks(c)%field%top - chunks(c)%field%bottom)+1
    !write(6,*) 'ly(',cy,')=',ly(cy)
  ENDDO

 !write(6,*) 'Px=',px
  !write(6,*) 'Py=',py
  !write(6,*) 'grid%x_cells=',grid%x_cells
  !write(6,*) 'grid%y_cells=',grid%y_cells


  call DMDACreate2D(PETSC_COMM_WORLD, &
                    DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE, &
                    DMDA_STENCIL_STAR, &
				            grid%x_cells,grid%y_cells, &
                    px,py, &
				            1,1, &
                    lx(1:px),ly(1:py), &
                    petscDA,perr)

  ! Setup the KSP Solver
  call KSPCreate(MPI_COMM_WORLD,kspObj,perr)
  call KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,max_iters,perr)
  call KSPSetFromOptions(kspObj,perr)

  write(6,*)'PETSc Setup:'
  write(6,*)'Absolute Tolerance set to', eps
  write(6,*)'max_iters set to', max_iters

  call MPI_Comm_Size(MPI_COMM_WORLD,mpisize,perr)

  if(mpisize .eq. 1) then
    call DMCreateMatrix(petscDA,'seqaij',A,perr)
  else
    call DMCreateMatrix(petscDA,'mpiaij',A,perr)
  endif

  ! Setup the initial solution vector
  call DMCreateGlobalVector(petscDA,X,perr)

  ! Duplicate Vector to setup B
  call DMCreateGlobalVector(petscDA,B,perr)

  ! Local Vector for RHS Vector
  call DMCreateLocalVector(petscDA,RHSLoc,perr)

  ! Local Vector for X Vector
  call DMCreateLocalVector(petscDA,XLoc,perr)

end subroutine setup_petsc

subroutine cleanup_petsc()

#include "finclude/petscsys.h"

  call PetscFinalize(perr)

end subroutine cleanup_petsc


subroutine setupSol_petsc(c,rx,ry)

    PetscErrorCode :: ierr
    INTEGER       :: c                                ! What chunk are we solving
    INTEGER       :: x_min,x_max,y_min,y_max
    INTEGER       :: g_xmin, g_xmax, g_ymin, g_ymax
    INTEGER       :: left,right,top,bottom
    INTEGER :: i,j,ilen,count
    REAL(kind=8),pointer,dimension(:) :: rowdata
    INTEGER,pointer,dimension(:) :: rowloc
    REAL :: val
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

    call VecZeroEntries(X,perr)

!    allocate(rowdata(0:((x_max-x_min)+1) * ((y_max-y_min)+1)))
!    allocate(rowloc(0:((x_max-x_min)+1) * ((y_max-y_min)+1)))

!    ilen = ((x_max+1)-(x_min-1))+1

!    if(parallel%task .eq. 1) write(6,*) 'x_min:', x_min
!    if(parallel%task .eq. 1) write(6,*) 'x_max:', x_max

!    if(parallel%task .eq. 1) write(6,*) 'y_min:', y_min
!    if(parallel%task .eq. 1) write(6,*) 'y_max:', y_max

!    do j = y_min, y_max
!      do i = x_min, x_max
!        rowdata(count) = chunks(c)%field%u(i,j)
!        ! -1 from indexes for x_min and y_min for ghosts - local vector has ghost rows as well so
!        ! row position has to have an extra offset from 0
!        rowloc(count) = ((i-1)-(x_min-1)) +(((j-1)-(y_min-1))  * ilen)
!        if(parallel%task .eq. 1) write(6,*)'rowloc:',rowloc(count)
!        count = count + 1
!      enddo
!    enddo

    !call VecSetValues(XLoc,count,rowloc,rowdata,INSERT_VALUES,perr);
    !call DMLocalToGlobalBegin(petscDA,XLoc,INSERT_VALUES,X,perr)
    !call DMLocalToGlobalEnd(petscDA,XLoc,INSERT_VALUES,X,perr)

!    deallocate(rowloc)
!    deallocate(rowdata)

    call DMDAVecGetArrayF90(petscDA,X,xv,perr)

    do j = bottom, top
      do i = left, right
          xv(i-1,j-1) = chunks(c)%field%u((i-left)+1,(j-bottom)+1)
      enddo
    enddo

    call DMDAVecRestoreArrayF90(petscDA,X,xv,perr)



end subroutine setupSol_petsc



subroutine setupRHS_petsc(c,rx,ry)

    PetscErrorCode :: ierr
    INTEGER       :: c                                ! What chunk are we solving
    INTEGER       :: left,right,top,bottom
    INTEGER       :: x_min,x_max,y_min,y_max
    INTEGER       :: g_xmin, g_xmax, g_ymin, g_ymax
    INTEGER :: i,j,ilen,count,i_local,j_local
    REAL(kind=8),pointer,dimension(:) :: rowdata
    INTEGER,pointer,dimension(:) :: rowloc
    REAL :: val
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

    call VecZeroEntries(B,perr)

!    allocate(rowdata(0:((x_max-x_min)+1) * ((y_max-y_min)+1)))
!    allocate(rowloc(0:((x_max-x_min)+1) * ((y_max-y_min)+1)))

!    ilen = ((x_max+1)-(x_min-1))+1

!    do j = bottom, top
!      do i = left, right

!        val = chunks(c)%field%u(i,j)

!        i_local = i - left
!        j_local = j - bottom

!        IF(i == g_xmin) THEN
!          val = val + rx * chunks(c)%field%u(i_local-1,j_local)
!        ELSE IF(i == g_xmax) THEN
!          val = val + rx * chunks(c)%field%u(i_local+1,j_local)
!        ENDIF

!        IF(j == g_ymin) THEN
!          val = val + ry * chunks(c)%field%u(i_local,j_local-1)
!        ELSE IF(j == g_ymax) THEN
!          val = val + ry * chunks(c)%field%u(i_local,j_local+1)
!        ENDIF

!        rowdata(count) = val
!        ! -1 from indexes for x_min and y_min for ghosts - local vector has ghost rows as well so
!        ! row position has to have an extra offset from 0
!        ! -1 from i and j to zero index them for the PETSc indexing system, which is zero indexed
!        rowloc(count) = ((i-1)-(x_min-1)) +(((j-1)-(y_min-1))  * ilen   )     
!        count = count + 1
!      enddo
!    enddo

    !call VecSetValues(RHSLoc,count,rowloc,rowdata,INSERT_VALUES,perr);

    !call DMLocalToGlobalBegin(petscDA,RHSLoc,INSERT_VALUES,B,perr)
    !call DMLocalToGlobalEnd(petscDA,RHSLoc,INSERT_VALUES,B,perr)

!    deallocate(rowloc)
!    deallocate(rowdata)


    call DMDAVecGetArrayF90(petscDA,B,bv,perr)

    do j = bottom, top
      do i = left, right

        i_local = (i - left)+1
        j_local = (j - bottom)+1

        val = chunks(c)%field%u(i_local,j_local)

        ! This bit now commented out in hypre version - boundary stuff moved to matrix

!        IF(i == g_xmin) THEN
!          val = val + rx * chunks(c)%field%work_array6(i_local,j_local) * chunks(c)%field%u(i_local-1,j_local)
!        ELSE IF(i == g_xmax) THEN
!          val = val + rx * chunks(c)%field%work_array6(i_local+1,j_local) * chunks(c)%field%u(i_local+1,j_local)
!        ENDIF

!        IF(j == g_ymin) THEN
!          val = val + ry *  chunks(c)%field%work_array7(i_local,j_local) * chunks(c)%field%u(i_local,j_local-1)
!        ELSE IF(j == g_ymax) THEN
!          val = val + ry * chunks(c)%field%work_array7(i_local,j_local+1) * chunks(c)%field%u(i_local,j_local+1)
!        ENDIF

        bv(i-1,j-1) = val

      enddo
    enddo

    call DMDAVecRestoreArrayF90(petscDA,B,bv,perr)


end subroutine setupRHS_petsc


subroutine getSolution_petsc(c)

    integer :: i,j,ilen,rowloc
    PetscScalar :: sol_v(1)
    PetscOffset iss
    integer :: count
    INTEGER       :: x_min,x_max,y_min,y_max
    INTEGER       :: c                                ! What chunk are we solving
    PetscScalar,pointer :: xv(:,:)
    INTEGER       :: left,right,top,bottom

!#   define sol_a(ib) sol_v(iss+(ib))

!!    call DMGlobalToLocalBegin(petscDA,X,INSERT_VALUES,XLoc,perr)
!!    call DMGlobalToLocalEnd(petscDA,X,INSERT_VALUES,XLoc,perr)

!!    call VecGetArray(XLoc,sol_v,iss,perr)

!    x_min = chunks(c)%field%x_min
!    x_max = chunks(c)%field%x_max
!    y_min = chunks(c)%field%y_min
!    y_max = chunks(c)%field%y_max

!!    ilen = ((x_max+1)-(x_min-1))+1

!!    do j = y_min, y_max
!!      do i = x_min, x_max

!!        ! Local Vector still has ghosts, calculate i,j index
!!        rowloc = ((i-1)-(x_min-1)) +(((j-1)-(y_min-1))  * ilen)
!!        chunks(c)%field%u(i,j) = sol_a(rowloc+1)    ! +1 since indexed from 1

!!      enddo
!!    enddo

!!    call VecRestoreArray(XLoc,sol_v,iss,perr)

!!    call DMLocalToGlobalBegin(petscDA,XLoc,INSERT_VALUES,X,perr)
!!    call DMLocalToGlobalEnd(petscDA,XLoc,INSERT_VALUES,X,perr)

    left = chunks(c)%field%left
    right = chunks(c)%field%right
    top = chunks(c)%field%top
    bottom = chunks(c)%field%bottom

    call DMDAVecGetArrayF90(petscDA,X,xv,perr)


    do j = bottom, top
      do i = left, right
        chunks(c)%field%u((i-left)+1,(j-bottom)+1) = xv(i-1,j-1)
      enddo
    enddo

    call DMDAVecRestoreArrayF90(petscDA,X,xv,perr)


end subroutine getSolution_petsc


subroutine setupMatA_petsc(c,rx,ry)

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
  call MatZeroEntries(A,perr)

  ! Kx = work_array6
  ! Ky = work_array7

  DO j_g = bottom,top
    DO i_g = left,right

        ! -1 is applied to all indexes to shift them to zero-index, presuming the fortran code maintains
        ! an index from 1 to max

        i = (i_g - left) + 1
        j = (j_g - bottom) + 1

        count = 1
        row(MatStencil_i,1) = i_g - 1
        row(MatStencil_j,1) = j_g - 1

        c1 = (1.0+(2.0*rx)+(2.0*ry))
        c2 = (-1*rx) * chunks(c)%field%work_array6(i,j)
        c3 = (-1*rx) * chunks(c)%field%work_array6(i+1,j)
        c4 = (-1*ry) * chunks(c)%field%work_array7(i,j)
        c5 = (-1*ry) * chunks(c)%field%work_array7(i,j+1)

        IF(i_g == 1) THEN ! Global X Low Boundary
          c3 = (-2*rx) * chunks(c)%field%work_array6(i+1,j)
        ENDIF

        IF(i_g == grid%x_cells) THEN
          c2 = (-2*rx) * chunks(c)%field%work_array6(i,j)
        ENDIF
        
        IF(j_g == 1) THEN
          c5 = (-2*ry) * chunks(c)%field%work_array7(i,j+1)
        ENDIF

        IF(j_g == grid%y_cells) THEN
          c4 = (-2*ry) * chunks(c)%field%work_array7(i,j)
        ENDIF


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

        call MatSetValuesStencil(A,1,row,count-1,column,stencil,INSERT_VALUES,perr)

    ENDDO
  ENDDO

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,perr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,perr)

end subroutine setupMatA_petsc


subroutine solve_petsc(numit)

    use MPI

    INTEGER,INTENT(INOUT) :: numit
    INTEGER :: errcode, mpierr
    KSPConvergedReason :: reason

    call KSPSetOperators(kspObj,A,A,SAME_NONZERO_PATTERN,perr)
    call KSPSolve(kspObj,B,X,perr)
    call KSPGetIterationNumber(kspObj, numit, perr)

    call KSPGetConvergedReason(kspObj,reason,perr)
    if(reason < 0) then
      write(6,*) ' Error: Did not Converge. Calling MPI_Abort'
	  write(6,*) ' Divergence Reason:'

	  if(reason == -2) then
		write(6,*) ' Diverged Null'
	  else if(reason == -3) then
		write(6,*) ' Diverged ITS'
	  else if(reason == -4) then
		write(6,*) ' Diverged DTOL'
	  else if(reason == -4) then
		write(6,*) ' Diverged Breakdown'
	  else if(reason == -5) then
		write(6,*) ' Diverged Breakdown BCGS'
	  else if(reason == -6) then
		write(6,*) ' Diverged NonSymmetric'
	  else if(reason == -7) then
		write(6,*) ' Diverged Indefinte PC'
	  else if(reason == -8) then
		write(6,*) ' Diverged nanorinf'
	  else if(reason == -9) then
		write(6,*) ' Diverged indefinite mat'
	 endif

      call MPI_Abort(MPI_COMM_WORLD,errcode,mpierr)
    endif 

end subroutine solve_petsc



! Apply Paul Garrett Approach
! (1) Use CG to Execute 10 iteration, retrieve eigenvalues
! (2) Set Solver to Chebyshev type, set eigenvalues from CG execution
! (3) Set Residual for Chebyshev to reduce iteration count
! (4) Execute Chebyshev Solve
subroutine solve_petsc_pgcg(eps,max_iters,numit)

!    use MPI

#include "finclude/petscsys.h"

    INTEGER,INTENT(INOUT) :: numit
    INTEGER :: errcode, mpierr
    KSPConvergedReason :: reason
    PC :: tPC
    REAL(kind=8) :: eps
    INTEGER :: max_iters
    PetscReal :: r(0:9)    ! Real Component of EigenValue Array
    PetscReal :: c(0:9)    ! Complex Component of EigenValue Array
    PetscInt  :: neig      ! Number of Eigenvalues computed
    PetscReal :: emax      ! Max EigenValue
    PetscReal :: emin      ! Min EigenValue

    call KSPSetComputeEigenValues(kspObj,PETSC_TRUE,perr) ! Must be called before KSPSetup Routines to enable eigenvalue calculation

    call KSPSetType(kspObj,KSPCG,perr)
    call KSPGetPC(kspObj,tPC,perr)
    call PCSetType(tPC,PCNONE,perr)

    call KSPSetOperators(kspObj,A,A,SAME_NONZERO_PATTERN,perr)
    call KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,10,perr)

    call KSPSolve(kspObj,B,X,perr)

    ! Don't check convergence reasons or iteration count here - unlikely to have converged yet

    ! Calculate the EigenValues from CG Solve after 10 iterations. Array should be returned sorted.
    call KSPComputeEigenValues(kspObj,10,r,c,neig,perr)
    emax = r(neig-1)
    emin = r(0)

    if(parallel%task .eq. 0) write(6,*) ' Emax:', emax, ', Emin:', emin

    ! Main Solve is in Next section within Chebyshev

    call KSPSetComputeEigenValues(kspObj,PETSC_FALSE,perr) ! Disable EigenValue Calculation

    call KSPSetType(kspObj,KSPCHEBYSHEV,perr)
    call KSPGetPC(kspObj,tPC,perr)
    call PCSetType(tPC,PCNONE,perr)

    !call KSPSetInitialGuessNonzero(kspObj,PETSC_TRUE,perr)
    call KSPSetOperators(kspObj,A,A,SAME_NONZERO_PATTERN,perr)
    call KSPSetTolerances(kspObj,eps,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION,max_iters,perr)
    call KSPChebyshevSetEigenValues(kspObj,emax,emin,perr)

    call KSPSolve(kspObj,B,X,perr)

    call KSPGetIterationNumber(kspObj, numit, perr)

    call KSPGetConvergedReason(kspObj,reason,perr)
    if(reason < 0) then
      write(6,*) ' Error: Did not Converge. Calling MPI_Abort'
	    write(6,*) ' Divergence Reason:'

	    if(reason == -2) then
		    write(6,*) ' Diverged Null'
	      else if(reason == -3) then
		    write(6,*) ' Diverged ITS, took ', numit, ' iterations'
	      else if(reason == -4) then
		    write(6,*) ' Diverged DTOL'
	      else if(reason == -4) then
		    write(6,*) ' Diverged Breakdown'
	      else if(reason == -5) then
		    write(6,*) ' Diverged Breakdown BCGS'
	      else if(reason == -6) then
		    write(6,*) ' Diverged NonSymmetric'
	      else if(reason == -7) then
		    write(6,*) ' Diverged Indefinte PC'
	      else if(reason == -8) then
		    write(6,*) ' Diverged nanorinf'
	      else if(reason == -9) then
		    write(6,*) ' Diverged indefinite mat'
	    endif

      call MPI_Abort(MPI_COMM_WORLD,errcode,mpierr)
    endif 

end subroutine solve_petsc_pgcg


subroutine printXVec(fileName)
    use MPI

    implicit none

#   include "finclude/petscviewer.h"

    character(len=*) :: fileName
    character(len=8) :: id
    PetscViewer :: viewer

    ! Write Global Solution Vector
    call PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    call PetscViewerSetType(viewer,"ascii",perr)
    call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX,perr)
    call PetscViewerFileSetName(viewer,fileName,perr)
    call VecView(X,viewer,perr)
    call PetscViewerDestroy(viewer,perr)

end subroutine printXVec



subroutine printBVec(fileName)
    use MPI

    implicit none

#   include "finclude/petscviewer.h"

    character(len=*) :: fileName
    character(len=8) :: id
    PetscViewer :: viewer

    call PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    call PetscViewerSetType(viewer,"ascii",perr)
    call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_INDEX,perr)
    call PetscViewerFileSetName(viewer,fileName,perr)
    call VecView(B,viewer,perr)
    call PetscViewerDestroy(viewer,perr)

end subroutine printBVec

subroutine printMatA(fileName)

    use MPI

    implicit none

#   include "finclude/petscviewer.h"

    character(len=*) :: fileName
    character(len=8) :: id
    PetscViewer :: viewer

    call PetscViewerCreate(MPI_COMM_WORLD,viewer,perr)
    call PetscViewerSetType(viewer,"ascii",perr)
    call PetscViewerSetFormat(viewer,PETSC_VIEWER_DEFAULT,perr)
    call PetscViewerFileSetName(viewer,fileName,perr)
    call MatView(A,viewer,perr)
    call PetscViewerDestroy(viewer,perr)

end subroutine printMatA

end module PETScTeaLeaf
