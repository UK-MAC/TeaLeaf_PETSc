MODULE PETScTeaLeaf

  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE setup_petsc(eps,max_iters)


  REAL(kind=8) :: eps
  INTEGER :: max_iters

END SUBROUTINE setup_petsc

SUBROUTINE cleanup_petsc()

END SUBROUTINE cleanup_petsc

SUBROUTINE setupSol_petsc(c,rx,ry)

    INTEGER       :: c                                ! What chunk are we solving
    REAL(KIND=8),INTENT(IN)   :: rx,ry

END SUBROUTINE setupSol_petsc

SUBROUTINE setupRHS_petsc(c,rx,ry)

    INTEGER       :: c                                ! What chunk are we solving
    REAL(KIND=8),INTENT(IN)   :: rx,ry

END SUBROUTINE setupRHS_petsc

SUBROUTINE getSolution_petsc(c)

    INTEGER       :: c                                ! What chunk are we solving

END SUBROUTINE getSolution_petsc

SUBROUTINE setupMatA_petsc(c,rx,ry)

  INTEGER       :: c                                ! What chunk are we solving
  REAL(KIND=8),INTENT(IN)   :: rx,ry

END SUBROUTINE setupMatA_petsc

SUBROUTINE solve_petsc(numit,error)

    INTEGER,INTENT(INOUT) :: numit
    REAL(KIND=8),INTENT(INOUT) :: error

END SUBROUTINE solve_petsc

SUBROUTINE solve_petsc_pgcg(eps,max_iters,numit_cg,numit_cheby,error)


    INTEGER,INTENT(INOUT) :: numit_cg, numit_cheby
    INTEGER :: errcode, mpierr
    REAL(kind=8) :: eps,error
    INTEGER :: max_iters

END SUBROUTINE solve_petsc_pgcg

SUBROUTINE printXVec(fileName)

    IMPLICIT NONE

    CHARACTER(LEN=*) :: fileName

END SUBROUTINE printXVec

SUBROUTINE printBVec(fileName)

    IMPLICIT NONE

    CHARACTER(LEN=*) :: fileName

END SUBROUTINE printBVec

SUBROUTINE printMatA(fileName)

    IMPLICIT NONE


    CHARACTER(LEN=*) :: fileName

END SUBROUTINE printMatA

END MODULE PETScTeaLeaf
