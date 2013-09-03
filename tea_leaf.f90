MODULE tea_leaf_module

CONTAINS

SUBROUTINE tea_leaf()
 
  USE clover_module
  USE tea_leaf_kernel_module
  USE update_halo_module

  IMPLICIT NONE

  INTEGER :: c, n, j,k
  REAL(KIND=8) :: ry,rx, error

  INTEGER :: fields(NUM_FIELDS)
  INTEGER :: numit,numit_cg,numit_cheby                ! number of iterations taken
  CHARACTER(len=30) :: n_char

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      fields=0
      fields(FIELD_ENERGY1) = 1
      fields(FIELD_DENSITY1) = 1
      CALL update_halo(fields,2)

      ! INIT - set up coefficients based on material temperature
      IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_init(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%celldx,                      &
              chunks(c)%field%celldy,                      &
              chunks(c)%field%volume,                      &
              chunks(c)%field%density1,                    &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%work_array1,                 &
              chunks(c)%field%u,                           &
              chunks(c)%field%work_array2,                 &
              chunks(c)%field%work_array3,                 &
              chunks(c)%field%work_array4,                 &
              chunks(c)%field%work_array5,                 &
              chunks(c)%field%work_array6,                 &
              chunks(c)%field%work_array7,                 &
              coefficient)
      ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_init_c(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%celldx,                      &
              chunks(c)%field%celldy,                      &
              chunks(c)%field%volume,                      &
              chunks(c)%field%density1,                    &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%work_array1,                 &
              chunks(c)%field%u,                           &
              chunks(c)%field%work_array2,                 &
              chunks(c)%field%work_array3,                 &
              chunks(c)%field%work_array4,                 &
              chunks(c)%field%work_array5,                 &
              chunks(c)%field%work_array6,                 &
              chunks(c)%field%work_array7,                 &
              coefficient)
      ENDIF


      ! JACOBI

      rx = dt/(chunks(c)%field%celldx(chunks(c)%field%x_min)**2);
      ry = dt/(chunks(c)%field%celldy(chunks(c)%field%y_min)**2);

      IF(.not. use_PETSC_kernels) THEN
        DO n=1,max_iters

          IF(use_fortran_kernels) THEN
              CALL tea_leaf_kernel_solve(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  rx,                                          &
                  ry,                                          &
                  chunks(c)%field%work_array6,                 &
                  chunks(c)%field%work_array7,                 &
                  error,                                       &
                  chunks(c)%field%work_array1,                 &
                  chunks(c)%field%u,                           &
                  chunks(c)%field%work_array2)
          ELSEIF(use_C_kernels) THEN
              CALL tea_leaf_kernel_solve_c(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  rx,                                          &
                  ry,                                          &
                  chunks(c)%field%work_array6,                 &
                  chunks(c)%field%work_array7,                 &
                  error,                                       &
                  chunks(c)%field%work_array1,                 &
                  chunks(c)%field%u,                           &
                  chunks(c)%field%work_array2)
          ENDIF

          ! CALL update_halo
          fields=0
          fields(FIELD_U) = 1
          CALL update_halo(fields,2)

          CALL clover_max(error)

          IF (error .LT. eps) EXIT

        ENDDO
      ELSEIF (use_PETSC_kernels) THEN

        ! Substitute for PETSc Solve

		    ! write(6,*) 'Calling PETSc Object Population and Solve'

        CALL setupMatA_petsc(c,rx,ry)
        CALL setupRHS_petsc(c,rx,ry)
        CALL setupSol_petsc(c,rx,ry)

        write(n_char,'(I30)'),step

        !CALL printXVec('XVec.' // TRIM(adjustl(n_char)) // '.iter')
        !CALL printBVec('BVec.' // TRIM(adjustl(n_char)) // '.iter')
        !if(step .eq. 1) CALL printMatA('MatA.' // TRIM(adjustl(n_char)) // '.iter')

        
        if(use_pgcg) then    
          !if(parallel%task .eq. 0) write(0,*) ' Using PGCG'
          CALL solve_petsc_pgcg(eps,max_iters,numit_cg,numit_cheby)  ! Use Paul Garrett's Approach
          if(parallel%task .eq. 0) write(6,*) 'Achieved convergence in ', numit_cg ,' CG iterations and ', numit_cheby, ' Cheby Iterations'
          if(parallel%task .eq. 0) write(6,*) 'Current Total Iterations is : ',  total_cg_iter, ' CG Iterations and ', total_cheby_iter, ' Chebyshev Iterations'

        else 
          CALL solve_petsc(numit)    ! Use Command Line Specified Approach
          if(parallel%task .eq. 0) write(6,*) 'Achieved convergence in ', numit ,' iterations'
          if(parallel%task .eq. 0) write(6,*) 'Current Total Iterations: ',  total_petsc_iter
        endif
   
        CALL getSolution_petsc(c)

        !CALL printXVec('XVec.After.' // TRIM(adjustl(n_char)) // '.iter')
      ENDIF

      ! RESET
      IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_finalise(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                           &
              chunks(c)%field%y_min,                           &
              chunks(c)%field%y_max,                           &
              chunks(c)%field%energy1,                         &
              chunks(c)%field%density1,                        &
              chunks(c)%field%u)
      ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_finalise_c(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                           &
              chunks(c)%field%y_min,                           &
              chunks(c)%field%y_max,                           &
              chunks(c)%field%energy1,                         &
              chunks(c)%field%density1,                        &
              chunks(c)%field%u)
      ENDIF

      fields=0
      fields(FIELD_ENERGY1) = 1
      CALL update_halo(fields,1)

    ENDIF

  ENDDO

END SUBROUTINE tea_leaf

END MODULE tea_leaf_module
