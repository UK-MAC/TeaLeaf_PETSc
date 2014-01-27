MODULE generate_chunk_kernel_module

CONTAINS

SUBROUTINE generate_chunk_kernel(x_min,x_max,y_min,y_max, &
                                 vertexx,                 &
                                 vertexy,                 &
                                 cellx,                   &
                                 celly,                   &
                                 density0,                &
                                 energy0,                 &
                                 xvel0,                   &
                                 yvel0,                   &
                                 u0,                      &
                                 number_of_states,        &
                                 state_density,           &
                                 state_energy,            &
                                 state_xvel,              &
                                 state_yvel,              &
                                 state_xmin,              &
                                 state_xmax,              &
                                 state_ymin,              &
                                 state_ymax,              &
                                 state_radius,            &
                                 state_geometry,          &
                                 g_rect,                  &
                                 g_circ                   )

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexy
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: cellx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celly
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u0
  INTEGER      :: number_of_states
  REAL(KIND=8), DIMENSION(number_of_states) :: state_density
  REAL(KIND=8), DIMENSION(number_of_states) :: state_energy
  REAL(KIND=8), DIMENSION(number_of_states) :: state_xvel
  REAL(KIND=8), DIMENSION(number_of_states) :: state_yvel
  REAL(KIND=8), DIMENSION(number_of_states) :: state_xmin
  REAL(KIND=8), DIMENSION(number_of_states) :: state_xmax
  REAL(KIND=8), DIMENSION(number_of_states) :: state_ymin
  REAL(KIND=8), DIMENSION(number_of_states) :: state_ymax
  REAL(KIND=8), DIMENSION(number_of_states) :: state_radius
  INTEGER     , DIMENSION(number_of_states) :: state_geometry
  INTEGER      :: g_rect
  INTEGER      :: g_circ

  REAL(KIND=8) :: radius
  INTEGER      :: state

  INTEGER      :: j,k,jt,kt

  ! State 1 is always the background state

!$OMP PARALLEL
!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
      energy0(j,k)=state_energy(1)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
      density0(j,k)=state_density(1)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
      xvel0(j,k)=state_xvel(1)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
      yvel0(j,k)=state_yvel(1)
    ENDDO
  ENDDO
!$OMP END DO

  DO state=2,number_of_states

! Could the velocity setting be thread unsafe?

!$OMP DO PRIVATE(radius)
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        IF(state_geometry(state).EQ.g_rect ) THEN
          IF(vertexx(j).GE.state_xmin(state).AND.vertexx(j).LT.state_xmax(state)) THEN
            IF(vertexy(k).GE.state_ymin(state).AND.vertexy(k).LT.state_ymax(state)) THEN
              energy0(j,k)=state_energy(state)
              density0(j,k)=state_density(state)
              DO kt=k,k+1
                DO jt=j,j+1
                  xvel0(jt,kt)=state_xvel(state)
                  yvel0(jt,kt)=state_yvel(state)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ELSEIF(state_geometry(state).EQ.g_circ ) THEN
          radius=SQRT(cellx(j)*cellx(j)+celly(k)*celly(k))
          IF(radius.LE.state_radius(state))THEN
            energy0(j,k)=state_energy(state)
            density0(j,k)=state_density(state)
            DO kt=k,k+1
              DO jt=j,j+1
                xvel0(jt,kt)=state_xvel(state)
                yvel0(jt,kt)=state_yvel(state)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO

  ENDDO

!$OMP DO 
  DO k=y_min-1, y_max+1
    DO j=x_min-1, x_max+1
      u0(j,k) =  energy0(j,k) * density0(j,k)
    ENDDO
  ENDDO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE generate_chunk_kernel

END MODULE generate_chunk_kernel_module
