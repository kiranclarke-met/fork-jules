! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE phi_m_h_ml_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHI_M_H_ML_MOD'

CONTAINS
SUBROUTINE phi_m_h (                                                           &
 points,surft_pts,surft_index,                                                 &
 tl_1,lw_down,sw_surft,z_tq,tstar,vshr,z0m,z0h,phi_m,phi_h                      &
)

USE model_runner_mod, ONLY: phi_mh_model_runner

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER ::                                                                     &
 points                                                                        &
                      ! IN Number of points.
,surft_pts                                                                     &
                      ! IN Number of tile points.
,surft_index(points)  ! IN Index of tile points.

REAL(KIND=real_jlslsm) ::                                                      &
 tl_1(points)                                                                  &
                  ! IN Air temperature at the lowest model level (K).
,lw_down(points)                                                               &
                  ! IN Downward longwave radiation (W/m2).
,sw_surft(points)                                                              &
                  ! IN Downward shortwave radiation on tiles (W/m2).
,z_tq(points)                                                                  &
                  ! IN Height of temperature, moisture and scalar
                  !    level above the roughness height (m).
,tstar(points)                                                                 &
                  ! IN Surface temperature (K).
,vshr(points)                                                                   &
                  ! IN Wind speed difference between the surface and
!                       !    the lowest wind level in the atmosphere (m/s).
,z0m(points)                                                                    &
                  ! IN Roughness length for momentum (m).
,z0h(points)    ! IN Roughness length for heat/moisture/scalars (m).

REAL(KIND=real_jlslsm) ::                                                      &
 phi_m(points)                                                                 &
                  ! OUT Stability function for momentum.
,phi_h(points)  ! OUT Stability function for
!                       !     heat/moisture/scalars.

!    Workspace usage----------------------------------------------------
!    No work areas are required.

!  Define local variables.

INTEGER :: k,l   ! Loop counter; tile point index.

REAL(KIND=real_jlslsm) :: input(6)   ! Inputs to the neural network.
REAL(KIND=real_jlslsm) :: output(2)  ! Stability corrections from the neural network.

REAL(KIND=real_jlslsm) ::                                                      &
 phi_mn                                                                        &
                ! Neutral value of stability function for momentum
,phi_hn       ! Neutral value of stability function for scalars.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHI_M_H'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                               &
!$OMP  PRIVATE(k,l,input,output,phi_mn,phi_hn)                                 &
!$OMP  SHARED(surft_pts,surft_index,tl_1,lw_down,sw_surft,z_tq,tstar,vshr,     &
!$OMP  z0m,z0h,phi_m,phi_h) IF(surft_pts > 1)
DO k = 1,surft_pts
  l = surft_index(k)

  !-----------------------------------------------------------------------
  ! 1. Calculate neutral values of PHI_M and PHI_H.
  !-----------------------------------------------------------------------

  phi_mn = LOG( (z_tq(l) + z0m(l)) / z0m(l) )
  phi_hn = LOG( (z_tq(l) + z0m(l)) / z0h(l) )

  input(1) = lw_down(l)
  input(2) = sw_surft(l)
  input(3) = vshr(l)
  input(4) = tl_1(l)
  input(5) = tstar(l)
  input(6) = z_tq(l)

  CALL phi_mh_model_runner(input, output)

  phi_m(l) = phi_mn + output(1)
  phi_h(l) = phi_hn + output(2)

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE phi_m_h
END MODULE phi_m_h_ml_mod
