! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE phi_m_h_ml_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PHI_M_H_ML_MOD'

CONTAINS

SUBROUTINE phi_m_h_ml (                                                       &
 points,surft_pts,surft_index,                                               &
 tl_1,lw_down,sw_down,z1_uv,tstar,wind,z0m,                                  &
 phi_m,phi_h,phi_m_corr,phi_h_corr                                          &
)

USE model_runner_mod, ONLY: phi_mh_model_runner
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: points
INTEGER, INTENT(IN) :: surft_pts
INTEGER, INTENT(IN) :: surft_index(points)

REAL(KIND=real_jlslsm), INTENT(IN) :: tl_1(points)
REAL(KIND=real_jlslsm), INTENT(IN) :: lw_down(points)
REAL(KIND=real_jlslsm), INTENT(IN) :: sw_down(points)
REAL(KIND=real_jlslsm), INTENT(IN) :: z1_uv(points)
REAL(KIND=real_jlslsm), INTENT(IN) :: tstar(points)
REAL(KIND=real_jlslsm), INTENT(IN) :: wind(points)
REAL(KIND=real_jlslsm), INTENT(IN) :: z0m(points)

REAL(KIND=real_jlslsm), INTENT(OUT) :: phi_m(points)
REAL(KIND=real_jlslsm), INTENT(OUT) :: phi_h(points)
REAL(KIND=real_jlslsm), INTENT(OUT) :: phi_m_corr(points)
REAL(KIND=real_jlslsm), INTENT(OUT) :: phi_h_corr(points)

INTEGER :: k, l
REAL(KIND=real_jlslsm) :: input(6)
REAL(KIND=real_jlslsm) :: output(1)
REAL(KIND=real_jlslsm) :: phi_mn

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PHI_M_H_ML'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

phi_m = 0.0
phi_h = 0.0
phi_m_corr = 0.0
phi_h_corr = 0.0

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                               &
!$OMP PRIVATE(k,l,input,output,phi_mn)                                  &
!$OMP SHARED(surft_pts,surft_index,tl_1,lw_down,sw_down,z1_uv,tstar,wind,z0m,  &
!$OMP        phi_m,phi_m_corr) IF(surft_pts > 1)
DO k = 1,surft_pts
  l = surft_index(k)

  input(1) = tl_1(l)
  input(2) = lw_down(l)
  input(3) = sw_down(l)
  input(4) = z1_uv(l)
  input(5) = tstar(l)
  input(6) = wind(l)

  CALL phi_mh_model_runner(input, output)

  phi_mn = LOG((z1_uv(l) + z0m(l)) / z0m(l))
  phi_m_corr(l) = output(1)
  
  phi_m(l) = phi_m_corr(l) + phi_mn
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE phi_m_h_ml

END MODULE phi_m_h_ml_mod
