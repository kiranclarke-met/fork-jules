! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE phi_m_h_nn_debug_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CONTAINS

SUBROUTINE phi_m_h_nn_debug(                                                   &
 points,surft_pts,surft_index,pts_index,                                       &
 tl_1,lw_down,sw_surft,vshr_land,z1_tq,z0m                                     &
)

USE ennuf_phi_m_model_runner_mod, ONLY: ennuf_phi_m_model_runner

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 points,surft_pts,surft_index(points),pts_index(points)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 tl_1(points),lw_down(points),sw_surft(points),vshr_land(points),              &
 z1_tq(points),z0m(points)

INTEGER :: k,l
REAL(KIND=real_jlslsm) :: input(5), phi_cap, phi_m_nn

DO k = 1,surft_pts
  l = surft_index(k)

  input(1) = tl_1(l)
  input(2) = lw_down(l)
  input(3) = sw_surft(l)
  input(4) = vshr_land(l)
  input(5) = z1_tq(l)

  phi_cap = 0.0_real_jlslsm
  IF (z0m(l) > TINY(1.0_real_jlslsm) .AND. z1_tq(l) > TINY(1.0_real_jlslsm)) THEN
    phi_cap = LOG(z1_tq(l) / z0m(l)) - 0.000001_real_jlslsm
  END IF

  CALL ennuf_phi_m_model_runner(input, phi_cap, phi_m_nn)

  WRITE(6,'(A,I0,A,I0,A,ES14.6)') 'phi_m_nn_debug l=',l,                      &
       ' pts_index=',pts_index(l),' value=',phi_m_nn
END DO

END SUBROUTINE phi_m_h_nn_debug

END MODULE phi_m_h_nn_debug_mod
