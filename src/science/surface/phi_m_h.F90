! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE phi_m_h_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHI_M_H_MOD'

CONTAINS
!   SUBROUTINE PHI_M_H ----------------------------------------------

!  Purpose: Calculate the integrated froms of the Monin-Obukhov
!           stability functions for surface exchanges.

!  Documentation: UM Documentation Paper No 24.

!--- Arguments:---------------------------------------------------------
SUBROUTINE phi_m_h (                                                           &
 points,surft_pts,surft_index,pts_index,                                       &
 tl_1,lw_down,sw_surft,vshr_land,z1_tq,                                        &
 recip_l_mo,z_uv,z_tq,z0m,z0h,phi_m,phi_h                                      &
)

USE atm_fields_bounds_mod, ONLY: tdims
USE theta_field_sizes, ONLY: t_i_length

USE jules_surface_mod, ONLY: a,b,d,c_over_d

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE ennuf_phi_m_model_runner_mod, ONLY: ennuf_phi_m_model_runner

IMPLICIT NONE

INTEGER ::                                                                     &
 points                                                                        &
                      ! IN Number of points.
,surft_pts                                                                     &
                      ! IN Number of tile points.
,surft_index(points)                                                           &
                      ! IN Index of tile points.
,pts_index(points)    ! IN Index of points.

REAL(KIND=real_jlslsm) ::                                                      &
tl_1(points)                                                                   &
!                       ! IN Air temperature forcing.
,lw_down(points)                                                               &
!                       ! IN Downward longwave forcing.
,sw_surft(points)                                                              &
!                       ! IN Downward shortwave forcing.
,vshr_land(points)                                                             &
!                       ! IN Wind speed forcing.
,z1_tq(points)                                                                 &
!                       ! IN Height forcing.
,recip_l_mo(points)                                                            &
!                       ! IN Reciprocal of the Monin-Obukhov length
!                       !     (m^-1).
,z_uv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
!                       ! IN Height of wind level above roughness
!                       !    height (m)
,z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     &
!                       ! IN Height of temperature, moisture and scalar
!                       !    lev above the roughness height (m).
,z0m(points)                                                                   &
                  ! IN Roughness length for momentum (m).
,z0h(points)    ! IN Roughness length for heat/moisture/scalars
!                       !    (m)

REAL(KIND=real_jlslsm) ::                                                      &
 phi_m(points)                                                                 &
                  ! OUT Stability function for momentum.
,phi_h(points)  ! OUT Stability function for
!                       !     heat/moisture/scalars.

!    Workspace usage----------------------------------------------------
!    No work areas are required.

!  Define local variables.

INTEGER :: i,j,k,l   ! Loop counter; horizontal field index.

REAL(KIND=real_jlslsm) ::                                                      &
 phi_mn                                                                        &
                ! Neutral value of stability function for momentum
,phi_hn                                                                        &
                ! Neutral value of stability function for scalars.
,zeta_uv                                                                       &
                ! Temporary in calculation of PHI_M.
,zeta_0m                                                                       &
                ! Temporary in calculation of PHI_M.
,zeta_tq                                                                       &
                ! Temporary in calculation of PHI_H.
,zeta_0h                                                                       &
                ! Temporary in calculation of PHI_H.
,x_uv_sq                                                                       &
                ! Temporary in calculation of PHI_M.
,x_0m_sq                                                                       &
                ! Temporary in calculation of PHI_M.
,x_uv                                                                          &
                ! Temporary in calculation of PHI_M.
,x_0m                                                                          &
                ! Temporary in calculation of PHI_M.
,y_tq                                                                          &
                ! Temporary in calculation of PHI_H.
,y_0h                                                                          &
                ! Temporary in calculation of PHI_H.
,phi_h_fz1                                                                     &
                ! Temporary in calculation of PHI_H.
,phi_h_fz0      ! Temporary in calculation of PHI_H.

REAL(KIND=real_jlslsm) ::                                                      &
 forcing_data(5, points),phi_cap(points)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHI_M_H'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

forcing_data(1,:) = tl_1
forcing_data(2,:) = lw_down
forcing_data(3,:) = sw_surft
forcing_data(4,:) = vshr_land
forcing_data(5,:) = z1_tq

phi_cap = LOG(z1_tq / z0m) - 0.000001

!CDIR NODEP
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                               &
!$OMP  PRIVATE(k,l,j,i,phi_mn, phi_hn,zeta_uv,zeta_tq,zeta_0m,zeta_0h,         &
!$OMP  phi_h_fz1,phi_h_fz0,x_uv_sq,x_0m_sq,x_uv,x_0m,y_tq,y_0h)                &
!$OMP  SHARED(surft_pts,surft_index,pts_index,t_i_length,z_uv,z0m,z0h,z_tq,    &
!$OMP  recip_l_mo,phi_m,phi_h,forcing_data,phi_cap) IF(surft_pts > 1)
DO k = 1,surft_pts
  l = surft_index(k)
  j=(pts_index(l) - 1) / t_i_length + 1
  i = pts_index(l) - (j-1) * t_i_length

  CALL ennuf_phi_m_model_runner(forcing_data(:,l), phi_cap(l), phi_m(l))

  !-----------------------------------------------------------------------
  ! 1. Calculate neutral values of PHI_H.
  !-----------------------------------------------------------------------

  phi_hn = LOG( (z_tq(i,j) + z0m(l)) / z0h(l) )

  !-----------------------------------------------------------------------
  ! 2. Calculate stability parameters.
  !-----------------------------------------------------------------------

  zeta_tq = (z_tq(i,j) + z0m(l)) * recip_l_mo(l)
  zeta_0h = z0h(l) * recip_l_mo(l)

  !-----------------------------------------------------------------------
  ! 3. Calculate PHI_H for neutral and stable conditions.
  !    Formulation of Beljaars and Holtslag (1991).
  !-----------------------------------------------------------------------

  IF (recip_l_mo(l)  >=  0.0) THEN
    phi_h_fz1 = SQRT(1.0 + (2.0 / 3.0) * a * zeta_tq)
    phi_h_fz0 = SQRT(1.0 + (2.0 / 3.0) * a * zeta_0h)
    phi_h(l) = phi_hn +                                                        &
                 phi_h_fz1 * phi_h_fz1 * phi_h_fz1                             &
               - phi_h_fz0 * phi_h_fz0 * phi_h_fz0                             &
               + b * ( (zeta_tq - c_over_d) * EXP(-d * zeta_tq)                &
                      -(zeta_0h - c_over_d) * EXP(-d * zeta_0h) )

    !-----------------------------------------------------------------------
    ! 4. Calculate PHI_H for unstable conditions.
    !-----------------------------------------------------------------------

  ELSE

    y_tq = SQRT(1.0-16.0 * zeta_tq)
    y_0h = SQRT(1.0-16.0 * zeta_0h)
    phi_h(l) = phi_hn - 2.0 * LOG( (1.0 + y_tq) / (1.0 + y_0h) )

  END IF

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE phi_m_h
END MODULE phi_m_h_mod
