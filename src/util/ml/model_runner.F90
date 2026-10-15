! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE model_runner_mod
IMPLICIT NONE
CONTAINS
    SUBROUTINE phi_mh_model_runner(input, output)
        USE phi_mh_model_mod, ONLY: phi_mh_model
        USE um_types, ONLY: real_jlslsm
        IMPLICIT NONE
        
        REAL(KIND=real_jlslsm), INTENT(IN) :: input(6)
        REAL(KIND=real_jlslsm), INTENT(OUT) :: output(2)

        REAL(KIND=4) :: input4(1,6)
        REAL(KIND=4) :: output4(1,2)
        REAL(KIND=4) :: output_c4(2)

        ! ['LWdown', 'SWdown', 'Wind', 'Tair', 't_surf', 'height']
        REAL(KIND=4), PARAMETER :: means_x(6) = (/ &
            367.83865788178264, 399.68296285347157, 2.594808933008636, &
            296.82003569854203, 300.07291348062955, 19.198102852852422 /)
        REAL(KIND=4), PARAMETER :: scales_x(6) = (/ &
            52.12088722063538, 338.8578708959273, 1.4291875098862172, &
            7.538753776866312, 9.825139431636444, 15.760030818609607 /)
        
        ! ['stability_momentum', 'stability_heat']
        REAL(KIND=4), PARAMETER :: means_y(2) = (/ &
            0.47816303684798217, 0.5976246151661571 /)
        REAL(KIND=4), PARAMETER :: scales_y(2) = (/ &
            2.9266061442735984, 10.545749676198394 /)
        
        ! cast to kind=4
        input4(1,:) = REAL(input(:), KIND=4)

        ! standardise
        input4(1,:) = (input4(1,:) - means_x) / scales_x

        CALL phi_mh_model(input4, output4)

        output_c4 = (output4(1,:) * scales_y) + means_y
        
        output = REAL(output_c4, KIND=real_jlslsm)
    END SUBROUTINE phi_mh_model_runner
END MODULE model_runner_mod