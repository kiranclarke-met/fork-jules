MODULE model_runner_mod
IMPLICIT NONE
CONTAINS
    SUBROUTINE phi_mh_model_runner(input, output)
        USE phi_mh_model_mod, ONLY: phi_mh_model
        USE um_types, ONLY: real_jlslsm
        IMPLICIT NONE
        
        REAL(KIND=real_jlslsm), INTENT(IN) :: input(5)
        REAL(KIND=real_jlslsm), INTENT(OUT) :: output(1)

        REAL(KIND=4) :: input4(1,5)
        REAL(KIND=4) :: output4(1,1)
        REAL(KIND=4) :: output_c4(1)

        ! AU trained version (2026-11-23_12:04:26_53593736)
        ! ['Tair', 'LWdown', 'SWdown', 'height', 'Wind']
        REAL(KIND=4), PARAMETER :: means_x(5) = (/ &
            296.6138234958321, 372.05691444219667, 364.82169288258785, &
            20.47087904053021, 2.5105288454760135 /)
        REAL(KIND=4), PARAMETER :: scales_x(5) = (/ &
            7.589321871379171, 52.221798489739776, 334.54418376243177, &
            16.451662441824737, 1.4398640888433654 /)
        
        ! ['stability_momentum']
        REAL(KIND=4), PARAMETER :: means_y(1) = (/ &
            0.30427611588046116 /)
        REAL(KIND=4), PARAMETER :: scales_y(1) = (/ &
            2.9615663486605337 /)
        
        ! cast to kind=4
        input4(1,:) = REAL(input(:), KIND=4)

        ! standardise
        input4(1,:) = (input4(1,:) - means_x) / scales_x

        CALL phi_mh_model(input4, output4)

        output_c4 = (output4(1,:) * scales_y) + means_y

        output = REAL(output_c4, KIND=real_jlslsm)
    END SUBROUTINE phi_mh_model_runner

    SUBROUTINE model_runner_debug()
        USE um_types, ONLY: real_jlslsm
        IMPLICIT NONE

        INTEGER :: unit = 1234
        REAL(KIND=real_jlslsm) :: input(5)
        REAL(KIND=real_jlslsm) :: output(1)

        ! read in input data from file to input array 
        OPEN(unit, FILE="ennuf/test_data/input.dat",FORM="UNFORMATTED", STATUS="OLD", ACTION="READ", ACCESS="STREAM")
            READ(unit) input
        CLOSE(unit)

        ! call NN model 
        CALL phi_mh_model_runner(input, output)

        ! write output data to file from output array 
        OPEN(unit, FILE="ennuf/test_data/output.dat",FORM="UNFORMATTED", STATUS="REPLACE", ACTION="WRITE", ACCESS="STREAM")
            WRITE(unit) output
        CLOSE(unit)
    END SUBROUTINE model_runner_debug
END MODULE model_runner_mod