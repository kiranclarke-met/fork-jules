MODULE model_runner_mod
IMPLICIT NONE
CONTAINS
    SUBROUTINE phi_mh_model_runner(input, output)
        USE phi_mh_model_mod, ONLY: phi_mh_model
        USE um_types, ONLY: real_jlslsm
        IMPLICIT NONE
        
        REAL(KIND=real_jlslsm), INTENT(IN) :: input(6)
        REAL(KIND=real_jlslsm), INTENT(OUT) :: output(1)

        REAL(KIND=4) :: input4(1,6)
        REAL(KIND=4) :: output4(1,1)
        REAL(KIND=4) :: output_c4(1)

        ! AU trained version (2026-11-04_10:12:16_49310581)
        ! ['Tair', 'LWdown', 'SWdown', 'height', 't_surf', 'Wind']
        REAL(KIND=4), PARAMETER :: means_x(6) = (/ &
            296.7774522931202, 369.1881082483503, 376.2386889028187, &
            19.49643008174675, 299.8030091373082, 2.5449246299614416 /)
        REAL(KIND=4), PARAMETER :: scales_x(6) = (/ &
            7.504076025198663, 51.97930654914925, 336.3099266207799, &
            15.860279753622542, 9.708808361873714, 1.429862896548756 /)
        
        ! ['stability_momentum']
        REAL(KIND=4), PARAMETER :: means_y(1) = (/ &
            0.4537166975220352 /)
        REAL(KIND=4), PARAMETER :: scales_y(1) = (/ &
            3.1674286203689337 /)
        
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
        REAL(KIND=real_jlslsm) :: input(6)
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