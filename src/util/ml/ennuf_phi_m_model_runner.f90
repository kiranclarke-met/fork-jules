MODULE ennuf_phi_m_model_runner_mod
IMPLICIT NONE
CONTAINS
    SUBROUTINE ennuf_phi_m_model_runner(input, phi_cap, output)
        ! ennuf_phi_m_model_runner
        ! estimates boundary layer stability function for momentum from air temperature, incoming longwave and shortwave, wind speed, and height
        ! standardises the input, runs the model, un-standardises, and caps the output
        !! input: real[5] - (/ tair, lwdown, swdown, wind, height /)
        !! phi_cap: real - log(z/z0)-1e-6
        !! output: real - stability function for momentum

        USE ennuf_phi_m_model_mod, ONLY: ennuf_phi_m_model ! NN model
        USE um_types, ONLY: real_jlslsm
        IMPLICIT NONE

        REAL(KIND=real_jlslsm), INTENT(IN) :: input(5)
        REAL(KIND=real_jlslsm), INTENT(IN) :: phi_cap  ! cap = log(z/z0)-1e-6 to ensure ustar > 0
        REAL(KIND=real_jlslsm), INTENT(OUT) :: output

        ! input4 etc. written by GH Copilot
        REAL(KIND=4) :: input4(1,5)
        REAL(KIND=4) :: output4(1,1)
        REAL(KIND=4) :: phi_cap4
        REAL(KIND=4) :: output_scalar4

        ! X_scaler, learned by sklearn
        REAL(KIND=4), PARAMETER :: means_x(5) = (/ &
            287.17465998290277, 325.506758326137, 220.01192477255702, &
            2.8093898848877683, 16.997653763205694 /)
        REAL(KIND=4), PARAMETER :: scales_x(5) = (/ &
            11.039222332840785, 59.861978403610635, 293.12886703490466, &
            1.8839326925133042, 16.146385206066146 /)
        
        ! y_scaler, learned by sklearn
        REAL(KIND=4), PARAMETER :: mean_y = -0.2011566915146441
        REAL(KIND=4), PARAMETER :: scale_y = 1.0931165851019948

        input4(1,:) = REAL(input(:), KIND=4)
        phi_cap4 = REAL(phi_cap, KIND=4)

        input4(1,:) = (input4(1,:) - means_x) / scales_x  ! standardise array

        CALL ennuf_phi_m_model(input4, output4)

        output_scalar4 = output4(1,1)
        output_scalar4 = (output_scalar4 * scale_y) + mean_y  ! un-standardise output
        output_scalar4 = SIGN(EXP(ABS(output_scalar4)) - 1.0, output_scalar4)  ! un-transform output
        output_scalar4 = MIN(output_scalar4, phi_cap4)  ! cap output

        output = REAL(output_scalar4, KIND=real_jlslsm)
        
    END SUBROUTINE ennuf_phi_m_model_runner

    SUBROUTINE ennuf_phi_m_model_runner_debug()
        USE um_types, ONLY: real_jlslsm
        IMPLICIT NONE

        INTEGER :: unit = 1234
        REAL(KIND=real_jlslsm) :: input(5) ! array for input data
        REAL(KIND=real_jlslsm) :: output ! scalar output data

        REAL(KIND=real_jlslsm), PARAMETER :: phi_cap = 4.0_real_jlslsm

        ! read in input data from file to input array 
        OPEN(unit, FILE="data/input_sample.dat",FORM="UNFORMATTED", STATUS="OLD", ACTION="READ", ACCESS="STREAM")
            READ(unit) input
        CLOSE(unit)

        ! call NN model 
        CALL ennuf_phi_m_model_runner(input, phi_cap, output)

        ! write output data to file from output array 
        OPEN(unit, FILE="data/output_prediction_fortran.dat",FORM="UNFORMATTED", STATUS="REPLACE", ACTION="WRITE", ACCESS="STREAM")
            WRITE(unit) output
        CLOSE(unit)

    END SUBROUTINE ennuf_phi_m_model_runner_debug
END MODULE ennuf_phi_m_model_runner_mod