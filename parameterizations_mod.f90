! module parameterizations_mod
! PURPOSE: Compute parametrizations (at this stage, chemistry emissions and depositions only)
module parameterizations_mod
    use radiation_mod
    use chemistry_mod
    use prognostics_mod
    use time_mod
    implicit none
    ! TODO: Tallenna parametrisointitendenssit
contains

    subroutine parameterizations_init()
        implicit none

    end subroutine

    ! subroutine compute_parameterizations(progn)
    ! PURPOSE: Main driver for parameterization computations, called in the main time loop.
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_parameterizations(progn)
        use time_mod
        implicit none
        type(prognostics_type), intent(inout) :: progn

        ! Compute chemistry parameterizations (emission and deposition) every dt_chem.
        IF ( time >= time_start_chemistry .and. MOD( NINT((time - time_start)*10.0), NINT(dt_chem*10.0)) == 0 ) THEN
            call compute_chemistry(progn, -1.0_dp) ! A value < 0 indicates that the PAR be computed using the radiation routine.
            !print *, 'Chemistry called'
        end if

    end subroutine

    subroutine compute_chemistry(progn, PAR_val)
        implicit none
        type(prognostics_type), intent(inout) :: progn
        real(kind = 8), intent(in) :: PAR_val ! Predefined PAR value for testing. If PAR_val < 0, the value is computed using a routine.
        integer :: level, chemical
        real(kind = 8) :: concentrations(num_chemical_elements), concentration_tendencies(num_chemical_elements), t
        integer :: input_array(2), istate = 1

        if (PAR_val < 0) then
            call compute_radiation() ! Compute PAR
        else
            PAR = PAR_val
        end if

        call compute_aerodynamic_resistance(progn) ! r_a

        input_array(1) = num_chemical_elements
        !print*, time
        do level = 2, nz-1 ! TESTAUS
            istate = 1
            input_array(2) = level
            concentration_tendencies = 0

            concentrations(1) = progn%O3%concentration(level)
            concentrations(2) = progn%O1D%concentration(level)
            concentrations(3) = progn%OH%concentration(level)
            concentrations(4) = progn%REST%concentration(level)
            concentrations(5) = progn%NO2%concentration(level)
            concentrations(6) = progn%NO%concentration(level)
            concentrations(7) = progn%CH2O%concentration(level)
            concentrations(8) = progn%HO2%concentration(level)
            concentrations(9) = progn%CO%concentration(level)
            concentrations(10) = progn%CO2%concentration(level)
            concentrations(11) = progn%CH4%concentration(level)
            concentrations(12) = progn%CH3O2%concentration(level)
            concentrations(13) = progn%isoprene%concentration(level)
            concentrations(14) = progn%RO2%concentration(level)
            concentrations(15) = progn%MVK%concentration(level)
            concentrations(16) = progn%H2O2%concentration(level)
            concentrations(17) = progn%HNO3%concentration(level)
            concentrations(18) = progn%NO3%concentration(level)
            concentrations(19) = progn%N2O5%concentration(level)
            concentrations(20) = progn%SO2%concentration(level)
            concentrations(21) = progn%H2SO4_P%concentration(level)
            concentrations(22) = progn%H2SO4%concentration(level)
            concentrations(23) = progn%alpha_pinene%concentration(level)
            concentrations(24) = progn%HNO3_P%concentration(level)
            concentrations(25) = progn%ELVOC%concentration(level)

            t = time
            call dlsode(f, input_array, concentrations, t, t+dt_chem, itol, rtol, atol, itask, &
                        istate, iopt, rwork, lrw, iwork, liw, jac, mf)
            !call f(input_array, time, concentrations, concentration_tendencies)
            !print *, concentration_tendencies
            !stop

            progn%O3%concentration(level) = concentrations(1)
            progn%O1D%concentration(level) = concentrations(2)
            progn%OH%concentration(level) = concentrations(3)
            progn%REST%concentration(level) = concentrations(4)
            progn%NO2%concentration(level) = concentrations(5)
            progn%NO%concentration(level) = concentrations(6)
            progn%CH2O%concentration(level) = concentrations(7)
            progn%HO2%concentration(level) = concentrations(8)
            progn%CO%concentration(level) = concentrations(9)
            progn%CO2%concentration(level) = concentrations(10)
            progn%CH4%concentration(level) = concentrations(11)
            progn%CH3O2%concentration(level) = concentrations(12)
            progn%isoprene%concentration(level) = concentrations(13)
            progn%RO2%concentration(level) = concentrations(14)
            progn%MVK%concentration(level) = concentrations(15)
            progn%H2O2%concentration(level) = concentrations(16)
            progn%HNO3%concentration(level) = concentrations(17)
            progn%NO3%concentration(level) = concentrations(18)
            progn%N2O5%concentration(level) = concentrations(19)
            progn%SO2%concentration(level) = concentrations(20)
            progn%H2SO4_P%concentration(level) = concentrations(21)
            progn%H2SO4%concentration(level) = concentrations(22)
            progn%alpha_pinene%concentration(level) = concentrations(23)
            progn%HNO3_P%concentration(level) = concentrations(24)
            progn%ELVOC%concentration(level) = concentrations(25)

        end do
        !print*, time
        !print *, concentrations
        !print *, get_exp_coszen(time, daynumber, latitude)


    contains


        SUBROUTINE  jac (NEQ, T, Y, ML, MU, PD, NRPD)
        INTEGER  NEQ, ML, MU, NRPD
        DOUBLE PRECISION  T, Y(3), PD(NRPD,3)
        RETURN
        END subroutine

        subroutine f(iarray, t, y, ydot)
            use parameters_mod
            implicit none
            double precision t, y(*), ydot(*)
            integer :: iarray(2), level

            level = iarray(2)
            !print *, level

            progn%O3%concentration(level) = 24E-9 * progn%M(level)
            progn%O1D%concentration(level) = y(2)
            progn%OH%concentration(level) = y(3)
            progn%REST%concentration(level) = y(4)
            progn%NO2%concentration(level) = 0.2E-9 * progn%M(level)
            progn%NO%concentration(level) = 0.07E-9 * progn%M(level)
            progn%CH2O%concentration(level) = y(7)
            progn%HO2%concentration(level) = y(8)
            progn%CO%concentration(level) = 100E-9 * progn%M(level)
            progn%CO2%concentration(level) = y(10)
            progn%CH4%concentration(level) = 1759E-9 * progn%M(level)
            progn%CH3O2%concentration(level) = y(12)
            if (box) then
                progn%isoprene%concentration(level) = 2.2E-9 * progn%M(level)
            else
                progn%isoprene%concentration(level) = y(13)
            end if
            progn%RO2%concentration(level) = y(14)
            progn%MVK%concentration(level) = y(15)
            progn%H2O2%concentration(level) = y(16)
            progn%HNO3%concentration(level) = y(17)
            progn%NO3%concentration(level) = y(18)
            progn%N2O5%concentration(level) = y(19)
            progn%SO2%concentration(level) = 0.5E-9 *progn%M(level)
            progn%H2SO4_P%concentration(level) = y(21)
            progn%H2SO4%concentration(level) = y(22)
            if (box) then
                progn%alpha_pinene%concentration(level) = 2.2E-9 * progn%M(level)
            else
                progn%alpha_pinene%concentration(level) = y(23)
            end if
            progn%HNO3_P%concentration(level) = y(24)
            progn%ELVOC%concentration(level) = y(25)

            call compute_reactions(progn, level)

            call progn%O3%compute_parameterized_tendency(progn, level)
            call progn%O1D%compute_parameterized_tendency(progn, level)
            call progn%OH%compute_parameterized_tendency(progn, level)
            call progn%REST%compute_parameterized_tendency(progn, level)
            call progn%NO2%compute_parameterized_tendency(progn, level)
            call progn%NO%compute_parameterized_tendency(progn, level)
            call progn%CH2O%compute_parameterized_tendency(progn, level)
            call progn%HO2%compute_parameterized_tendency(progn, level)
            call progn%CO%compute_parameterized_tendency(progn, level)
            call progn%CO2%compute_parameterized_tendency(progn, level)
            call progn%CH4%compute_parameterized_tendency(progn, level)
            call progn%CH3O2%compute_parameterized_tendency(progn, level)
            call progn%isoprene%compute_parameterized_tendency(progn, level)
            call progn%RO2%compute_parameterized_tendency(progn, level)
            call progn%MVK%compute_parameterized_tendency(progn, level)
            call progn%H2O2%compute_parameterized_tendency(progn, level)
            call progn%HNO3%compute_parameterized_tendency(progn, level)
            call progn%NO3%compute_parameterized_tendency(progn, level)
            call progn%N2O5%compute_parameterized_tendency(progn, level)
            call progn%SO2%compute_parameterized_tendency(progn, level)
            call progn%H2SO4_P%compute_parameterized_tendency(progn, level)
            call progn%H2SO4%compute_parameterized_tendency(progn, level)
            call progn%alpha_pinene%compute_parameterized_tendency(progn, level)
            call progn%HNO3_P%compute_parameterized_tendency(progn, level)
            call progn%ELVOC%compute_parameterized_tendency(progn, level)

            ydot(1) = 0
            ydot(2) = progn%O1D%parameterized_tendency(level)
            ydot(3) = progn%OH%parameterized_tendency(level)
            ydot(4) = 0 !progn%REST%parameterized_tendency(level)
            ydot(5) = 0
            ydot(6) = 0
            ydot(7) = progn%CH2O%parameterized_tendency(level)
            ydot(8) = progn%HO2%parameterized_tendency(level)
            ydot(9) = 0
            ydot(10) = 0 !progn%CO2%parameterized_tendency(level)
            ydot(11) = 0
            ydot(12) = progn%CH3O2%parameterized_tendency(level)
            if (box) then
                ydot(13) = 0
            else
                ydot(13) = progn%isoprene%parameterized_tendency(level)
            end if
            ydot(14) = progn%RO2%parameterized_tendency(level)
            ydot(15) = progn%MVK%parameterized_tendency(level)
            ydot(16) = progn%H2O2%parameterized_tendency(level)
            ydot(17) = progn%HNO3%parameterized_tendency(level)
            ydot(18) = progn%NO3%parameterized_tendency(level)
            ydot(19) = progn%N2O5%parameterized_tendency(level)
            ydot(20) = 0
            ydot(21) = progn%H2SO4_P%parameterized_tendency(level)
            ydot(22) = progn%H2SO4%parameterized_tendency(level)
            if (box) then
                ydot(23) = 0
            else
                ydot(23) = progn%alpha_pinene%parameterized_tendency(level)
            end if
            ydot(24) = progn%HNO3_P%parameterized_tendency(level)
            ydot(25) = progn%ELVOC%parameterized_tendency(level)

        end subroutine


    end subroutine

    subroutine compute_reactions(progn, level)
        implicit none
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        ! TODO: Naitahan ei tarvits laskea joka f:n iteraatiolla!
        call O3_diss%compute_rate_coefficient(progn, level)
        call O1D_H2O%compute_rate_coefficient(progn, level)
        call O1D_N2%compute_rate_coefficient(progn, level)
        call O1D_O2%compute_rate_coefficient(progn, level)
        call NO2_diss%compute_rate_coefficient(progn, level)
        call CH2O_DISS%compute_rate_coefficient(progn, level)
        call OH_CO%compute_rate_coefficient(progn, level)
        call OH_CH4%compute_rate_coefficient(progn, level)
        call OH_isoprene%compute_rate_coefficient(progn, level)
        call OH_MVK%compute_rate_coefficient(progn, level)
        call HO2_NO%compute_rate_coefficient(progn, level)
        call CH3O2_NO%compute_rate_coefficient(progn, level)
        call RO2_NO%compute_rate_coefficient(progn, level)
        call OH_CH2O%compute_rate_coefficient(progn, level)
        call HO2_diss%compute_rate_coefficient(progn, level)
        call CH3O2_HO2%compute_rate_coefficient(progn, level)
        call RO2_HO2%compute_rate_coefficient(progn, level)
        call OH_NO2%compute_rate_coefficient(progn, level)
        call NO_O3%compute_rate_coefficient(progn, level)
        call OH_HO2%compute_rate_coefficient(progn, level)
        call OH_H2O2%compute_rate_coefficient(progn, level)
        call NO_NO3%compute_rate_coefficient(progn, level)
        call NO2_O3%compute_rate_coefficient(progn, level)
        call NO2_NO3%compute_rate_coefficient(progn, level)
        call N2O5_diss%compute_rate_coefficient(progn, level)
        call N2O5_H2O%compute_rate_coefficient(progn, level)
        call N2O5_2H2O%compute_rate_coefficient(progn, level)
        call HO2_O3%compute_rate_coefficient(progn, level)
        call OH_SO2%compute_rate_coefficient(progn, level)
        call H2SO4_diss%compute_rate_coefficient(progn, level)
        call HNO3_diss%compute_rate_coefficient(progn, level)
        call alpha_pinene_OH%compute_rate_coefficient(progn, level)
        call alpha_pinene_O3%compute_rate_coefficient(progn, level)
        call isoprene_O3%compute_rate_coefficient(progn, level)
        call ELVOC_form%compute_rate_coefficient(progn, level)

        call O3_diss%compute_rate(progn, level)
        call O1D_H2O%compute_rate(progn, level)
        call O1D_N2%compute_rate(progn, level)
        call O1D_O2%compute_rate(progn, level)
        call NO2_diss%compute_rate(progn, level)
        call CH2O_DISS%compute_rate(progn, level)
        call OH_CO%compute_rate(progn, level)
        call OH_CH4%compute_rate(progn, level)
        call OH_isoprene%compute_rate(progn, level)
        call OH_MVK%compute_rate(progn, level)
        call HO2_NO%compute_rate(progn, level)
        call CH3O2_NO%compute_rate(progn, level)
        call RO2_NO%compute_rate(progn, level)
        call OH_CH2O%compute_rate(progn, level)
        call HO2_diss%compute_rate(progn, level)
        call CH3O2_HO2%compute_rate(progn, level)
        call RO2_HO2%compute_rate(progn, level)
        call OH_NO2%compute_rate(progn, level)
        call NO_O3%compute_rate(progn, level)
        call OH_HO2%compute_rate(progn, level)
        call OH_H2O2%compute_rate(progn, level)
        call NO_NO3%compute_rate(progn, level)
        call NO2_O3%compute_rate(progn, level)
        call NO2_NO3%compute_rate(progn, level)
        call N2O5_diss%compute_rate(progn, level)
        call N2O5_H2O%compute_rate(progn, level)
        call N2O5_2H2O%compute_rate(progn, level)
        call HO2_O3%compute_rate(progn, level)
        call OH_SO2%compute_rate(progn, level)
        call H2SO4_diss%compute_rate(progn, level)
        call HNO3_diss%compute_rate(progn, level)
        call alpha_pinene_OH%compute_rate(progn, level)
        call alpha_pinene_O3%compute_rate(progn, level)
        call isoprene_O3%compute_rate(progn, level)
        call ELVOC_form%compute_rate(progn, level)

    end subroutine


    subroutine parameterizations_output(progn)
        implicit none
        type(prognostics_type), intent(in) :: progn

    end subroutine

end module parameterizations_mod
