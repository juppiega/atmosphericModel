! module prognostics_mod
! PURPOSE: Contain values and governing equations for all of the variables advanced forwards in time.
module prognostics_mod

    use grid_mod
    use time_mod
    use radiation_mod
    use chemistry_mod
    use aerosol_mod

    implicit none

    public
    !public prognostics_type, chemical_reaction
    !public chemical_element, compute_aerodynamic_resistance

    ! type, abstract :: chemical_element
    ! PURPOSE: Base type, from which all chemical elements (e.g. isoprene) are derived.
    !          Most of the variables and governing equations are the same, but the
    !          emission routine has to be written separately for every element.
    ! TO ADD NEW CHEMICAL ELEMENT:
    ! 1. Derive a new type from chemical_element and write reaction equations
    ! 2. Add it to prognostics_type
    ! 3. Add it to prognostics_type::init_chemical_elements
    ! 4. Add it to prognostics_type::euler_next and dynamics_mod
    ! 5. Add it to set_boundary_conditions
    ! 6. Add it to parameterizations_mod
    ! 7. Add it to write_files
    type, abstract :: chemical_element
        real(kind = 8) :: concentration(nz) ! Concentration in the atmospheric column [# / cm^3]
        real(kind = 8) :: parameterized_tendency(2:nz-1)    ! Tendency from emission, deposition, chemical reactions [# / cm^3 / s]
        real(kind = 8) :: dynamical_tendency(nz-2)    ! Tendency from atmospheric turbulence (computed in dynamics_mod) [# / cm^3 / s]
        real(kind = 8) :: reaction_rate(2:nz-1)         ! Reaction rate.
        real(kind = 8) :: molar_mass                  ! Molar mass [g / mol]
        real(kind = 8) :: emission = 0.0    ! Emission at 10 meters [# / cm^3 / s]
        real(kind = 8) :: deposition = 0.0  ! Deposition at the ground [# / cm^3 / s]
        real(kind = 8) :: Henry_const       ! Henry constant for resistance calculations [M / atm] (from chemistry_mod)
        real(kind = 8) :: reactivity        ! Reactivity (f0) for resistance calculations [unitless] (from chemistry_mod)
        integer :: file_unit
        logical :: file_created
        real(kind = 8) :: r_st, r_m, r_lu, r_dc, r_cl, r_gs, r_b, r_c ! Molecule specific resistances [s/m]. r_a in chemistry_mod
        character(len = 16) :: name
    contains
        procedure(compute_emission_interface), deferred :: compute_emission ! Emission computation routine (written in the derived class)
        procedure(compute_reaction_rate_interface), deferred :: compute_reaction_rate ! Reaction rate computation routine (written in the derived class)
        ! Resistance and deposition computations
        procedure :: compute_quasi_laminar_resistance
        procedure :: compute_canopy_resistance
        procedure :: compute_deposition_velocity
        procedure :: compute_deposition
        procedure :: output
        procedure :: close_file

        ! Tendency computations
        procedure :: compute_parameterized_tendency
        procedure :: compute_dynamical_tendency

        ! Initializator (because gfortran 5 lacks type parameters)
        procedure :: init
    end type

    ! type, extends(chemical_element) :: alpha_pinene_type
    ! PURPOSE: The alpha pinene type derived from chemical_element
    type, extends(chemical_element) :: alpha_pinene_type
    contains
        procedure :: compute_emission => compute_emission_alpha_pinene
        procedure :: compute_reaction_rate => reaction_rate_alpha_pinene
    end type

    ! type, extends(chemical_element) :: isoprene_type
    ! PURPOSE: The isoprene type derived from chemical_element
    type, extends(chemical_element) :: isoprene_type
    contains
        procedure :: compute_emission => compute_emission_isoprene
        procedure :: compute_reaction_rate => reaction_rate_isoprene
    end type

    ! type, extends(chemical_element) :: OH_type
    ! PURPOSE: The OH type derived from chemical_element
    type, extends(chemical_element) :: OH_type
    contains
        procedure :: compute_emission => compute_emission_OH
        procedure :: compute_reaction_rate => reaction_rate_OH
    end type

    ! type, extends(chemical_element) :: O3_type
    ! PURPOSE: The O3 type derived from chemical_element
    type, extends(chemical_element) :: O3_type
    contains
        procedure :: compute_emission => compute_emission_O3
        procedure :: compute_reaction_rate => reaction_rate_O3
    end type

    ! type, extends(chemical_element) :: O1D_type
    ! PURPOSE: The O1D type derived from chemical_element
    type, extends(chemical_element) :: O1D_type
    contains
        procedure :: compute_emission => compute_emission_O1D
        procedure :: compute_reaction_rate => reaction_rate_O1D
    end type

    ! type, extends(chemical_element) :: REST_type
    ! PURPOSE: The REST type derived from chemical_element
    type, extends(chemical_element) :: REST_type
    contains
        procedure :: compute_emission => compute_emission_REST
        procedure :: compute_reaction_rate => reaction_rate_REST
    end type

    ! type, extends(chemical_element) :: NO_type
    ! PURPOSE: The NO type derived from chemical_element
    type, extends(chemical_element) :: NO_type
    contains
        procedure :: compute_emission => compute_emission_NO
        procedure :: compute_reaction_rate => reaction_rate_NO
    end type

    ! type, extends(chemical_element) :: NO2_type
    ! PURPOSE: The NO2 type derived from chemical_element
    type, extends(chemical_element) :: NO2_type
    contains
        procedure :: compute_emission => compute_emission_NO2
        procedure :: compute_reaction_rate => reaction_rate_NO2
    end type

    ! type, extends(chemical_element) :: CH2O_type
    ! PURPOSE: The CH2O type derived from chemical_element
    type, extends(chemical_element) :: CH2O_type
    contains
        procedure :: compute_emission => compute_emission_CH2O
        procedure :: compute_reaction_rate => reaction_rate_CH2O
    end type

    ! type, extends(chemical_element) :: HO2_type
    ! PURPOSE: The HO2 type derived from chemical_element
    type, extends(chemical_element) :: HO2_type
    contains
        procedure :: compute_emission => compute_emission_HO2
        procedure :: compute_reaction_rate => reaction_rate_HO2
    end type

    ! type, extends(chemical_element) :: CO_type
    ! PURPOSE: The CO type derived from chemical_element
    type, extends(chemical_element) :: CO_type
    contains
        procedure :: compute_emission => compute_emission_CO
        procedure :: compute_reaction_rate => reaction_rate_CO
    end type

    ! type, extends(chemical_element) :: CO2_type
    ! PURPOSE: The CO2 type derived from chemical_element
    type, extends(chemical_element) :: CO2_type
    contains
        procedure :: compute_emission => compute_emission_CO2
        procedure :: compute_reaction_rate => reaction_rate_CO2
    end type

    ! type, extends(chemical_element) :: CH4_type
    ! PURPOSE: The CH4 type derived from chemical_element
    type, extends(chemical_element) :: CH4_type
    contains
        procedure :: compute_emission => compute_emission_CH4
        procedure :: compute_reaction_rate => reaction_rate_CH4
    end type

    ! type, extends(chemical_element) :: CH3O2_type
    ! PURPOSE: The CH3O2 type derived from chemical_element
    type, extends(chemical_element) :: CH3O2_type
    contains
        procedure :: compute_emission => compute_emission_CH3O2
        procedure :: compute_reaction_rate => reaction_rate_CH3O2
    end type

    ! type, extends(chemical_element) :: RO2_type
    ! PURPOSE: The RO2 type derived from chemical_element
    type, extends(chemical_element) :: RO2_type
    contains
        procedure :: compute_emission => compute_emission_RO2
        procedure :: compute_reaction_rate => reaction_rate_RO2
    end type

    ! type, extends(chemical_element) :: MVK_type
    ! PURPOSE: The MVK type derived from chemical_element
    type, extends(chemical_element) :: MVK_type
    contains
        procedure :: compute_emission => compute_emission_MVK
        procedure :: compute_reaction_rate => reaction_rate_MVK
    end type

    ! type, extends(chemical_element) :: H2O2_type
    ! PURPOSE: The H2O2 type derived from chemical_element
    type, extends(chemical_element) :: H2O2_type
    contains
        procedure :: compute_emission => compute_emission_H2O2
        procedure :: compute_reaction_rate => reaction_rate_H2O2
    end type

    ! type, extends(chemical_element) :: HNO3_type
    ! PURPOSE: The HNO3 type derived from chemical_element
    type, extends(chemical_element) :: HNO3_type
    contains
        procedure :: compute_emission => compute_emission_HNO3
        procedure :: compute_reaction_rate => reaction_rate_HNO3
    end type

    ! type, extends(chemical_element) :: NO3_type
    ! PURPOSE: The NO3 type derived from chemical_element
    type, extends(chemical_element) :: NO3_type
    contains
        procedure :: compute_emission => compute_emission_NO3
        procedure :: compute_reaction_rate => reaction_rate_NO3
    end type

    ! type, extends(chemical_element) :: N2O5_type
    ! PURPOSE: The N2O5 type derived from chemical_element
    type, extends(chemical_element) :: N2O5_type
    contains
        procedure :: compute_emission => compute_emission_N2O5
        procedure :: compute_reaction_rate => reaction_rate_N2O5
    end type

    ! type, extends(chemical_element) :: SO2_type
    ! PURPOSE: The SO2 type derived from chemical_element
    type, extends(chemical_element) :: SO2_type
    contains
        procedure :: compute_emission => compute_emission_SO2
        procedure :: compute_reaction_rate => reaction_rate_SO2
    end type

    ! type, extends(chemical_element) :: H2SO4_type
    ! PURPOSE: The H2SO4 type derived from chemical_element
    type, extends(chemical_element) :: H2SO4_type
    contains
        procedure :: compute_emission => compute_emission_H2SO4
        procedure :: compute_reaction_rate => reaction_rate_H2SO4
    end type

    ! type, extends(chemical_element) :: H2SO4_P_type
    ! PURPOSE: The H2SO4_P type derived from chemical_element
    type, extends(chemical_element) :: H2SO4_P_type
    contains
        procedure :: compute_emission => compute_emission_H2SO4_P
        procedure :: compute_reaction_rate => reaction_rate_H2SO4_P
    end type

    ! type, extends(chemical_element) :: HNO3_P_type
    ! PURPOSE: The HNO3_P type derived from chemical_element
    type, extends(chemical_element) :: HNO3_P_type
    contains
        procedure :: compute_emission => compute_emission_HNO3_P
        procedure :: compute_reaction_rate => reaction_rate_HNO3_P
    end type


    ! type, extends(chemical_element) :: ELVOC_type
    ! PURPOSE: The ELVOC type derived from chemical_element
    type, extends(chemical_element) :: ELVOC_type
    contains
        procedure :: compute_emission => compute_emission_ELVOC
        procedure :: compute_reaction_rate => reaction_rate_ELVOC
    end type






    ! type prognostics_type
    ! PURPOSE: Contains all of the prognostics (ua, va, theta and chemical components)
    type prognostics_type
        real(kind = 8), dimension(nz) :: ua, ua_mid, va, va_mid, theta, theta_mid ! Prognostics (*_mid are for leapfrog)
        real(kind = 8), dimension(nz-2) :: dudt, dvdt, dThetaDt ! Tendencies [m/s^2] and [K/s]
        real(kind = 8), dimension(nz) :: T, M, N2, O2, pressure
        type(alpha_pinene_type) :: alpha_pinene                 ! Instantation of alpha_pinene
        type(isoprene_type) :: isoprene                         ! Instantation of isoprene
        type(OH_type) :: OH
        type(O3_type) :: O3
        type(O1D_type) :: O1D
        type(REST_type) :: REST
        type(NO2_type) :: NO2
        type(NO_type) :: NO
        type(CH2O_type) :: CH2O
        type(HO2_type) :: HO2
        type(CO_type) :: CO
        type(CO2_type) :: CO2
        type(CH4_type) :: CH4
        type(CH3O2_type) :: CH3O2
        type(RO2_type) :: RO2
        type(MVK_type) :: MVK
        type(H2O2_type) :: H2O2
        type(HNO3_type) :: HNO3
        type(NO3_type) :: NO3
        type(N2O5_type) :: N2O5
        type(SO2_type) :: SO2
        type(H2SO4_P_type) :: H2SO4_P
        type(H2SO4_type) :: H2SO4
        type(HNO3_P_type) :: HNO3_P
        type(ELVOC_type) :: ELVOC
        real(kind = 8), dimension(2, 1:nz) :: cond_sink
        real(kind = 8), dimension(nr_bins, 1:nz) :: size_distribution
        real(kind = 8), dimension(nz) :: PN, PM, PV
    contains
        procedure :: leapfrog_middle                            ! Leapfrog to middle (NOT FULLY IMPLEMENTED)
        procedure :: leapfrog_next                              ! Leapfrog to next timestep (NOT FULLY IMPLEMENTED)
        procedure :: euler_next                                 ! Take an Euler step
        procedure :: init_chemical_elements                     ! Initialize the chemical components
        procedure :: compute_diagnostics
    end type prognostics_type







    ! type, abstract :: chemical_reaction
    ! PURPOSE: base type for chemical reactions.
    type, abstract :: chemical_reaction
        real(kind = 8) :: rate_coefficient(2:nz-1)
        real(kind = 8) :: rate(2:nz-1)
    contains
        procedure(compute_chemistry_interface), deferred, private :: compute_rate_coefficient
        procedure(compute_chemistry_interface), deferred :: compute_rate
    end type

    type, extends(chemical_reaction) :: O3_diss_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_O3_diss
        procedure :: compute_rate => rate_O3_diss
    end type

    type, extends(chemical_reaction) :: O1D_H2O_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_O1D_H2O
        procedure :: compute_rate => rate_O1D_H2O
    end type

    type, extends(chemical_reaction) :: O1D_N2_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_O1D_N2
        procedure :: compute_rate => rate_O1D_N2
    end type

    type, extends(chemical_reaction) :: O1D_O2_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_O1D_O2
        procedure :: compute_rate => rate_O1D_O2
    end type

    type, extends(chemical_reaction) :: NO2_diss_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_NO2_diss
        procedure :: compute_rate => rate_NO2_diss
    end type

    type, extends(chemical_reaction) :: CH2O_diss_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_CH2O_diss
        procedure :: compute_rate => rate_CH2O_diss
    end type

    type, extends(chemical_reaction) :: OH_CO_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_CO
        procedure :: compute_rate => rate_OH_CO
    end type

    type, extends(chemical_reaction) :: OH_CH4_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_CH4
        procedure :: compute_rate => rate_OH_CH4
    end type

    type, extends(chemical_reaction) :: OH_isoprene_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_isoprene
        procedure :: compute_rate => rate_OH_isoprene
    end type

    type, extends(chemical_reaction) :: OH_MVK_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_MVK
        procedure :: compute_rate => rate_OH_MVK
    end type

    type, extends(chemical_reaction) :: HO2_NO_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_HO2_NO
        procedure :: compute_rate => rate_HO2_NO
    end type

    type, extends(chemical_reaction) :: CH3O2_NO_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_CH3O2_NO
        procedure :: compute_rate => rate_CH3O2_NO
    end type

    type, extends(chemical_reaction) :: RO2_NO_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_RO2_NO
        procedure :: compute_rate => rate_RO2_NO
    end type

    type, extends(chemical_reaction) :: OH_CH2O_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_CH2O
        procedure :: compute_rate => rate_OH_CH2O
    end type

    type, extends(chemical_reaction) :: HO2_diss_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_HO2_diss
        procedure :: compute_rate => rate_HO2_diss
    end type

    type, extends(chemical_reaction) :: CH3O2_HO2_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_CH3O2_HO2
        procedure :: compute_rate => rate_CH3O2_HO2
    end type

    type, extends(chemical_reaction) :: RO2_HO2_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_RO2_HO2
        procedure :: compute_rate => rate_RO2_HO2
    end type

    type, extends(chemical_reaction) :: OH_NO2_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_NO2
        procedure :: compute_rate => rate_OH_NO2
    end type

    type, extends(chemical_reaction) :: NO_O3_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_NO_O3
        procedure :: compute_rate => rate_NO_O3
    end type

    type, extends(chemical_reaction) :: OH_HO2_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_HO2
        procedure :: compute_rate => rate_OH_HO2
    end type

    type, extends(chemical_reaction) :: OH_H2O2_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_H2O2
        procedure :: compute_rate => rate_OH_H2O2
    end type

    type, extends(chemical_reaction) :: NO_NO3_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_NO_NO3
        procedure :: compute_rate => rate_NO_NO3
    end type

    type, extends(chemical_reaction) :: NO2_O3_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_NO2_O3
        procedure :: compute_rate => rate_NO2_O3
    end type

    type, extends(chemical_reaction) :: NO2_NO3_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_NO2_NO3
        procedure :: compute_rate => rate_NO2_NO3
    end type

    type, extends(chemical_reaction) :: N2O5_diss_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_N2O5_diss
        procedure :: compute_rate => rate_N2O5_diss
    end type

    type, extends(chemical_reaction) :: N2O5_H2O_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_N2O5_H2O
        procedure :: compute_rate => rate_N2O5_H2O
    end type

    type, extends(chemical_reaction) :: N2O5_2H2O_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_N2O5_2H2O
        procedure :: compute_rate => rate_N2O5_2H2O
    end type

    type, extends(chemical_reaction) :: HO2_O3_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_HO2_O3
        procedure :: compute_rate => rate_HO2_O3
    end type

    type, extends(chemical_reaction) :: OH_SO2_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_OH_SO2
        procedure :: compute_rate => rate_OH_SO2
    end type

    type, extends(chemical_reaction) :: H2SO4_diss_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_H2SO4_diss
        procedure :: compute_rate => rate_H2SO4_diss
    end type

    type, extends(chemical_reaction) :: HNO3_diss_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_HNO3_diss
        procedure :: compute_rate => rate_HNO3_diss
    end type

    type, extends(chemical_reaction) :: alpha_pinene_OH_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_alpha_pinene_OH
        procedure :: compute_rate => rate_alpha_pinene_OH
    end type

    type, extends(chemical_reaction) :: alpha_pinene_O3_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_alpha_pinene_O3
        procedure :: compute_rate => rate_alpha_pinene_O3
    end type

    type, extends(chemical_reaction) :: isoprene_O3_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_isoprene_O3
        procedure :: compute_rate => rate_isoprene_O3
    end type

    type, extends(chemical_reaction) :: ELVOC_form_type
    contains
        procedure :: compute_rate_coefficient => rate_coefficient_ELVOC_form
        procedure :: compute_rate => rate_ELVOC_form
    end type

    type(O3_diss_type) :: O3_diss
    type(O1D_H2O_type) :: O1D_H2O
    type(O1D_N2_type) :: O1D_N2
    type(O1D_O2_type) :: O1D_O2
    type(NO2_diss_type) :: NO2_diss
    type(CH2O_diss_type) :: CH2O_DISS
    type(OH_CO_type) :: OH_CO
    type(OH_CH4_type) :: OH_CH4
    type(OH_isoprene_type) :: OH_isoprene
    type(OH_MVK_type) :: OH_MVK
    type(HO2_NO_type) :: HO2_NO
    type(CH3O2_NO_type) :: CH3O2_NO
    type(RO2_NO_type) :: RO2_NO
    type(OH_CH2O_type) :: OH_CH2O
    type(HO2_diss_type) :: HO2_diss
    type(CH3O2_HO2_type) :: CH3O2_HO2
    type(RO2_HO2_type) :: RO2_HO2
    type(OH_NO2_type) :: OH_NO2
    type(NO_O3_type) :: NO_O3
    type(OH_HO2_type) :: OH_HO2
    type(OH_H2O2_type) :: OH_H2O2
    type(NO_NO3_type) :: NO_NO3
    type(NO2_O3_type) :: NO2_O3
    type(NO2_NO3_type) :: NO2_NO3
    type(N2O5_diss_type) :: N2O5_diss
    type(N2O5_H2O_type) :: N2O5_H2O
    type(N2O5_2H2O_type) :: N2O5_2H2O
    type(HO2_O3_type) :: HO2_O3
    type(OH_SO2_type) :: OH_SO2
    type(H2SO4_diss_type) :: H2SO4_diss
    type(HNO3_diss_type) :: HNO3_diss
    type(alpha_pinene_OH_type) :: alpha_pinene_OH
    type(alpha_pinene_O3_type) :: alpha_pinene_O3
    type(isoprene_O3_type) :: isoprene_O3
    type(ELVOC_form_type) :: ELVOC_form


    ! Interface for the emission routie. Required for the deferred routine compute_emission in chemical_element.
    interface
        subroutine compute_emission_interface(this, progn)
            import
            class(chemical_element), intent(inout) :: this
            type(prognostics_type), intent(in) :: progn
        end subroutine
    end interface

    interface
        subroutine compute_chemistry_interface(this, progn, level)
            import
            class(chemical_reaction), intent(inout) :: this
            type(prognostics_type), intent(inout) :: progn
            integer, intent(in) :: level ! Height level [1, nz]. Only range [2, nz-1] should be used.
        end subroutine
    end interface

    interface
        subroutine compute_reaction_rate_interface(this, level)
            import
            class(chemical_element), intent(inout) :: this
            integer, intent(in) :: level ! Height level [1, nz]. Only range [2, nz-1] should be used.
        end subroutine
    end interface

contains

    ! subroutine leapfrog_middle(this)
    ! PURPOSE: Leapfrog step to middle point. (NOT FULY IMPLEMENTED)
    subroutine leapfrog_middle(this)
        implicit none
        class(prognostics_type), intent(inout) :: this

        this%ua_mid(updInd) = this%ua_mid(updInd) + dt * this%dudt
        this%va_mid(updInd) = this%va_mid(updInd) + dt * this%dvdt
        this%theta_mid(updInd) = this%theta_mid(updInd) + dt * this%dThetaDt
    end subroutine

    ! subroutine leapfrog_next(this)
    ! PURPOSE: Leapfrog to next point.
    subroutine leapfrog_next(this)
        implicit none
        class(prognostics_type), intent(inout) :: this
        call this%euler_next() ! Same form, however the tendencies refer to the midpoints
    end subroutine

    ! subroutine prognostics_type::euler_next(this)
    ! PURPOSE: Update prognostics to next time using previously computed tendencies.
    subroutine euler_next(this)
        implicit none
        class(prognostics_type), intent(inout) :: this

        ! Update u, v, theta
        this%ua(updInd) = this%ua(updInd) + dt * this%dudt
        this%va(updInd) = this%va(updInd) + dt * this%dvdt
        this%theta(updInd) = this%theta(updInd) + dt * this%dThetaDt

        ! Update chemical components
        this%O3%concentration(updInd) = this%O3%concentration(updInd) +  dt * (this%O3%dynamical_tendency)
        this%O1D%concentration(updInd) = this%O1D%concentration(updInd) +  dt * (this%O1D%dynamical_tendency)
        this%OH%concentration(updInd) = this%OH%concentration(updInd) +  dt * (this%OH%dynamical_tendency)
        this%REST%concentration(updInd) = this%REST%concentration(updInd) +  dt * (this%REST%dynamical_tendency)
        this%NO2%concentration(updInd) = this%NO2%concentration(updInd) +  dt * (this%NO2%dynamical_tendency)
        this%NO%concentration(updInd) = this%NO%concentration(updInd) +  dt * (this%NO%dynamical_tendency)
        this%CH2O%concentration(updInd) = this%CH2O%concentration(updInd) +  dt * (this%CH2O%dynamical_tendency)
        this%HO2%concentration(updInd) = this%HO2%concentration(updInd) +  dt * (this%HO2%dynamical_tendency)
        this%CO%concentration(updInd) = this%CO%concentration(updInd) +  dt * (this%CO%dynamical_tendency)
        this%CO2%concentration(updInd) = this%CO2%concentration(updInd) +  dt * (this%CO2%dynamical_tendency)
        this%CH4%concentration(updInd) = this%CH4%concentration(updInd) +  dt * (this%CH4%dynamical_tendency)
        this%CH3O2%concentration(updInd) = this%CH3O2%concentration(updInd) +  dt * (this%CH3O2%dynamical_tendency)
        this%isoprene%concentration(updInd) = this%isoprene%concentration(updInd) +  dt * (this%isoprene%dynamical_tendency)
        this%RO2%concentration(updInd) = this%RO2%concentration(updInd) +  dt * (this%RO2%dynamical_tendency)
        this%MVK%concentration(updInd) = this%MVK%concentration(updInd) +  dt * (this%MVK%dynamical_tendency)
        this%H2O2%concentration(updInd) = this%H2O2%concentration(updInd) +  dt * (this%H2O2%dynamical_tendency)
        this%HNO3%concentration(updInd) = this%HNO3%concentration(updInd) +  dt * (this%HNO3%dynamical_tendency)
        this%NO3%concentration(updInd) = this%NO3%concentration(updInd) +  dt * (this%NO3%dynamical_tendency)
        this%N2O5%concentration(updInd) = this%N2O5%concentration(updInd) +  dt * (this%N2O5%dynamical_tendency)
        this%SO2%concentration(updInd) = this%SO2%concentration(updInd) +  dt * (this%SO2%dynamical_tendency)
        this%H2SO4%concentration(updInd) = this%H2SO4%concentration(updInd) +  dt * (this%H2SO4%dynamical_tendency)
        this%H2SO4_P%concentration(updInd) = this%H2SO4_P%concentration(updInd) +  dt * (this%H2SO4_P%dynamical_tendency)
        this%alpha_pinene%concentration(updInd) = this%alpha_pinene%concentration(updInd)+dt*this%alpha_pinene%dynamical_tendency
        this%HNO3_P%concentration(updInd) = this%HNO3_P%concentration(updInd) +  dt * (this%HNO3_P%dynamical_tendency)
        this%ELVOC%concentration(updInd) = this%ELVOC%concentration(updInd) +  dt * (this%ELVOC%dynamical_tendency)


    end subroutine

    ! subroutine prognostics_type::init_chemical_element(this)
    ! PURPOSE: gfortran 5 still lacks type parameters, so initialize the member constants (molar_mass, Henry_const, reactivity) here.
    subroutine init_chemical_elements(this)
        implicit none
        class(prognostics_type), intent(inout) :: this
        integer :: i


        call this%alpha_pinene%init(M_alpha_pinene, H_alpha_pinene, f_alpha_pinene,'alpha_pinene')
        call this%isoprene%init(M_isoprene, H_isoprene, f_isoprene,'isoprene')
        call this%OH%init(17.01_dp, 0.0_dp, 0.0_dp, 'OH')
        call this%O3%init(0.0_dp, 0.0_dp, 0.0_dp, 'O3')
        call this%O1D%init(0.0_dp, 0.0_dp, 0.0_dp, 'O1D')
        call this%REST%init(0.0_dp, 0.0_dp, 0.0_dp, 'REST')
        call this%NO2%init(0.0_dp, 0.0_dp, 0.0_dp, 'NO2')
        call this%NO%init(0.0_dp, 0.0_dp, 0.0_dp, 'NO')
        call this%CH2O%init(0.0_dp, 0.0_dp, 0.0_dp, 'CH2O')
        call this%HO2%init(0.0_dp, 0.0_dp, 0.0_dp, 'HO2')
        call this%CO%init(0.0_dp, 0.0_dp, 0.0_dp, 'CO')
        call this%CO2%init(0.0_dp, 0.0_dp, 0.0_dp, 'CO2')
        call this%CH4%init(0.0_dp, 0.0_dp, 0.0_dp, 'CH4')
        call this%CH3O2%init(0.0_dp, 0.0_dp, 0.0_dp, 'CH3O2')
        call this%RO2%init(0.0_dp, 0.0_dp, 0.0_dp, 'RO2')
        call this%MVK%init(0.0_dp, 0.0_dp, 0.0_dp, 'MVK')
        call this%H2O2%init(0.0_dp, 0.0_dp, 0.0_dp, 'H2O2')
        call this%HNO3%init(0.0_dp, 0.0_dp, 0.0_dp, 'HNO3')
        call this%NO3%init(0.0_dp, 0.0_dp, 0.0_dp, 'NO3')
        call this%N2O5%init(0.0_dp, 0.0_dp, 0.0_dp, 'N2O5')
        call this%SO2%init(0.0_dp, 0.0_dp, 0.0_dp, 'SO2')
        call this%H2SO4%init(0.0_dp, 0.0_dp, 0.0_dp, 'H2SO4')
        call this%H2SO4_P%init(0.0_dp, 0.0_dp, 0.0_dp, 'H2SO4_P')
        call this%HNO3_P%init(0.0_dp, 0.0_dp, 0.0_dp, 'HNO3_P')
        call this%ELVOC%init(0.0_dp, 0.0_dp, 0.0_dp, 'ELVOC')


        this%O3%concentration = O3_conc
        this%NO2%concentration = NO2_conc
        this%NO%concentration = NO_conc
        this%CO%concentration = CO_conc
        this%CH4%concentration = CH4_conc
        this%SO2%concentration = SO2_conc
        if (box) then
            this%isoprene%concentration = 4.8D10!2.2E-9 * N_air
            this%alpha_pinene%concentration = 4.8D10!2.2E-9 * N_air
        end if

        this%cond_sink = 0.0
        this%PN = 0
        this%PM = 0
        this%PV = 0
        do i = 1, nz
            CALL Aerosol_init(diameter, particle_mass, particle_volume, this%size_distribution(:,i), &
            particle_density, nucleation_coef, molecular_mass, molar_mass, &
            molecular_volume, molecular_dia, mass_accomm)
        end do

    end subroutine

    subroutine compute_diagnostics(this, P0)
        implicit none
        class(prognostics_type), intent(inout) :: this
        real(kind = 8), intent(in) :: P0 ! Surface pressure [Pa]
        integer :: i

        this%T = this%theta - (g/Cp) * z

        this%pressure(1) = P0
        do i = 2, nz
            this%pressure(i) = this%pressure(i-1) * exp(-(z(i)-z(i-1))*Mair*g / (R * (this%T(i-1)+this%T(i))/2))
        end do

        this%M = this%pressure / (R*this%T) * NA * 1D-6
        this%O2 = 0.21*this%M
        this%N2 = 0.78*this%M

    end subroutine

    ! suborutine chemical_element::init(this, molar_mass, Henry_const, reactivity)
    ! PURPOSE: Write initial values to the member variables
    ! INPUTS:
    !    (real*8) molar_mass [g/mol], Henry_const [M/atm], reactivity [unitless]
    subroutine init(this, molar_mass, Henry_const, reactivity, name)
        implicit none
        class(chemical_element), intent(inout) :: this
        real(kind = 8), intent(in) :: molar_mass, Henry_const, reactivity
        character(len = *), intent(in) :: name

        ! Initialize constants
        this%molar_mass = molar_mass
        this%Henry_const = Henry_const
        this%reactivity = reactivity
        this%name = name
        this%file_created = .false.

        ! Initialize columns to zero
        this%concentration = 0.0
        this%dynamical_tendency = 0.0
        this%parameterized_tendency = 0.0
        this%reaction_rate = 0.0

        ! Compute resistances, which do not depend on the atmospheric or radiative consitions:
        this%r_m = 1.0 / (3.3E-4 * this%Henry_const + 100.0 * this%reactivity)
        this%r_lu = 2000.0 / (1E-5*this%Henry_const + this%reactivity)
        this%r_cl = 1.0 / (1E-5*this%Henry_const/2000.0 + this%reactivity/1000.0)
        this%r_gs = 1.0 / (1E-5*this%Henry_const/500 + this%reactivity/200)

    end subroutine

    subroutine output(this)
        implicit none
        class(chemical_element), intent(inout) :: this
        CHARACTER(255) :: outfmt

        WRITE(outfmt, '(a, i3, a)') '(', nz, 'es25.16)'

        if (.not. this%file_created) then
            OPEN(FILE = TRIM(ADJUSTL(outdir))//adjustl('/')//adjustl(trim(this%name))//adjustl('.dat')  , &
                STATUS = 'REPLACE', ACTION = 'WRITE', newunit=this%file_unit)
            this%file_created = .true.
        end if

        WRITE(this%file_unit, outfmt) this%concentration

    end subroutine

    subroutine close_file(this)
        implicit none
        class(chemical_element), intent(inout) :: this

        close(this%file_unit)

    end subroutine

    !
    ! Emission routines
    !

    subroutine compute_emission(this, progn)
        implicit none
        class(chemical_element), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_OH(this, progn)
        implicit none
        class(OH_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_O3(this, progn)
        implicit none
        class(O3_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_O1D(this, progn)
        implicit none
        class(O1D_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_REST(this, progn)
        implicit none
        class(REST_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_NO2(this, progn)
        implicit none
        class(NO2_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_NO(this, progn)
        implicit none
        class(NO_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_CH2O(this, progn)
        implicit none
        class(CH2O_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_HO2(this, progn)
        implicit none
        class(HO2_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_CO(this, progn)
        implicit none
        class(CO_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_CO2(this, progn)
        implicit none
        class(CO2_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_CH4(this, progn)
        implicit none
        class(CH4_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_CH3O2(this, progn)
        implicit none
        class(CH3O2_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_RO2(this, progn)
        implicit none
        class(RO2_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_MVK(this, progn)
        implicit none
        class(MVK_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_H2O2(this, progn)
        implicit none
        class(H2O2_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_HNO3(this, progn)
        implicit none
        class(HNO3_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_NO3(this, progn)
        implicit none
        class(NO3_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_N2O5(this, progn)
        implicit none
        class(N2O5_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_SO2(this, progn)
        implicit none
        class(SO2_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_H2SO4(this, progn)
        implicit none
        class(H2SO4_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_H2SO4_P(this, progn)
        implicit none
        class(H2SO4_P_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_HNO3_P(this, progn)
        implicit none
        class(HNO3_P_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine

    subroutine compute_emission_ELVOC(this, progn)
        implicit none
        class(ELVOC_type), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
    end subroutine


    ! subroutine alpha_pinene_type::compute_emission_alpha_pinene(this, progn)
    ! PURPOSE: Conpute the emission of alpha pinene at level 2 (this%emission) [# / cm^3 / s]
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_emission_alpha_pinene(this, progn)
        implicit none
        class(alpha_pinene_type), intent(inout) :: this
        ! Inputs
        type(prognostics_type), intent(in) :: progn
        ! Local constants and variables
        real(kind = 8), parameter :: B = 0.09 ! [1/K]
        real(kind = 8), parameter :: Ts = 303.15 ! [K]
        real(kind = 8) :: adjust_factor, T

        T = progn%T(2) ! Temperature at level 2
        adjust_factor = exp(B * (T - Ts) ) ! Adjust emission depending on temperature

        ! Compute emission from vegetation.
        this%emission = compute_vegetation_emission(adjust_factor, this%molar_mass)

    end subroutine compute_emission_alpha_pinene

    ! subroutine isoprene_type::compute_emission_isoprene(this, progn)
    ! PURPOSE: Conpute the emission of isoprene at level 2 (this%emission) [# / cm^3 / s]
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_emission_isoprene(this, progn)
        implicit none
        class(isoprene_type), intent(inout) :: this
        ! Inputs
        type(prognostics_type), intent(in) :: progn
        ! Constants for isoprene computation.
        real(kind = 8), parameter :: a = 0.0027, c_L1 = 1.006, c_T1 = 95E3, c_T2 = 230E3
        real(kind = 8), parameter :: Ts = 303.15, Tm = 314
        ! Local variables
        real(kind = 8) :: adjust_factor, T, C_L, C_T

        T = progn%T(2) ! Temperature at level 2

        ! Adjust emission depending on radiation and temperature.
        C_L = (a * c_L1 * PAR) / sqrt(1 + a**2 * PAR**2)
        C_T = exp(c_T1*(T-Ts)/(R*T*Ts)) / (1 + exp(c_T2*(T-Tm) / R*T*Ts) )
        adjust_factor = C_L * C_T

        ! Compute emission from vegetation.
        this%emission = compute_vegetation_emission(adjust_factor, this%molar_mass)


    end subroutine compute_emission_isoprene

    ! function compute_vegetation_emission(adjust_factor, molar_mass)
    ! PURPOSE: Compute emission from vegetation
    ! INPUT:
    !   (real*8) adjust_factor, molar_mass
    function compute_vegetation_emission(adjust_factor, molar_mass) result(F_vegetation_num)
        use parameters_mod
        implicit none
        ! Inputs:
        real(kind = 8), intent(in) :: adjust_factor, molar_mass
        ! Local variables
        real(kind = 8) :: F_vegetation_num, F_vegetation_mass

        ! Compute mass emission flux
        F_vegetation_mass = foliar_density * adjust_factor * emission_factor * emission_activity * 1E-9 / one_hour ! [g / cm^2 / s]

        ! Convert mass flux per area to number flux per volume
        F_vegetation_num = F_vegetation_mass * NA / molar_mass / 1000 ! [# / cm^3 / s]

    end function

    !
    ! Deposition routines
    !

    ! subroutine compute_aerodynamic_resistance(progn)
    ! PURPOSE: Compute atmospheric resistance [called from parameterizaions_mod]
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_aerodynamic_resistance(progn)
        implicit none
        type(prognostics_type), intent(in) :: progn
        real(kind = 8), parameter :: z_0m = 0.1, z = 10.0, d = 2.0*10/3, z_0h = 0.1*z_0m, k = 0.4
        real(kind = 8) :: wind_speed ! [m/s] at 10 meters (level 2)

        wind_speed = sqrt(progn%ua(2)**2 + progn%va(2)**2) ! Wind speed at level 2
        r_a = log((z-d)/z_0m) * log((z-d)/z_0h) / (k**2 * wind_speed) ! Atmospheric resistance

    end subroutine

    ! subroutine chemical_element::compute_quasi_laminar_resistance(this, progn)
    ! PURPOSE: Compute r_b
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_quasi_laminar_resistance(this, progn)
        implicit none
        class(chemical_element), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
        ! Local variables and constants
        real(kind = 8) :: wind_speed ! [m/s] at 10 meters (level 2)
        real(kind = 8) :: u_star, D_i
        real(kind = 8), parameter :: M_water = 18.01528 ! [g/mol]
        real(kind = 8), parameter :: kin_visc = 1.59E-5, k = 0.4, z_0m = 0.1, z = 10.0

        wind_speed = sqrt(progn%ua(2)**2 + progn%va(2)**2) ! Wind speed at level 2
        D_i = 2.4E-5 * sqrt(M_water / this%molar_mass)     ! Molecular diffusivity
        u_star = wind_speed * k / log(z / z_0m)            ! Fricional velocity

        this%r_b = 5*(kin_visc / D_i)**(2.0/3.0) / u_star  ! quasi-laminar resistance

    end subroutine

    ! subroutine compute_canopy_resistance(this, progn)
    ! PURPOSE: compute r_c
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_canopy_resistance(this, progn)
        implicit none
        class(chemical_element), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
        ! Local variables and constants
        real(kind = 8) :: G, T
        real(kind = 8), parameter :: M_water = 18.01528 ! [g/mol]
        real(kind = 8), parameter :: r_ac = 2000 ! [s/m]

        T = progn%T(2) - 273.15 ! [C]
        G = PAR/0.45 * 0.219 ! Solar radiation [W/m^2]

        ! Compute component resistances, which depend on the solar radiation
        if (0 <= T .and. T <= 40) then
            this%r_st = 130 * sqrt(this%molar_mass / M_water) * (1 + (200/(G+0.1))**2 * (400/(T*(40-T))))
        else
            this%r_st = r_max
        end if
        this%r_dc = 100 * (1 + 1000/(G+10))

        ! Compute r_c
        this%r_c = 1 / (1/(this%r_st+this%r_m) + 1/this%r_lu + 1/(this%r_dc+this%r_cl) + 1/(this%r_gs+r_ac))

    end subroutine

    ! function chemical_element::compute_deposition_velocity(this, progn)
    ! PURPOSE: Compute deposition velocity by combining r_a, r_b and r_c
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    ! OUSPUT:
    !   (real*8) deposition_velocity [m/s]
    function compute_deposition_velocity(this, progn) result(deposition_velocity)
        implicit none
        class(chemical_element), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
        real(kind = 8) :: deposition_velocity ! [m/s]

        call this%compute_canopy_resistance(progn) ! Update this%r_c
        call this%compute_quasi_laminar_resistance(progn) ! Update this%r_b

        ! Compute deposition velocity:
        deposition_velocity = 1 / (r_a + this%r_b + this%r_c)

    end function

    ! subroutine chemical_element::compute_deposition(this, progn)
    ! PURPOSE: Compute deposition based on the deposition velocity
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_deposition(this, progn)
        implicit none
        class(chemical_element), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn

        ! Deposition flux F_dep = -v_dep * concentration
        this%deposition = -1.0 * this%compute_deposition_velocity(progn)*100 * this%concentration(2) / 1000 ! [# / cm^3 / s]

    end subroutine

    ! subroutine chemical_element::compute_parameterized_tendency(this, progn)
    ! PURPOSE: Compute tendency arising from the combination of emission and deposition.
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_parameterized_tendency(this, progn, level)
        implicit none
        class(chemical_element), intent(inout) :: this
        type(prognostics_type), intent(in) :: progn
        integer, intent(in) :: level

        this%parameterized_tendency(level) = 0.0

        if (.not. box .and. level == 2) then
            call this%compute_emission(progn) ! First, compute emissions
            call this%compute_deposition(progn) ! Second, compute the deposition.

            ! Tendency (at 10 meters)
            this%parameterized_tendency(level) = this%emission - this%deposition
        end if

        ! Chemical reactions.
        call this%compute_reaction_rate(level)
        this%parameterized_tendency(level) = this%parameterized_tendency(level) + this%reaction_rate(level)

    end subroutine

    ! subroutine chemical_element::compute_dynamical_tendency(this, Kh)
    ! PURPOSE: Compute tendency due to atmospheric turbulence.
    ! INPUT:
    !   (real*8) Kh(nz-1) [Turbulence coefficients at mid levels]
    subroutine compute_dynamical_tendency(this, Kh)
        use derivatives_mod
        implicit none
        class(chemical_element), intent(inout) :: this
        real(kind = 8), intent(in) :: Kh(:)

        ! (dC/dt)_dyn = d(Kh * dC/dz) / dz
        this%dynamical_tendency = zDerivMidLevel(Kh * zDeriv(this%concentration))

    end subroutine

    ! CHEMICAL REACTIONS

    function concentration(progn, level, species)
        implicit none
        type(prognostics_type), intent(in) :: progn
        integer, intent(in) :: level
        class(chemical_element), intent(in) :: species
        real(kind = 8) :: concentration

        concentration = species%concentration(level)

    end function

    function concentration_multiply(progn, level, species1, species2)
        implicit none
        type(prognostics_type), intent(in) :: progn
        integer, intent(in) :: level
        class(chemical_element), intent(in) :: species1, species2
        real(kind = 8) :: concentration_multiply

        concentration_multiply = species1%concentration(level) * &
                                 species2%concentration(level)

    end function

    ! REACTION RATE COEFFICIENTS

    subroutine rate_coefficient_O3_diss(this, progn, level)
        implicit none
        class(O3_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 3.83D-5 * get_exp_coszen(time, daynumber, latitude)

    end subroutine

    subroutine rate_coefficient_O1D_H2O(this, progn, level)
        implicit none
        class(O1D_H2O_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.63D-10 * exp(60.0_dp / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_O1D_N2(this, progn, level)
        implicit none
        class(O1D_N2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 2.15D-11 * exp(110.0_dp / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_O1D_O2(this, progn, level)
        implicit none
        class(O1D_O2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 3.3D-11 * exp(55.0_dp / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_NO2_diss(this, progn, level)
        implicit none
        class(NO2_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.67D-2 * get_exp_coszen(time, daynumber, latitude)

    end subroutine

    subroutine rate_coefficient_CH2O_diss(this, progn, level)
        implicit none
        class(CH2O_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.47D-4 * get_exp_coszen(time, daynumber, latitude)

    end subroutine

    subroutine rate_coefficient_OH_CO(this, progn, level)
        implicit none
        class(OH_CO_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 2.4D-13

    end subroutine

    subroutine rate_coefficient_OH_CH4(this, progn, level)
        implicit none
        class(OH_CH4_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 2.45D-12 * exp(-1775 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_OH_isoprene(this, progn, level)
        implicit none
        class(OH_isoprene_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1D-10

    end subroutine

    subroutine rate_coefficient_OH_MVK(this, progn, level)
        implicit none
        class(OH_MVK_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 2.4D-11

    end subroutine

    subroutine rate_coefficient_HO2_NO(this, progn, level)
        implicit none
        class(HO2_NO_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 3.5D-12 * exp(250 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_CH3O2_NO(this, progn, level)
        implicit none
        class(CH3O2_NO_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 2.8D-12 * exp(300 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_RO2_NO(this, progn, level)
        implicit none
        class(RO2_NO_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1D-11

    end subroutine

    subroutine rate_coefficient_OH_CH2O(this, progn, level)
        implicit none
        class(OH_CH2O_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 5.5D-12 * exp(125 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_HO2_diss(this, progn, level)
        implicit none
        class(HO2_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = ((2.2D-13*EXP(600/progn%T(level)))+ &
                                        (1.9D-33*EXP(980/progn%T(level))*progn%M(level)))* &
                                        (1+(1+1.4D-21*EXP(2200/progn%T(level))*H2O)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_CH3O2_HO2(this, progn, level)
        implicit none
        class(CH3O2_HO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 4.1D-13 * exp(750 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_RO2_HO2(this, progn, level)
        implicit none
        class(RO2_HO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.5D-11

    end subroutine

    subroutine rate_coefficient_OH_NO2(this, progn, level)
        implicit none
        class(OH_NO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 3.5D-12 * exp(340 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_NO_O3(this, progn, level)
        implicit none
        class(NO_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 3D-12 * exp(-1500 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_OH_HO2(this, progn, level)
        implicit none
        class(OH_HO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 4.8D-11 * exp(250 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_OH_H2O2(this, progn, level)
        implicit none
        class(OH_H2O2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 2.9D-12 * exp(-160 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_NO_NO3(this, progn, level)
        implicit none
        class(NO_NO3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.8D-11 * exp(110 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_NO2_O3(this, progn, level)
        implicit none
        class(NO2_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.4D-13 * exp(-2470 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_NO2_NO3(this, progn, level)
        implicit none
        class(NO2_NO3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = (0.35*(3.6D-30*(progn%T(level)/300)**(-4.1)*progn%M(level))* &
                                        (1.9D-12*(progn%T(level)/300)**0.2)) / &
                                        ((3.6D-30*(progn%T(level)/300)**(-4.1)*progn%M(level))+ &
                                        (1.9D-12*(progn%T(level)/300)**0.2))  ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_N2O5_diss(this, progn, level)
        implicit none
        class(N2O5_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = (0.35*(1.3D-3*(progn%T(level)/300)**(-3.5)*EXP(-11000/progn%T(level))*progn%M(level))* &
                                        (9.7D14*(progn%T(level)/300)**0.1*EXP(- 11080/progn%T(level)))) / &
                                        ((1.3D-3*(progn%T(level)/300)**(-3.5)*EXP(-11000/progn%T(level))*progn%M(level))+ &
                                        (9.7D14*(progn%T(level)/300)**0.1*EXP(- 11080/progn%T(level)))) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_N2O5_H2O(this, progn, level)
        implicit none
        class(N2O5_H2O_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 2.5D-22

    end subroutine

    subroutine rate_coefficient_N2O5_2H2O(this, progn, level)
        implicit none
        class(N2O5_2H2O_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.8D-39

    end subroutine

    subroutine rate_coefficient_HO2_O3(this, progn, level)
        implicit none
        class(HO2_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 2.03D-16*(progn%T(level)/300)**4.57*EXP(693/progn%T(level)) ! TODO

    end subroutine

    subroutine rate_coefficient_OH_SO2(this, progn, level)
        implicit none
        class(OH_SO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.5D-12

    end subroutine

    subroutine rate_coefficient_H2SO4_diss(this, progn, level)
        implicit none
        class(H2SO4_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level


        if (box .or. .not. aerosol_coupling) then
            this%rate_coefficient(level) = 1D-3!progn%cond_sink(1,level)
        else
            this%rate_coefficient(level) = progn%cond_sink(1,level)
        end if
        this%rate_coefficient(level) = 1D-3 ! TESTI
        !print *, this%rate_coefficient(level)

    end subroutine

    subroutine rate_coefficient_HNO3_diss(this, progn, level)
        implicit none
        class(HNO3_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1D-3

    end subroutine

    subroutine rate_coefficient_alpha_pinene_OH(this, progn, level)
        implicit none
        class(alpha_pinene_OH_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.2D-11 * exp(440 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_alpha_pinene_O3(this, progn, level)
        implicit none
        class(alpha_pinene_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 6.3D-16 * exp(-580 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_isoprene_O3(this, progn, level)
        implicit none
        class(isoprene_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate_coefficient(level) = 1.03D-14 * exp(-1995 / progn%T(level)) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

    end subroutine

    subroutine rate_coefficient_ELVOC_form(this, progn, level)
        implicit none
        class(ELVOC_form_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        if (box .or. .not. aerosol_coupling) then
            this%rate_coefficient(level) = 1D-3!progn%cond_sink(2,level)
        else
            this%rate_coefficient(level) = progn%cond_sink(2,level)
        end if

    end subroutine

    ! REACTION RATES

    subroutine rate_O3_diss(this, progn, level)
        implicit none
        class(O3_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%O3)

    end subroutine

    subroutine rate_O1D_H2O(this, progn, level)
        implicit none
        class(O1D_H2O_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%O1D) * H2O

    end subroutine

    subroutine rate_O1D_N2(this, progn, level)
        implicit none
        class(O1D_N2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%O1D) * progn%N2(level)

    end subroutine

    subroutine rate_O1D_O2(this, progn, level)
        implicit none
        class(O1D_O2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%O1D) * progn%O2(level)

    end subroutine

    subroutine rate_NO2_diss(this, progn, level)
        implicit none
        class(NO2_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%NO2)

    end subroutine

    subroutine rate_CH2O_diss(this, progn, level)
        implicit none
        class(CH2O_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%CH2O)

    end subroutine

    subroutine rate_OH_CO(this, progn, level)
        implicit none
        class(OH_CO_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level,progn%OH, progn%CO)

    end subroutine

    subroutine rate_OH_CH4(this, progn, level)
        implicit none
        class(OH_CH4_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%OH, progn%CH4)

    end subroutine

    subroutine rate_OH_isoprene(this, progn, level)
        implicit none
        class(OH_isoprene_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%OH, progn%isoprene)

    end subroutine

    subroutine rate_OH_MVK(this, progn, level)
        implicit none
        class(OH_MVK_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%OH, progn%MVK)

    end subroutine

    subroutine rate_HO2_NO(this, progn, level)
        implicit none
        class(HO2_NO_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%HO2, progn%NO)

    end subroutine

    subroutine rate_CH3O2_NO(this, progn, level)
        implicit none
        class(CH3O2_NO_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%CH3O2, progn%NO)

    end subroutine

    subroutine rate_RO2_NO(this, progn, level)
        implicit none
        class(RO2_NO_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%RO2, progn%NO)

    end subroutine

    subroutine rate_OH_CH2O(this, progn, level)
        implicit none
        class(OH_CH2O_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%OH, progn%CH2O)

    end subroutine

    subroutine rate_HO2_diss(this, progn, level)
        implicit none
        class(HO2_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%HO2, progn%HO2)

    end subroutine

    subroutine rate_CH3O2_HO2(this, progn, level)
        implicit none
        class(CH3O2_HO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%CH3O2, progn%HO2)

    end subroutine

    subroutine rate_RO2_HO2(this, progn, level)
        implicit none
        class(RO2_HO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%RO2, progn%HO2)

    end subroutine

    subroutine rate_OH_NO2(this, progn, level)
        implicit none
        class(OH_NO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%OH, progn%NO2)

    end subroutine

    subroutine rate_NO_O3(this, progn, level)
        implicit none
        class(NO_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%NO, progn%O3)

    end subroutine

    subroutine rate_OH_HO2(this, progn, level)
        implicit none
        class(OH_HO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%OH, progn%HO2)

    end subroutine

    subroutine rate_OH_H2O2(this, progn, level)
        implicit none
        class(OH_H2O2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%OH, progn%H2O2)

    end subroutine

    subroutine rate_NO_NO3(this, progn, level)
        implicit none
        class(NO_NO3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%NO, progn%NO3)

    end subroutine

    subroutine rate_NO2_O3(this, progn, level)
        implicit none
        class(NO2_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%NO2, progn%O3)

    end subroutine

    subroutine rate_NO2_NO3(this, progn, level)
        implicit none
        class(NO2_NO3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%NO2, progn%NO3)

    end subroutine

    subroutine rate_N2O5_diss(this, progn, level)
        implicit none
        class(N2O5_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%N2O5)

    end subroutine

    subroutine rate_N2O5_H2O(this, progn, level)
        implicit none
        class(N2O5_H2O_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%N2O5) * H2O

    end subroutine

    subroutine rate_N2O5_2H2O(this, progn, level)
        implicit none
        class(N2O5_2H2O_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%N2O5) * H2O * H2O

    end subroutine

    subroutine rate_HO2_O3(this, progn, level)
        implicit none
        class(HO2_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%HO2, progn%O3)

    end subroutine

    subroutine rate_OH_SO2(this, progn, level)
        implicit none
        class(OH_SO2_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%OH, progn%SO2)

    end subroutine

    subroutine rate_H2SO4_diss(this, progn, level)
        implicit none
        class(H2SO4_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%H2SO4)

    end subroutine

    subroutine rate_HNO3_diss(this, progn, level)
        implicit none
        class(HNO3_diss_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration(progn, level, progn%HNO3)

    end subroutine

    subroutine rate_alpha_pinene_OH(this, progn, level)
        implicit none
        class(alpha_pinene_OH_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%alpha_pinene, progn%OH)

    end subroutine

    subroutine rate_alpha_pinene_O3(this, progn, level)
        implicit none
        class(alpha_pinene_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%alpha_pinene, progn%O3)

    end subroutine

    subroutine rate_isoprene_O3(this, progn, level)
        implicit none
        class(isoprene_O3_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = this%rate_coefficient(level) * concentration_multiply(progn, level, progn%isoprene, progn%O3)

    end subroutine

    subroutine rate_ELVOC_form(this, progn, level)
        implicit none
        class(ELVOC_form_type), intent(inout) :: this
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: level

        this%rate(level) = (0.05 * alpha_pinene_OH%rate(level) + 0.1  * alpha_pinene_O3%rate(level)) - &
                            this%rate_coefficient(level) * concentration(progn, level, progn%ELVOC)

    end subroutine

    ! CONCENTRATION TENDENCIES DUE TO CHEMICAL REACTIONS

    subroutine compute_reaction_rate(this, level)
        implicit none
        class(chemical_element), intent(inout) :: this
        integer, intent(in) :: level
    end subroutine

    subroutine reaction_rate_O3(this, level)
        implicit none
        class(O3_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 0
    end subroutine

    subroutine reaction_rate_O1D(this, level)
        implicit none
        class(O1D_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = O3_diss%rate(level) - O1D_H2O%rate(level) - O1D_N2%rate(level) - O1D_O2%rate(level)
        !print *, O3_diss%rate(level)

    end subroutine

    subroutine reaction_rate_OH(this, level)
        implicit none
        class(OH_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 2*O1D_H2O%rate(level) - OH_CO%rate(level) - OH_CH4%rate(level) &
                                    -OH_isoprene%rate(level) - OH_MVK%rate(level) + HO2_NO%rate(level) &
                                    -OH_CH2O%rate(level) - OH_NO2%rate(level) - OH_HO2%rate(level) &
                                    -OH_H2O2%rate(level) + HO2_O3%rate(level) - OH_SO2%rate(level) &
                                    -alpha_pinene_OH%rate(level)

    end subroutine

    subroutine reaction_rate_REST(this, level)
        implicit none
        class(REST_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 0

    end subroutine

    subroutine reaction_rate_NO(this, level)
        implicit none
        class(NO_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 0

    end subroutine

    subroutine reaction_rate_NO2(this, level)
        implicit none
        class(NO2_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 0

    end subroutine

    subroutine reaction_rate_CH2O(this, level)
        implicit none
        class(CH2O_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = -CH2O_diss%rate(level) + OH_MVK%rate(level) + CH3O2_NO%rate(level) &
                                    + RO2_NO%rate(level) - OH_CH2O%rate(level)

    end subroutine

    subroutine reaction_rate_HO2(this, level)
        implicit none
        class(HO2_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = CH2O_diss%rate(level) + OH_CO%rate(level) + OH_MVK%rate(level)  &
                                    - HO2_NO%rate(level) + CH3O2_NO%rate(level) + RO2_NO%rate(level) &
                                    + OH_CH2O%rate(level) - 2*HO2_diss%rate(level) - CH3O2_HO2%rate(level) &
                                    - RO2_HO2%rate(level) - OH_HO2%rate(level) + OH_H2O2%rate(level) &
                                    - HO2_O3%rate(level)

    end subroutine

    subroutine reaction_rate_CO(this, level)
        implicit none
        class(CO_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 0

    end subroutine

    subroutine reaction_rate_CO2(this, level)
        implicit none
        class(CO2_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 0

    end subroutine

    subroutine reaction_rate_CH4(this, level)
        implicit none
        class(CH4_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 0

    end subroutine

    subroutine reaction_rate_CH3O2(this, level)
        implicit none
        class(CH3O2_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = OH_CH4%rate(level) - CH3O2_NO%rate(level) - CH3O2_HO2%rate(level)

    end subroutine

    subroutine reaction_rate_isoprene(this, level)
        implicit none
        class(isoprene_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = -OH_isoprene%rate(level) - isoprene_O3%rate(level)

    end subroutine

    subroutine reaction_rate_RO2(this, level)
        implicit none
        class(RO2_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = OH_isoprene%rate(level) - RO2_NO%rate(level) - RO2_HO2%rate(level)

    end subroutine

    subroutine reaction_rate_MVK(this, level)
        implicit none
        class(MVK_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = -OH_MVK%rate(level) + RO2_NO%rate(level)

    end subroutine

    subroutine reaction_rate_H2O2(this, level)
        implicit none
        class(H2O2_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = HO2_diss%rate(level) - OH_H2O2%rate(level)

    end subroutine

    subroutine reaction_rate_HNO3(this, level)
        implicit none
        class(HNO3_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = OH_NO2%rate(level) + 2*N2O5_H2O%rate(level) + 2*N2O5_2H2O%rate(level) - HNO3_diss%rate(level)

    end subroutine

    subroutine reaction_rate_NO3(this, level)
        implicit none
        class(NO3_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = -NO_NO3%rate(level) + NO2_O3%rate(level) - NO2_NO3%rate(level) + N2O5_diss%rate(level)

    end subroutine

    subroutine reaction_rate_N2O5(this, level)
        implicit none
        class(N2O5_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = NO2_NO3%rate(level) - N2O5_diss%rate(level) - N2O5_H2O%rate(level) - N2O5_2H2O%rate(level)

    end subroutine

    subroutine reaction_rate_SO2(this, level)
        implicit none
        class(SO2_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = 0

    end subroutine

    subroutine reaction_rate_H2SO4(this, level)
        implicit none
        class(H2SO4_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = OH_SO2%rate(level) - H2SO4_diss%rate(level)

    end subroutine

    subroutine reaction_rate_H2SO4_P(this, level)
        implicit none
        class(H2SO4_P_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = H2SO4_diss%rate(level)

    end subroutine

    subroutine reaction_rate_HNO3_P(this, level)
        implicit none
        class(HNO3_P_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = HNO3_diss%rate(level)

    end subroutine

    subroutine reaction_rate_ELVOC(this, level)
        implicit none
        class(ELVOC_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = ELVOC_form%rate(level)

    end subroutine

    subroutine reaction_rate_alpha_pinene(this, level)
        implicit none
        class(alpha_pinene_type), intent(inout) :: this
        integer, intent(in) :: level

        this%reaction_rate(level) = -alpha_pinene_OH%rate(level) - alpha_pinene_O3%rate(level)
    end subroutine

end module prognostics_mod
