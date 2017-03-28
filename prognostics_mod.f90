! module prognostics_mod
! PURPOSE: Contain values and governing equations for all of the variables advanced forwards in time.
module prognostics_mod

    use grid_mod
    use time_mod
    use radiation_mod
    use chemistry_mod

    implicit none

    private
    public prognostics_type
    public chemical_element, compute_aerodynamic_resistance

    ! type, abstract :: chemical_element
    ! PURPOSE: Base type, from which all chemical elements (e.g. isoprene) are derived.
    !          Most of the variables and governing equations are the same, but the
    !          emission routine has to be written separately for every element.
    ! TO ADD NEW CHEMICAL ELEMENT:
    ! 1. Derive a new type from chemical_element and write reaction equations
    ! 2. Add it to prognostics_type
    ! 3. Add it to prognostics_type::init_chemical_elements
    ! 4. Add it to prognostics_type::euler_next
    ! 5. Add it to set_boundary_conditions
    ! 6. Add it to parameterizations_mod
    ! 7. Add it to write_files
    type, abstract :: chemical_element
        real(kind = 8) :: concentration(nz) = 0.0 ! Concentration in the atmospheric column [# / cm^3]
        real(kind = 8) :: parameterized_tendency(nz-2) = 0.0    ! Tendency from emission, deposition, chemical reactions [# / cm^3 / s]
        real(kind = 8) :: dynamical_tendency(nz-2) = 0.0    ! Tendency from atmospheric turbulence (computed in dynamics_mod) [# / cm^3 / s]
        real(kind = 8) :: reaction_rate(nz-2) = 0.0         ! Reaction rate.
        real(kind = 8) :: molar_mass                  ! Molar mass [g / mol]
        real(kind = 8) :: emission = 0.0    ! Emission at 10 meters [# / cm^3 / s]
        real(kind = 8) :: deposition = 0.0  ! Deposition at the ground [# / cm^3 / s]
        real(kind = 8) :: Henry_const       ! Henry constant for resistance calculations [M / atm] (from chemistry_mod)
        real(kind = 8) :: reactivity        ! Reactivity (f0) for resistance calculations [unitless] (from chemistry_mod)
        real(kind = 8) :: r_st, r_m, r_lu, r_dc, r_cl, r_gs, r_b, r_c ! Molecule specific resistances [s/m]. r_a in chemistry_mod
    contains
        procedure(compute_emission_interface), deferred :: compute_emission ! Emission computation routine (written in the derived class)
        procedure(compute_reaction_rate_interface), deferred :: compute_reaction_rate ! Reaction rate computation routine (written in the derived class)
        ! Resistance and deposition computations
        procedure :: compute_quasi_laminar_resistance
        procedure :: compute_canopy_resistance
        procedure :: compute_deposition_velocity
        procedure :: compute_deposition

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
        procedure :: compute_reaction_rate => no_reaction
    end type

    ! type, extends(chemical_element) :: OH_type
    ! PURPOSE: The OH type derived from chemical_element
    type, extends(chemical_element) :: OH_type
    contains
        procedure :: compute_emission => no_emission
        procedure :: compute_reaction_rate => no_reaction
    end type






    ! type prognostics_type
    ! PURPOSE: Contains all of the prognostics (ua, va, theta and chemical components)
    type prognostics_type
        real(kind = 8), dimension(nz) :: ua, ua_mid, va, va_mid, theta, theta_mid ! Prognostics (*_mid are for leapfrog)
        real(kind = 8), dimension(nz-2) :: dudt, dvdt, dThetaDt ! Tendencies [m/s^2] and [K/s]
        type(alpha_pinene_type) :: alpha_pinene                 ! Instantation of alpha_pinene
        type(isoprene_type) :: isoprene                         ! Instantation of isoprene
        type(OH_type) :: OH                         ! Instantation of OH
    contains
        procedure :: leapfrog_middle                            ! Leapfrog to middle (NOT FULLY IMPLEMENTED)
        procedure :: leapfrog_next                              ! Leapfrog to next timestep (NOT FULLY IMPLEMENTED)
        procedure :: euler_next                                 ! Take an Euler step
        procedure :: init_chemical_elements                     ! Initialize the chemical components
    end type prognostics_type







    ! type, abstract :: chemical_reaction
    ! PURPOSE: base type for chemical reactions.
    type, abstract :: chemical_reaction
        real(kind = 8) :: rate_coefficient(nz-2)
        real(kind = 8) :: rate(nz-2)
    contains
        procedure(compute_chemistry_interface), deferred, private :: compute_rate_coefficient
        procedure(compute_chemistry_interface), deferred :: compute_rate
    end type

    type, extends(chemical_reaction) :: alpha_pinene_OH_type
    contains
        procedure :: compute_rate_coefficient => compute_rate_coefficient_alpha_pinene_OH
        procedure :: compute_rate => compute_rate_alpha_pinene_OH
    end type

    type(alpha_pinene_OH_type) :: alpha_pinene_OH






    ! Interface for the emission routie. Required for the deferred routine compute_emission in chemical_element.
    interface
        subroutine compute_emission_interface(this, progn)
            import
            class(chemical_element), intent(inout) :: this
            type(prognostics_type), intent(in) :: progn
        end subroutine
    end interface

    interface
        subroutine compute_chemistry_interface(this, progn)
            import
            class(chemical_reaction), intent(inout) :: this
            type(prognostics_type), intent(in) :: progn
        end subroutine
    end interface

    interface
        subroutine compute_reaction_rate_interface(this)
            import
            class(chemical_element), intent(inout) :: this
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
    this%alpha_pinene%concentration(updInd) = this%alpha_pinene%concentration(updInd) + &
    dt * (this%alpha_pinene%dynamical_tendency + this%alpha_pinene%parameterized_tendency)

    this%isoprene%concentration(updInd) = this%isoprene%concentration(updInd) + &
    dt * (this%isoprene%dynamical_tendency + this%isoprene%parameterized_tendency)

    this%OH%concentration(updInd) = this%OH%concentration(updInd) + &
    dt * (this%OH%dynamical_tendency + this%OH%parameterized_tendency)

end subroutine

! subroutine prognostics_type::init_chemical_element(this)
! PURPOSE: gfortran 5 still lacks type parameters, so initialize the member constants (molar_mass, Henry_const, reactivity) here.
subroutine init_chemical_elements(this)
    implicit none
    class(prognostics_type), intent(inout) :: this

    call this%alpha_pinene%init(M_alpha_pinene, H_alpha_pinene, f_alpha_pinene)
    call this%isoprene%init(M_isoprene, H_isoprene, f_isoprene)
    call this%OH%init(17.01_dp, 0.0_dp, 0.0_dp)

    this%OH%concentration = 1E6

end subroutine

! suborutine chemical_element::init(this, molar_mass, Henry_const, reactivity)
! PURPOSE: Write initial values to the member variables
! INPUTS:
!    (real*8) molar_mass [g/mol], Henry_const [M/atm], reactivity [unitless]
subroutine init(this, molar_mass, Henry_const, reactivity)
    implicit none
    class(chemical_element), intent(inout) :: this
    real(kind = 8), intent(in) :: molar_mass, Henry_const, reactivity

    ! Initialize constants
    this%molar_mass = molar_mass
    this%Henry_const = Henry_const
    this%reactivity = reactivity

    ! Initialize columns to zero
    this%concentration = 0.0
    this%dynamical_tendency = 0.0
    this%parameterized_tendency = 0.0

    ! Compute resistances, which do not depend on the atmospheric or radiative consitions:
    this%r_m = 1.0 / (3.3E-4 * this%Henry_const + 100.0 * this%reactivity)
    this%r_lu = 2000.0 / (1E-5*this%Henry_const + this%reactivity)
    this%r_cl = 1.0 / (1E-5*this%Henry_const/2000.0 + this%reactivity/1000.0)
    this%r_gs = 1.0 / (1E-5*this%Henry_const/500 + this%reactivity/200)

end subroutine

!
! Emission routines
!

subroutine no_emission(this, progn)
    implicit none
    class(isoprene_type), intent(inout) :: this
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

    T = progn%theta(2) ! Temperature at level 2
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

    T = progn%theta(2) ! Temperature at level 2

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

    T = progn%theta(2) - 273.15 ! [C]
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
subroutine compute_parameterized_tendency(this, progn)
    implicit none
    class(chemical_element), intent(inout) :: this
    type(prognostics_type), intent(in) :: progn

    call this%compute_emission(progn) ! First, compute emissions
    call this%compute_deposition(progn) ! Second, compute the deposition

    ! Tendency (at 10 meters)
    this%parameterized_tendency(1) = this%emission - this%deposition

    ! Chemical reactions.
    call this%compute_reaction_rate()
    this%parameterized_tendency = this%parameterized_tendency + this%reaction_rate

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

subroutine compute_rate_coefficient_alpha_pinene_OH(this, progn)
    implicit none
    class(chemical_reaction), intent(inout) :: this
    type(prognostics_type), intent(in) :: progn

    this%rate_coefficient = 1.2E-11 * exp(440 / progn%theta) ! TODO: Yksikko? Muuta potentiaalista todelliseksi.

end subroutine

subroutine compute_rate_alpha_pinene_OH(this, progn)
    implicit none
    class(chemical_reaction), intent(inout) :: this
    type(prognostics_type), intent(in) :: progn

    call this%compute_rate_coefficient(progn)
    this%rate = this%rate_coefficient * progn%alpha_pinene%concentration(2:nz-1) * progn%OH%concentration(2:nz-1)

end subroutine

! CONCENTRATION TENDENCIES DUE TO CHEMICAL REACTIONS
subroutine no_reaction(this)
    class(chemical_element), intent(inout) :: this
end subroutine

subroutine reaction_rate_alpha_pinene(this)
    implicit none
    class(chemical_element), intent(inout) :: this

    this%reaction_rate = -alpha_pinene_OH%rate

end subroutine

end module prognostics_mod
