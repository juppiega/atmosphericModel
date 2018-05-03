!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Parameters module
!
! - All the parameters are defined in this module.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE parameters_mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE

PUBLIC

!
! Declare numbers in 64bit floating point
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format
!
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,300)

REAL(dp), PARAMETER :: PI   = 2*ASIN(1.0_dp)   ! [-], the constant pi
REAL(dp), PARAMETER :: g    = 9.81_dp          ! [m s-2], gravitational acceleration
REAL(dp), PARAMETER :: R    = 8.3144598_dp     ! [J mol-1 K-1], universal gas constant
REAL(dp), PARAMETER :: NA   = 6.022140857e23_dp  ! [molec mol-1], Avogadro's number 
REAL(dp), PARAMETER :: Mair = 28.96e-3_dp        ! [kg mol-1], mean molar mass of air
REAL(dp), PARAMETER :: N_air = 2.4E19_dp        ! [1/cm^3], Number density of air molecules at surface pressure.
REAL(dp), PARAMETER :: kb   = 1.38064852e-23_dp  ! [m2 kg s-2 K-1], Boltzmann constant
REAL(dp), PARAMETER :: Cp = 1012.0_dp          ! [J kg-1 K-1], air specific heat at constant pressure, in reality, it has slight temperature dependency
REAL(dp), PARAMETER :: Omega = 2*PI/(24.0_dp*60.0_dp*60.0_dp)  ! [rad s-1], Earth angular speed
real(kind = 8), parameter :: surf_pressure = 101325 ![Pa]
real(kind = 8), parameter :: z0 = 0.0002 ! Roughness length [m]

real(kind = 8) :: ug = 10, vg = 0, thetaUpperBoundary = 303.15 ! Boundary conditions

integer, parameter :: leapfrog_atMidpoint = 1, leapfrog_atPreviousFullTime = 2, euler_step = 3
integer, parameter :: euler = 1, leapfrog = 2, RK4 = 3

logical, parameter :: output_chemistry = .true.
logical, parameter :: model_chemistry = .true.
logical, parameter :: box = .false.
logical, parameter :: model_aerosols = .true.
CHARACTER(255), PARAMETER :: outdir = 'output'
integer :: scheme = 2 ! 1 = K_theory, 2 = TEMF

INTEGER, PARAMETER ::  nr_bins = 100           ! Number of particle size bins
INTEGER, PARAMETER ::  nr_cond = 2             ! Number of condensable vapours

!
! Latitude and longitude of Hyytiälä:
!
REAL(dp), PARAMETER :: latitude_deg = 30D0!61.8455_dp  ! [degN]
REAL(dp), PARAMETER :: longitude_deg = 24.2833_dp  ! [degE]
REAL(dp), PARAMETER :: latitude = latitude_deg * PI/180.0_dp  ! [rad]
REAL(dp), PARAMETER :: longitude = longitude_deg * PI/180.0_dp  ! [rad]

REAL(dp), PARAMETER :: fcor = 2*Omega*SIN(latitude)  ! Coriolis parameter

END MODULE parameters_mod
