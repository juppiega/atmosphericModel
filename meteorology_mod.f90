!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Meteorology module
! - Compute eddy diffusivity coefficients
! - Compute wind fields
! - Compute theta field
! - Mix chemicals and aerosols
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE meteorology_mod

!-----------------------------------------------------------------------------------------
! Load modules
!-----------------------------------------------------------------------------------------
USE parameters_mod
USE time_mod
USE grid_mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE  ! This applies in all subprograms inside the module

PRIVATE
PUBLIC ua, va, theta  ! basic meteorology variables
PUBLIC meteorology_init, surface_values  ! functions

!
! Some constants
!
REAL(dp), PARAMETER :: lambda = 300.0_dp  ! maximum mixing length, meters
REAL(dp), PARAMETER :: vonk = 0.4_dp      ! von Karman constant, dimensionless

!
! Meteorological variables
!
REAL(dp), DIMENSION(nz) :: ua, va, theta

!
! For convenient
!
INTEGER :: I, J  ! used for loop


CONTAINS


!-----------------------------------------------------------------------------------------
! Interface
!
SUBROUTINE meteorology_init()
!
! Declaration
!

!
! Description
!
! Set initial values for ua, va, theta.
!-----------------------------------------------------------------------------------------
  !
  ! Wind velocity
  !
  ua = 0.0_dp
  ua(nz) = 10.0_dp
  ua(2:nz-1) = ua(nz) * z(2:nz-1)/z(nz) ! Linear profile

  va = 0.0_dp

  !
  ! Potential temperature
  !
  theta = 273.15_dp + 25.0_dp
  theta(nz) = 273.15_dp + 30.0_dp
END SUBROUTINE meteorology_init


!-----------------------------------------------------------------------------------------
! Interface
!
SUBROUTINE surface_values(temperature, time)
!
! Declaration
!
  REAL(dp), INTENT(in)            :: time ! input, in seconds
  REAL(dp), INTENT(out)           :: temperature ! output, in Kelvin
  LOGICAL, SAVE                   :: first_time = .TRUE.
  REAL(dp), DIMENSION(8,50), SAVE :: surface_data
  REAL(dp), DIMENSION(50), SAVE   :: temperature_data
  REAL(dp), PARAMETER             :: seconds_in_day = 24*60*60
  REAL(dp), PARAMETER             :: seconds_in_30min = 30*60
  INTEGER                         :: index
  REAL(dp) :: time24h, time30min, time24plus15, temp1, temp2, x
!
! Description
!
! Get surface temperature at specific time.
!
! Note: can also get water concentrantion, in ppt, if modify this
! subroutine to also use column 8.
!
! Data is taken from:
! http://www.atm.helsinki.fi/~junninen/smartSearch/smartSearch.php
!-----------------------------------------------------------------------------------------
  !
  ! Only when called for the first time, read in data from file
  ! With this trick, we don't need to open the file in the main program
  !
  IF (first_time) THEN
     OPEN(30, file='input/hyytiala_2011_8_10_t_h2o.dat', status='old')
     READ(30, *) surface_data
     temperature_data(1:50) = surface_data(7,1:50) ! in Celcius
     first_time = .FALSE.
  END IF

  time24h = MODULO(time, seconds_in_day) ! time modulo 24 hours
  time24plus15 = time24h + 15*60 ! time since 23:45 previous day
  time30min = MODULO(time24plus15, seconds_in_30min)
  index = 1 + FLOOR(time24plus15/seconds_in_30min)

  temp1 = temperature_data(index)
  temp2 = temperature_data(index + 1)
  x = time30min/seconds_in_30min

  !
  ! linear interpolation between previous and next temperature data value
  !
  temperature = temp1 + x*(temp2 - temp1) + 273.15_dp ! now in Kelvin

END SUBROUTINE surface_values


END MODULE meteorology_mod
