! module boundary_conditions_mod
! PURPOSE: Update boundary contitions
module boundary_conditions_mod
    use prognostics_mod
    implicit none
    private
    public set_boundary_conditions
contains

! subroutine set_boundary_conditions(progn, time)
! PURPOSE: Called every time step. Updates upper and lower boundary conditions.
! INPUT:
!   (prognostics_type) : progn [contains the full atmospheric state]
!   (real*8)           : time  [simulation time]
subroutine set_boundary_conditions(progn, time)
    implicit none
    type(prognostics_type), intent(inout) :: progn
    real(kind = 8), intent(in) :: time
    ! Locals:
    real(kind = 8) :: surfaceTheta

    ! Read surface theta from file.
    call surface_values(surfaceTheta, time)
    progn%theta(1) = surfaceTheta
    progn%theta_mid(1) = surfaceTheta ! For leapfrog, not fully implemented.

    ! Update chemical component lower boundary such that dc/dz = 0.
    progn%alpha_pinene%concentration(1) = progn%alpha_pinene%concentration(2)
    progn%isoprene%concentration(1) = progn%isoprene%concentration(2)

end subroutine

!-----------------------------------------------------------------------------------------
! Interface
!
SUBROUTINE surface_values(temperature, time)
!
! Declaration
!
  use parameters_mod, only: dp
  implicit none
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

end module boundary_conditions_mod
