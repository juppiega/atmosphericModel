! module radiation_mod
! PURPOSE: Compute solar radiation.
module radiation_mod
    use parameters_mod
    use time_mod
    implicit none

    real(kind = 8), public :: PAR

contains

! subroutine compute_radiation()
! PURPOSE: Compute PAR. Called in the parameterizations_mod.
subroutine compute_radiation()
    implicit none

    PAR = 1000.0 * get_exp_coszen(time, daynumber, latitude)

end subroutine

! You can add this function in the module meteorology
REAL(dp) FUNCTION get_exp_coszen(time,daynumber,latitude)
  REAL(dp), INTENT(in) :: time,latitude
  INTEGER, INTENT(in) :: daynumber
  REAL(dp) :: hourangle,zenith,coszen
  hourangle = get_hourangle(time)
  zenith = solar_zenith_angle(hourangle,daynumber,latitude)
  coszen = COS(zenith)
  IF (coszen > 0) THEN  ! sun is above horizon
     get_exp_coszen = EXP(-0.575/coszen)
  ELSE
     get_exp_coszen = 0
  ENDIF
END FUNCTION get_exp_coszen


REAL(dp) FUNCTION get_hourangle(time)
  REAL(dp), INTENT(in) :: time
  REAL(dp), PARAMETER :: one_day = 24*one_hour
  get_hourangle = MODULO(time,one_day)/one_day * 2 * pi - pi
END FUNCTION get_hourangle


REAL(dp) FUNCTION solar_zenith_angle(hourangle,daynumber,latitude)
  ! http://en.wikipedia.org/wiki/Solar_elevation_angle
  ! http://en.wikipedia.org/wiki/Position_of_the_Sun
  INTEGER, INTENT(in) :: daynumber
  REAL(dp), INTENT(in) :: hourangle,latitude
  REAL(dp) :: declination,elevation
  REAL(dp), PARAMETER :: to_rad = pi/180.

  declination = -23.44 * to_rad * COS(2 * pi * (daynumber + 10)/365.)
  elevation = COS(hourangle)*COS(declination)*COS(latitude) &
       + SIN(declination)*SIN(latitude)
  solar_zenith_angle = pi/2. - elevation
  ! Notes:
  ! - Not tested near equador or on the southern hemisphere.
  ! - solar_zenith_angle can be larger than pi/2, it just means
  !   the sun is below horizon.
  ! - solar_zenith_angle assumes time is in local solar time, which
  !   is usually not exactly true
END FUNCTION solar_zenith_angle

end module radiation_mod
