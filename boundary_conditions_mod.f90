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

    progn%E_tot(nz-1) = 1D-2

    ! Read surface theta from file.
    call surface_values(surfaceTheta, time)
    progn%q(1) = 0.622 * 288513966.0*exp(-4302.645/(surfaceTheta-29.65)) / 1013
    progn%theta(1) = surfaceTheta
    progn%theta_mid(1) = surfaceTheta ! For leapfrog, not fully implemented.

    ! Update chemical component lower boundary such that dc/dz = 0.
    progn%alpha_pinene%concentration(1) = progn%alpha_pinene%concentration(2)
    progn%isoprene%concentration(1) = progn%isoprene%concentration(2)
    progn%O1D%concentration(1) = progn%O1D%concentration(2)
    progn%OH%concentration(1) = progn%OH%concentration(2)
    progn%REST%concentration(1) = progn%REST%concentration(2)
    progn%CH2O%concentration(1) = progn%CH2O%concentration(2)
    progn%HO2%concentration(1) = progn%HO2%concentration(2)
    progn%CO2%concentration(1) = progn%CO2%concentration(2)
    progn%CH3O2%concentration(1) = progn%CH3O2%concentration(2)
    progn%RO2%concentration(1) = progn%RO2%concentration(2)
    progn%MVK%concentration(1) = progn%MVK%concentration(2)
    progn%H2O2%concentration(1) = progn%H2O2%concentration(2)
    progn%HNO3%concentration(1) = progn%HNO3%concentration(2)
    progn%NO3%concentration(1) = progn%NO3%concentration(2)
    progn%N2O5%concentration(1) = progn%N2O5%concentration(2)
    progn%H2SO4%concentration(1) = progn%H2SO4%concentration(2)
    progn%H2SO4_P%concentration(1) = progn%H2SO4_P%concentration(2)
    progn%HNO3_P%concentration(1) = progn%HNO3_P%concentration(2)
    progn%ELVOC%concentration(1) = progn%ELVOC%concentration(2)

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
!-----------------------------------------------------------------------------------------

  temperature = sin(2*PI*time/86400 + PI)*1 + 19 + 273.15 ! Marine PBL SST daily cycle


END SUBROUTINE surface_values

end module boundary_conditions_mod
