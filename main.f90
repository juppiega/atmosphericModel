!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Main program
!
! - Simulate emissions and chemical reactions of gases, aerosol processes as well as 
!   transport of gases and aerosol particles within the planetary boundary layer with a
!   column model.
! - Check Fortran conventions at http://www.fortran90.org/src/best-practices.html
! - Check code conventions at
!   http://www.cesm.ucar.edu/working_groups/Software/dev_guide/dev_guide/node7.html
!
! BUILDING:
! -Execute "make all" in the Debug folder.
!
! RUNNING:
! -Execute command
!    $ Debug/atmosphericModel
! in atmosphericModel/ folder (where the output/ folder resides).
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PROGRAM main

!-----------------------------------------------------------------------------------------
! Load modules
!-----------------------------------------------------------------------------------------
USE parameters_mod
USE time_mod
USE grid_mod

USE boundary_conditions_mod
use dynamics_mod
use parameterizations_mod
use prognostics_mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE
type(prognostics_type) :: progn
integer :: solution_method = euler
!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------
CALL time_init()         ! initialize time, time step, date

CALL dynamics_init()     ! Initialize dynamics

call parameterizations_init() ! Initialize parameterizations

CALL open_files()        ! open files for future use

call prognostics_init(progn) ! Write initial conditions to ua, va, theta, chemistry

CALL write_files(progn, time)   ! write initial values


!-----------------------------------------------------------------------------------------
! Start main loop
!-----------------------------------------------------------------------------------------
DO WHILE (time <= time_end)

  ! Set prognostics boundary conditions
  CALL set_boundary_conditions(progn, time)

  ! Compute values at next time step: u(n+1) = u(n) + dt * f(n), where f = du/dt
  call compute_dynamics(progn, euler_step) ! Compute dynamics tendencies (turbulent fluxes)
  call compute_parameterizations(progn)    ! Compute parameterizations (chemistry emissions, depositions)
  call progn%euler_next()                  ! Advance to the next timestep using forward Euler


  !---------------------------------------------------------------------------------------
  ! Ending loop actions
  !---------------------------------------------------------------------------------------

  !
  ! Write data to files and time infomation to screen
  !
  IF ( time >= time_start_output .and. MOD( NINT((time - time_start)*10.0), NINT(dt_output*10.0)) == 0 ) THEN
    WRITE(*, '(a8, f8.3, a6)') 'time = ', time/one_hour, '  hours'
    CALL write_files(progn, time)
  END IF

  !
  ! Advance to next time step
  !
  time = time + dt

  !
  ! Count loop number
  !
  counter = counter + 1

END DO

!-----------------------------------------------------------------------------------------
! Finalization
!-----------------------------------------------------------------------------------------
!
! Close all the opened files
!
CALL close_files()

!
! Count total time steps
!
WRITE(*,*) counter, 'time steps'


CONTAINS


!-----------------------------------------------------------------------------------------
! Interface
!
SUBROUTINE open_files()
!
! Declaration
!
  CHARACTER(255), PARAMETER :: outdir = 'output'
!
! Description
!
! Open files.
!-----------------------------------------------------------------------------------------

  OPEN(11, FILE = TRIM(ADJUSTL(outdir))//'/time.dat'  , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(12, FILE = TRIM(ADJUSTL(outdir))//'/h.dat'     , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(13, FILE = TRIM(ADJUSTL(outdir))//'/ua.dat'    , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(14, FILE = TRIM(ADJUSTL(outdir))//'/va.dat'    , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(15, FILE = TRIM(ADJUSTL(outdir))//'/theta.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(16, FILE = TRIM(ADJUSTL(outdir))//'/Km.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(17, FILE = TRIM(ADJUSTL(outdir))//'/Kh.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
  OPEN(18, FILE = TRIM(ADJUSTL(outdir))//'/Ri.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
  if (output_chemistry) then
    OPEN(19, FILE = TRIM(ADJUSTL(outdir))//'/alpha_pinene.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
    OPEN(20, FILE = TRIM(ADJUSTL(outdir))//'/isoprene.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
  end if
END SUBROUTINE open_files


!-----------------------------------------------------------------------------------------
! Interface
!
SUBROUTINE write_files(progn, time)
!
! Declaration
!
  type(prognostics_type), intent(in) :: progn
  REAL(DP), intent(in) :: TIME  ! current time
  CHARACTER(255) :: outfmt
!
! Description
!
! Save data to files.
!-----------------------------------------------------------------------------------------
  !
  ! Get output format for arrays with nz layers
  !
  WRITE(outfmt, '(a, i3, a)') '(', nz, 'es25.16)'

  !
  ! Only save h one time at the beginning
  !
  IF (time == time_start) THEN
    WRITE(12, *) z
  END IF

  WRITE(11, '(f8.4)'   ) time/(24*one_hour)  ! [day], time
  WRITE(13, outfmt) progn%ua                  ! [m s-1], u wind
  WRITE(14, outfmt) progn%va                  ! [m s-1], v wind
  WRITE(15, outfmt) progn%theta               ! [K], potential temperature
  if (output_chemistry) then
    write(19, outfmt) progn%alpha_pinene%concentration
    write(20, outfmt) progn%isoprene%concentration
    !write(19, *) progn%alpha_pinene%deposition
    !write(20, *) progn%isoprene%deposition
  end if

  call parameterizations_output(progn) ! TODO: empty as of yet
END SUBROUTINE write_files


!-----------------------------------------------------------------------------------------
! Interface
!
SUBROUTINE close_files()
!
! Declaration
!

!
! Description
!
! Close opened files
!-----------------------------------------------------------------------------------------

  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)
  CLOSE(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
END SUBROUTINE close_files

subroutine prognostics_init(progn)
    implicit none
    type(prognostics_type), intent(inout) :: progn

    ! Zonal wind: u(0) = 0, u(nz) = ug, linear profile
    progn%ua = 0.0
    progn%ua(nz) = ug
    progn%ua(2:nz-1) = progn%ua(nz) * z(2:nz-1)/z(nz) ! Linear profile
    progn%ua_mid = progn%ua ! For leapfrog

    ! Meridional wind: 0 everywhere
    progn%va = 0.0_dp
    progn%va_mid = progn%va

    ! Potential temperature: Almost constant, but theta(nz) = 30 C
    progn%theta = 273.15 + 25.0
    progn%theta(nz) = 273.15 + 30.0
    progn%theta_mid = progn%theta

    ! Compute middle tendencies for leapfrog
    call compute_dynamics(progn, leapfrog_atPreviousFullTime)
    call compute_parameterizations(progn)

    progn%ua_mid(updInd) = progn%ua(updInd) + 0.5 * dt * progn%dudt
    progn%va_mid(updInd) = progn%va(updInd) + 0.5 * dt * progn%dvdt
    progn%theta_mid(updInd) = progn%theta(updInd) + 0.5 * dt * progn%dThetaDt
    !write(*,*) progn%ua_mid - progn%ua

    ! Initialize chemistry
    call progn%init_chemical_elements()

end subroutine

END PROGRAM main
