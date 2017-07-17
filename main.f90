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
    use aerosol_mod

    !-----------------------------------------------------------------------------------------
    ! Variable declaration
    !-----------------------------------------------------------------------------------------
    IMPLICIT NONE
    type(prognostics_type) :: progn
    integer :: solution_method = euler
    real(kind = 8) :: previous_output_time = -1E30
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

        if (.not. box) then
            ! Set prognostics boundary conditions
            CALL set_boundary_conditions(progn, time)

            ! Compute values at next time step: u(n+1) = u(n) + dt * f(n), where f = du/dt
            call compute_dynamics(progn, euler_step) ! Compute dynamics tendencies (turbulent fluxes)
            call progn%compute_diagnostics(surf_pressure)
            call compute_parameterizations(progn)    ! Compute parameterizations (chemistry emissions, depositions)
            call progn%euler_next()                  ! Advance to the next timestep using forward Euler
        else
            call test_chemistry(progn)
         end if


        !---------------------------------------------------------------------------------------
        ! Ending loop actions
        !---------------------------------------------------------------------------------------

        !
        ! Write data to files and time infomation to screen
        !
        IF ( time >= time_start_output .and. time >= previous_output_time + dt_output ) THEN
            WRITE(*, '(a8, f12.8, a10)') 'time = ', time/one_hour, '  hours'
            CALL write_files(progn, time)
            previous_output_time = time
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
    CALL close_files(progn)

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
        OPEN(19, FILE = TRIM(ADJUSTL(outdir))//'/PN.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
        OPEN(20, FILE = TRIM(ADJUSTL(outdir))//'/PM.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
        OPEN(21, FILE = TRIM(ADJUSTL(outdir))//'/PV.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
        OPEN(22, FILE = TRIM(ADJUSTL(outdir))//'/size_distribution.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
    END SUBROUTINE open_files


    !-----------------------------------------------------------------------------------------
    ! Interface
    !
    SUBROUTINE write_files(progn, time)
        !
        ! Declaration
        !
        type(prognostics_type), intent(inout) :: progn
        REAL(DP), intent(in) :: TIME  ! current time
        CHARACTER(255) :: outfmt, outfmt_nrbins
        !
        ! Description
        !
        ! Save data to files.
        !-----------------------------------------------------------------------------------------
        !
        ! Get output format for arrays with nz layers
        !
        WRITE(outfmt, '(a, i3, a)') '(', nz, 'es25.16)'
        WRITE(outfmt_nrbins, '(a, i3, a)') '(', nr_bins, 'es25.16)'

        !
        ! Only save h one time at the beginning
        !
        IF (time == time_start) THEN
            WRITE(12, *) z
        END IF

        WRITE(11, '(f12.8)'   ) time/(24*one_hour)  ! [day], time
        WRITE(13, outfmt) progn%ua                  ! [m s-1], u wind
        WRITE(14, outfmt) progn%va                  ! [m s-1], v wind
        WRITE(15, outfmt) progn%theta               ! [K], potential temperature
        if (output_chemistry) then
            !print *, progn%OH%concentration(2)
            call progn%O3%output
            call progn%O1D%output
            call progn%OH%output
            call progn%REST%output
            call progn%NO2%output
            call progn%NO%output
            call progn%CH2O%output
            call progn%HO2%output
            call progn%CO%output
            call progn%CO2%output
            call progn%CH4%output
            call progn%CH3O2%output
            call progn%isoprene%output
            call progn%RO2%output
            call progn%MVK%output
            call progn%H2O2%output
            call progn%HNO3%output
            call progn%NO3%output
            call progn%N2O5%output
            call progn%SO2%output
            call progn%H2SO4_P%output
            call progn%H2SO4%output
            call progn%alpha_pinene%output
            call progn%HNO3_P%output
            call progn%ELVOC%output
        end if

        write(19, outfmt) progn%PN
        write(20, outfmt) progn%PM
        write(21, outfmt) progn%PV
        write(22, outfmt_nrbins) progn%size_distribution(:,2)


        call parameterizations_output(progn) ! TODO: empty as of yet
    END SUBROUTINE write_files


    !-----------------------------------------------------------------------------------------
    ! Interface
    !
    SUBROUTINE close_files(progn)
        type(prognostics_type), intent(inout) :: progn
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
        close(21)
        close(22)
        if (output_chemistry) then
            call progn%O3%close_file
            call progn%O1D%close_file
            call progn%OH%close_file
            call progn%REST%close_file
            call progn%NO2%close_file
            call progn%NO%close_file
            call progn%CH2O%close_file
            call progn%HO2%close_file
            call progn%CO%close_file
            call progn%CO2%close_file
            call progn%CH4%close_file
            call progn%CH3O2%close_file
            call progn%isoprene%close_file
            call progn%RO2%close_file
            call progn%MVK%close_file
            call progn%H2O2%close_file
            call progn%HNO3%close_file
            call progn%NO3%close_file
            call progn%N2O5%close_file
            call progn%SO2%close_file
            call progn%H2SO4_P%close_file
            call progn%H2SO4%close_file
            call progn%alpha_pinene%close_file
            call progn%HNO3_P%close_file
            call progn%ELVOC%close_file
         end if
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
        !call compute_dynamics(progn, leapfrog_atPreviousFullTime)
        !call compute_parameterizations(progn)

        progn%ua_mid(updInd) = progn%ua(updInd) + 0.5 * dt * progn%dudt
        progn%va_mid(updInd) = progn%va(updInd) + 0.5 * dt * progn%dvdt
        progn%theta_mid(updInd) = progn%theta(updInd) + 0.5 * dt * progn%dThetaDt
        !write(*,*) progn%ua_mid - progn%ua

        ! Initialize chemistry
        call progn%init_chemical_elements()

    end subroutine

    subroutine test_chemistry(progn)
        use radiation_mod
        implicit none
        type(prognostics_type), intent(inout) :: progn
        real(kind = 8) :: PAR_val

        progn%ua = 0.0
        progn%va = 0.0
        progn%theta = 300
        progn%T = 300
        progn%M = 2.4D19
        progn%O2 = 0.21*2.4D19
        progn%N2 = 0.78*2.4D19
        ! Chemical elements already initialized

        PAR_val = -1.0_dp!1000.0 * get_exp_coszen(0.0_dp, daynumber, latitude)
        IF ( time >= time_start_chemistry .and. MOD( NINT((time - time_start)*10.0), NINT(dt_chem*10.0)) == 0 ) THEN
            call compute_chemistry(progn, PAR_val)
        end if


    end subroutine

END PROGRAM main
