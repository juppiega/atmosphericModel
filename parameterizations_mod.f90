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

subroutine compute_parameterizations(progn)
    use time_mod
    implicit none
    type(prognostics_type), intent(inout) :: progn

    ! Compute chemistry parameterizations (emission and deposition) every dt_chem.
    IF ( time >= time_start_chemistry .and. MOD( NINT((time - time_start)*10.0), NINT(dt_chem*10.0)) == 0 ) THEN
        call compute_chemistry(progn)
        !print *, 'Chemistry called'
    end if

end subroutine

subroutine compute_chemistry(progn)
    implicit none
    type(prognostics_type), intent(inout) :: progn

    call compute_radiation() ! PAR
    call compute_aerodynamic_resistance(progn) ! r_a

    ! Compute parameterized tendencies for chemicel components. Q_emission + Q_deposition (equation (20) + (26))
    call progn % alpha_pinene % compute_parameterized_tendency(progn)
    call progn % isoprene     % compute_parameterized_tendency(progn)

end subroutine

subroutine parameterizations_output(progn)
    implicit none
    type(prognostics_type), intent(in) :: progn

end subroutine

end module parameterizations_mod
