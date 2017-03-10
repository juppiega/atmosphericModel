!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Dynamics module
! - Compute eddy diffusivity coefficients
! - Compute wind fields
! - Compute theta field
! - Mix chemicals and aerosols
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE dynamics_mod

!-----------------------------------------------------------------------------------------
! Load modules
!-----------------------------------------------------------------------------------------
USE parameters_mod
USE time_mod
USE grid_mod
use derivatives_mod
use prognostics_mod

!-----------------------------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------------------------
IMPLICIT NONE  ! This applies in all subprograms inside the module

PRIVATE
PUBLIC dynamics_init
PUBLIC compute_dynamics  ! functions

!
! Some constants
!
REAL(dp), PARAMETER :: lambda_ = 300.0_dp  ! maximum mixing length, meters
REAL(dp), PARAMETER :: vonk_ = 0.4_dp      ! von Karman constant, dimensionless

real(kind = 8), dimension(nz-1) :: richardsonNum_, Km_, Kh_, &
                  verticalWindShear_, mixingLength_, &
                  dThetaDz_, dUDz_, DVDz_
real(kind = 8), dimension(nz) :: u_, v_, theta_ ! TODO: Pointtereiksi, jotka osoittaa progn:n oikeisiin taulukoihin


CONTAINS

function compute_verticalWindShear() result(vws)
    implicit none
    real(kind = 8), allocatable :: vws(:) ! Vertical wind shear

    vws = sqrt(dUDz_**2 + DVDz_**2)

end function

function compute_richardsonNum() result(Ri)
    implicit none
    real(kind = 8), allocatable :: Ri(:), theta_halfLevel(:)

    theta_halfLevel = (theta_(2:nz) + theta_(1:nz-1)) / 2.0
    Ri = (g * dThetaDz_ / theta_halfLevel) / (verticalWindShear_**2)

end function

function compute_mixingLength() result(mixingLength) ! Aja initissa
    implicit none
    real(kind = 8), allocatable :: mixingLength(:)

    mixingLength = (vonk_ * z(2:size(z))) / (1 + vonk_ * z(2:size(z)) / lambda_)

end function

function compute_Km() result(Km)
    implicit none
    real(kind = 8) :: f(nz-1) ! Function of Richardson number and altitude
    real(kind = 8), allocatable :: Km(:)

    f = 0.1
    where (0 <= richardsonNum_ .and. richardsonNum_ < 0.2)
        f = max((1 - 5*richardsonNum_)**2, 0.1)
    end where
    where (richardsonNum_ < 0)
        f = (1 - 16*richardsonNum_)**0.5
    end where

    Km = mixingLength_**2 * verticalWindShear_ * f

end function

function compute_Kh() result(Kh)
    implicit none
    real(kind = 8) :: f(nz-1) ! Function of Richardson number and altitude
    real(kind = 8), allocatable :: Kh(:)

    f = 0.1
    where (0 <= richardsonNum_ .and. richardsonNum_ < 0.2)
        f = max((1 - 5*richardsonNum_)**2, 0.1)
    end where
    where (richardsonNum_ < 0)
        f = (1 - 16*richardsonNum_)**0.75
    end where

    Kh = mixingLength_**2 * verticalWindShear_ * f

end function

function compute_uTendency() result(dudt)
    implicit none
    real(kind = 8), allocatable :: dudt(:)

    dudt = fcor * (v_(updInd) - vg) + zDerivMidLevel(Km_ * dUDz_)

end function

function compute_vTendency() result(dvdt)
    implicit none
    real(kind = 8), allocatable :: dvdt(:)

    dvdt = -fcor * (u_(updInd) - ug) + zDerivMidLevel(Km_ * DVDz_)

end function

function compute_thetaTendency() result(dThetaDt)
    implicit none
    real(kind = 8), allocatable :: dThetaDt(:)

    dThetaDt = zDerivMidLevel(Kh_ * dThetaDz_)

end function

subroutine compute_dynamics(progn, stepType)
    implicit none
    type(prognostics_type), intent(inout) :: progn
    integer, intent(in) :: stepType
    CHARACTER(255) :: outfmt

    ! Save copies to member variables of dynamics_mod.
    u_ = progn%ua
    v_ = progn%va
    theta_ = progn%theta

    ! Compute derivatives and vertical shear: d(sqrt(u**2 + v**2))/dz
    dUDz_ = zDeriv(u_)
    DVDz_ = zDeriv(v_)
    dThetaDz_ = zDeriv(theta_)
    verticalWindShear_ = compute_verticalWindShear()

    ! Richardson number
    richardsonNum_ = compute_richardsonNum()

    ! Turbulent coefficients
    Km_ =  compute_Km()
    Kh_ =  compute_Kh()

    ! Compute tendencies for meteorology prognostics
    progn%dudt = compute_uTendency() ! Equation (1)
    progn%dvdt = compute_vTendency() ! Equation (2)
    progn%dThetaDt = compute_thetaTendency() ! Equation (3)

    ! Compute chemistry dynamical tendencies every timestep, but the parameterizations (in parameterizations_mod) only every dt_chem
    call progn % alpha_pinene % compute_dynamical_tendency(Kh_)
    call progn % isoprene     % compute_dynamical_tendency(Kh_)

    ! Output Km, Kh, Ri to files for diagnostics.
    IF ( time >= time_start_output .and. MOD( NINT((time - time_start)*10.0), NINT(dt_output*10.0)) == 0 ) THEN
        WRITE(outfmt, '(a, i3, a)') '(', nz-1, 'es25.16)'
!        WRITE(*, '(a8, f8.3, a6)') 'time = ', time/one_hour, '  hours'
        WRITE(16, outfmt) Km_
        WRITE(17, outfmt) Kh_
        WRITE(18, outfmt) richardsonNum_
    END IF

end subroutine


!-----------------------------------------------------------------------------------------
! Interface
!
SUBROUTINE dynamics_init()
    implicit none
    CHARACTER(255), PARAMETER :: outdir = 'output'
    mixingLength_ = compute_mixingLength() ! SAILYTA
    OPEN(16, FILE = TRIM(ADJUSTL(outdir))//'/Km.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
    OPEN(17, FILE = TRIM(ADJUSTL(outdir))//'/Kh.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
    OPEN(18, FILE = TRIM(ADJUSTL(outdir))//'/Ri.dat' , STATUS = 'REPLACE', ACTION = 'WRITE')
END SUBROUTINE dynamics_init




END MODULE dynamics_mod

