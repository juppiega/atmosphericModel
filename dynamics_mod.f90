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
    REAL(dp), PARAMETER :: lambda_ = 200.0_dp  ! maximum mixing length, meters
    REAL(dp), PARAMETER :: vonk_ = 0.4_dp      ! von Karman constant, dimensionless
    real(dp), parameter :: C_f = 0.185, C_n = 2.0, C_eps = (0.17)**1.5, C_r = 2.0
    real(dp), parameter :: Pr_0 = (0.17)**2/(2*(0.145)**2)
    real(dp), parameter :: hd_div = 2D0

    ! Variables for computing the boundary layer turbulent fluxes.
    real(kind = 8), dimension(nz-1) :: richardsonNum_, Km_, Kh_, &
        verticalWindShear_, mixingLength_, &
        dThetaDz_, dUDz_, DVDz_, DEDz_, dqdz_, E_tot_, E_pot_, E_kin_, N_squared_, &
        flux_uw_, flux_vw_, flux_ThetaW, flux_Etot_w_, flux_qw_, flux_qw_k, f_tau_, f_theta_,&
        dissipation_, shear_prod_, transport_, buoyancy_, E_u_

    real(kind = 8), dimension(nz) :: u_u_, v_u_, q_u_, theta_u_, w_kin_, mass_flux_, detrainment_
    real(kind = 8), dimension(nz) :: u_, v_, q_, theta_ ! TODO: Pointtereiksi, jotka osoittaa progn:n oikeisiin taulukoihin
    real(kind = 8) :: first_level_Ep_E_tot_, w_star_surf_, hd_, entr_


CONTAINS

    ! function compute_verticalWindShear() result(vws)
    ! PURPOSE: Compute vertical wind shear.
    ! OUTPUT:
    !     (real*8) : vws [vertical wind shear]
    function compute_verticalWindShear() result(vws)
        implicit none
        real(kind = 8), allocatable :: vws(:) ! Vertical wind shear

        vws = sqrt(dUDz_**2 + DVDz_**2)

    end function

    subroutine compute_richardsonNum_and_N()
        implicit none
        real(kind = 8), allocatable :: Ri(:), theta_halfLevel(:)

        theta_halfLevel = (theta_(2:nz) + theta_(1:nz-1)) / 2.0
        N_squared_ = g * dThetaDz_ / theta_halfLevel

        richardsonNum_ = (N_squared_) / (verticalWindShear_**2)

    end subroutine

    ! function compute_mixingLength() result(mixingLength)
    ! PURPOSE: Compute mixing length.
    ! OUTPUT:
    !     (real*8) : mixingLength
    function compute_mixingLength() result(mixingLength)
        implicit none
        real(kind = 8), allocatable :: mixingLength(:)

        mixingLength = (vonk_ * z(2:size(z))) / (1 + vonk_ * z(2:size(z)) / lambda_)

    end function

    ! function compute_mixingLength_TEMF() result(mixingLength)
    ! PURPOSE: Compute mixing length.
    ! OUTPUT:
    !     (real*8) : mixingLength
    subroutine compute_mixingLength_TEMF()
        implicit none
        real(kind = 8), allocatable :: reciproc(:)

        reciproc = 1 / (vonk_ * 0.5*(z(1:size(z)-1) + z(2:size(z)))) + &
            abs(fcor) / (C_f * sqrt(f_tau_*E_kin_)) + &
            1.0/lambda_
        where (N_squared_ > 0)
            reciproc = reciproc + sqrt(N_squared_) / (C_n * sqrt(f_tau_*E_kin_))
        end where

        mixingLength_ = 1 / reciproc

    end subroutine

    subroutine compute_normalized_stresses()
        implicit none

        ! Neutral values where unstable
        f_tau_ = 0.17
        f_theta_ = -0.145
        where (richardsonNum_ > 0)
            f_tau_ = 0.17 * (0.25 + 0.75 / (1 + 4*richardsonNum_))
            f_theta_ = -0.145 / (1 + 4*richardsonNum_)
        end where
    end subroutine

    subroutine compute_Ekin_Epot()
        implicit none
        real(kind = 8) :: E_pot_E_tot_ratio(nz-1)

        where (richardsonNum_ >= 0)
            E_pot_E_tot_ratio = richardsonNum_ / (Pr_0 + 3*richardsonNum_)
        elsewhere
            E_pot_E_tot_ratio = richardsonNum_ / (2*richardsonNum_ - Pr_0)
        end where

        E_pot_ = E_tot_ * E_pot_E_tot_ratio
        E_kin_ = E_tot_ - E_pot_
        first_level_Ep_E_tot_ = E_pot_E_tot_ratio(1)

    end subroutine

    ! function compute_Km() result(Km)
    ! PURPOSE: Compute turbulent coefficient for momentum transfer.
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

    ! function compute_Km_TEMF() result(Km)
    ! PURPOSE: Compute turbulent coefficient for momentum transfer.
    function compute_Km_TEMF() result(Km)
        implicit none
        real(kind = 8) :: theta_variance(nz-1), theta_mid(nz-1)
        real(kind = 8), allocatable :: Km(:), Km_conv(:), reciproc(:), mixingLength_conv(:)
        integer :: i

        theta_mid = (theta_(2:nz) + theta_(1:nz-1)) / 2.0
        theta_variance = max(2*E_pot_*N_squared_*(theta_mid / g)**2, 0.0D0)
        Km = f_tau_**2 * E_kin_**2 / (C_eps * E_kin_*sqrt(E_tot_)/mixingLength_ - &
            g*f_theta_*sqrt(E_kin_*theta_variance)/theta_mid)

        if (flux_ThetaW(1) > 0) then
            reciproc = 1 / (vonk_ * zmid) + 3 / (vonk_ * (hd_ - zmid))
            mixingLength_conv = 1 / reciproc
            Km_conv = (0.17**2)*mixingLength_conv*sqrt(E_kin_)/C_eps
            i = 1
            do while (zmid(i) <= hd_)
                if (zmid(i) <= hd_/hd_div .or. Km_conv(i) > Km(i)) then
                    Km(i) = Km_conv(i)
                else
                    Km(i) = ((zmid(i) - hd_/hd_div)*Km(i) + (hd_ - zmid(i))*Km_conv(i)) / (hd_ - hd_/hd_div)
                end if
                i = i + 1
            end do
        end if

        Km = max(Km, 1.57D-4)

    end function

    ! function compute_Kh() result(Kh)
    ! PURPOSE: Compute turbulent coefficient for heat and chemical transport.
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

    ! function compute_Kh_TEMF() result(Kh)
    ! PURPOSE: Compute turbulent coefficient for heat and chemical transport.
    function compute_Kh_TEMF() result(Kh)
        implicit none
        real(kind = 8), allocatable :: Kh(:), Kh_conv(:), reciproc(:), mixingLength_conv(:)
        integer :: i

        Kh = 2*F_theta_**2 * E_kin_ * mixingLength_ / (C_eps * sqrt(E_tot_))

        if (flux_ThetaW(1) > 0) then
            reciproc = 1 / (vonk_ * zmid) + 3 / (vonk_ * (hd_ - zmid))
            mixingLength_conv = 1 / reciproc
            Kh_conv = (0.17**2)*mixingLength_conv*sqrt(E_kin_)/C_eps/Pr_0
            i = 1
            do while (zmid(i) <= hd_)
                if (zmid(i) <= hd_/hd_div .or. Kh_conv(i) > Kh(i)) then
                    Kh(i) = Kh_conv(i)
                    mixingLength_(i) = mixingLength_conv(i)
                else
                    Kh(i) = ((zmid(i) - hd_/hd_div)*Kh(i) + (hd_ - zmid(i))*Kh_conv(i)) / (hd_ - hd_/hd_div)
                end if
                i = i + 1
            end do
        end if

        Kh = max(Kh, 1.57D-4 / 0.733)

    end function

    subroutine compute_e_flux()
        implicit none
        real(kind = 8) :: M_mid(nz-1)
        integer :: i

        flux_Etot_w_ = -Km_ * DEDz_
        if (flux_ThetaW(1) > 0) then
            M_mid = 0
            i = 1
            do while (zmid(i) < hd_)
                M_mid(i) = interp1_scalar(z, mass_flux_, zmid(i))
                i = i + 1
            end do
            flux_Etot_w_ = flux_Etot_w_ + M_mid * (E_u_ - E_tot_)
        end if

    end subroutine

    subroutine compute_fluxes(flux, K, gradient, updraft, environment)
        implicit none
        real(kind = 8), intent(in) :: K(:), gradient(:), updraft(:), environment(:)
        real(kind = 8), intent(inout) :: flux(:)
        real(kind = 8) :: midlevel_updraft(nz-1), updraft_flux(nz), midlevel_env(nz-1)
        integer :: i

        flux(2:size(flux)) = -K(2:size(flux)) * gradient(2:size(flux))
        if (flux_ThetaW(1) > 0) then
            updraft_flux = mass_flux_ * (updraft - environment)
            i = 1
            midlevel_updraft = 0
            do while(zmid(i) < hd_)
                midlevel_updraft(i) = interp1_scalar(z, updraft_flux, zmid(i))
                i = i + 1
            end do

            flux = flux + midlevel_updraft
        end if

    end subroutine

    subroutine compute_near_surface_turbulence()
        implicit none
        real(kind = 8) :: zt = 0.5*(z0 + z(2)), theta_t, &
            tau, factor, m, b, eps, buoyancy_first, flux_first
        integer, parameter :: maxIter = 1000

        theta_t = 0.5*(theta_(1) + theta_(2))
        flux_uw_(1) = -u_(2) / (zt * log(z(2)/z0)) * (f_tau_(1) / 0.17) * mixingLength_(1)
        flux_vw_(1) = -v_(2) / (zt * log(z(2)/z0)) * (f_tau_(1) / 0.17) * mixingLength_(1)

        tau = sqrt(flux_uw_(1)**2 + flux_vw_(1)**2)
        factor = (theta_(2) - theta_(1)) / (zt * log(z(2)/z0)) * (f_theta_(1) / 0.145)&
            * mixingLength_(1) / Pr_0
        m = g * hd_ / (theta_t * 5)
        b = tau**(1.5)
        eps = 1D-7

        !print *, '  Passed wind_stresses'
        !print *, theta_(2), E_tot_(1)
        if (theta_(2) < theta_(1)) then
            flux_first = regulaFalsi()
            if (flux_first < -1D8) then
                !eps = 1D-5
                flux_first = bisection()
                if (flux_first < -1D8) then
                    print *, factor, m, b, theta_(2)-theta_(1)
                    stop 'Could not find heat flux'
                end if
                !eps = 1D-7
                !print *, flux_first
            end if
            flux_ThetaW(1) = flux_first
        else
            flux_ThetaW(1) = factor * sqrt(tau)
        end if

        w_star_surf_ = (5*m*flux_ThetaW(1))**(1.0/3.0)
        factor = (q_(2) - q_(1)) / (zt * log(z(2)/z0)) * (f_theta_(1) / 0.145)&
            * mixingLength_(1) / Pr_0

        if (theta_(2) < theta_(1)) then
            flux_qw_(1) = factor * (b + 0.2*w_star_surf_**3)**(1.0/3.0)
        else
            flux_qw_(1) = factor * sqrt(tau)
        end if

        buoyancy_first = 2 * g * flux_ThetaW(1) / theta_t

        !print*, tau, mixingLength_(1), buoyancy_first

        E_tot_(1) = (1 + first_level_Ep_E_tot_)*(tau**(3.0/2.0) + mixingLength_(1)*buoyancy_first)**(2.0/3.0)/0.17
        E_pot_(1) = first_level_Ep_E_tot_ * E_tot_(1)
        E_kin_(1) = E_tot_(1) - E_pot_(1)
        !print *, '  Passed B and temp flux'

        !print *, E_tot_(1), E_tot_(2), first_level_Ep_E_tot_

    contains

        function f(x)
            real(kind = 8) :: x,f
            f = x - factor * (m * x + b)**(1.0/3.0)
        end function

        function regulaFalsi() result(flux)
            real(kind = 8) :: flux, x1, x2, f1, f2, c, fc
            integer :: numIter

            x1 = -b/m
            x2 = 1D4
            f1 = -b/m
            f2 = f(x2)
            numIter = 0

            do while (abs(f(x2)) > eps)
                c = (f2*x1 - f1*x2) / (f2 - f1);
                fc = f(c);
                if (sign(1D0,fc)*sign(1D0,f1) > 0) then
                    x1 = c;
                    f1 = fc;
                else
                    x2 = c;
                    f2 = fc;
                end if
                numIter = numIter + 1
                if (numIter > maxIter) then
                    flux = -1D10
                    return
                end if
            end do

            flux = x2
        end function


        function bisection() result(flux)
            real(kind = 8) :: flux, x1, x2, f1, f2, c, fc
            integer :: numIter

            !print *, 'Regula Falsi failed - resorting to bisection method'
            x1 = -b/m
            x2 = 1D4
            f1 = -b/m
            f2 = f(x2)
            numIter = 0

            do while (abs(f(x2)) > eps)
                c = 0.5D0*(x1+x2);
                fc = f(c);
                if (sign(1D0,fc)*sign(1D0,f1) > 0) then
                    x1 = c;
                    f1 = fc;
                else
                    x2 = c;
                    f2 = fc;
                end if
                numIter = numIter + 1
                if (numIter > maxIter) then
                    flux = -1D10
                    return
                end if
            end do

            flux = x2
        end function

    end subroutine



    ! function compute_E_tot_Tendency() result(dEdt)
    ! PURPOSE: Compute dE/dt.
    function compute_E_tot_Tendency() result(dEdt)
        implicit none
        real(kind = 8), allocatable :: dEdt(:)
        real(kind = 8) :: subsid(nz-1)

        dissipation_ = 0.07 * (E_tot_**(1.5) / mixingLength_)
        shear_prod_ = -flux_uw_*dUDz_ - flux_vw_*DVDz_
        transport_ = -zDeriv_E_tot(flux_Etot_w_)
        buoyancy_ = 2 * N_squared_ * flux_ThetaW / dThetaDz_
        subsid = 0.5*(w_subsidence(1:nz-1) + w_subsidence(2:nz)) * DEDz_
        dEdt = shear_prod_ + transport_ - dissipation_ - subsid
        where (N_squared_ < 0)
            dEdt = dEdt + buoyancy_
        end where

    end function

    ! function compute_uTendency() result(dudt)
    ! PURPOSE: Compute du/dt, Equation (1).
    function compute_uTendency() result(dudt)
        implicit none
        real(kind = 8), allocatable :: dudt(:)

        dudt = - zDerivMidLevel(flux_uw_) + fcor * (v_(updInd) - vg) - &
                 w_subsidence(updInd) * 0.5*(dUdz_(1:nz-2) + dUdz_(2:nz-1)) ! - (u_(updInd) - 5)/3600

    end function

    ! function compute_vTendency() result(dvdt)
    ! PURPOSE: Compute dv/dt, Equation (2).
    function compute_vTendency() result(dvdt)
        implicit none
        real(kind = 8), allocatable :: dvdt(:)

        dvdt = - zDerivMidLevel(flux_vw_) - fcor * (u_(updInd) - ug) - &
                w_subsidence(updInd) * 0.5*(dVdz_(1:nz-2) + dVdz_(2:nz-1))  ! - v_(updInd)/3600

    end function

    ! function compute_thetaTendency() result(dThetaDt)
    ! PURPOSE: Compute d Theta / dt. Equation (3).
    function compute_thetaTendency() result(dThetaDt)
        implicit none
        real(kind = 8), allocatable :: dThetaDt(:)

        dThetaDt = -zDerivMidLevel(flux_ThetaW) - w_subsidence(updInd) * 0.5*(dThetaDz_(1:nz-2) + dThetaDz_(2:nz-1)) &
                    -1.3/86400D0

    end function

    function compute_qTendency() result(dqDt)
        implicit none
        real(kind = 8), allocatable :: dqDt(:)

        dqDt = -zDerivMidLevel(flux_qw_) - w_subsidence(updInd) * 0.5*(dqdz_(1:nz-2) + dqdz_(2:nz-1))

    end function

    function compute_updraft(var) result(var_updraft)
        implicit none
        real(kind = 8), intent(in) :: var(nz)
        real(kind = 8) :: var_updraft(nz), integral(nz)
        real(kind = 8) :: var_u_mid, var_interp
        integer :: i

        !integral = cumtrapz(z, exp(entr_ * z) * var)
        !var_updraft = entr_ * exp(-entr_ * z) * integral + var(1) * exp(-entr_ * z)

        var_updraft = var
        var_updraft(1) = var(1)
        var_u_mid = var_updraft(1)
        var_interp = log(zmid(1)/z0) * (var(2) - var(1)) / log(z(2)/z0) + var(1)
        i = 2

        do while (i <= nz)
            var_updraft(i) = var_updraft(i-1) - (entr_ * (var_u_mid - var_interp))*(z(i)-z(i-1))
            if (i < nz) var_u_mid = var_u_mid - (entr_ * (var_updraft(i) - var(i))) * (zmid(i)-zmid(i-1))
            if (i < nz) var_interp = interp1_scalar(z, var, zmid(i))
            i = i + 1
        end do


    end function

    subroutine compute_E_u()
        implicit none
        real(kind = 8) :: integral(nz-1)
        integral = cumtrapz(zmid, exp(entr_* zmid) * E_tot_)
        E_u_ = entr_ * exp(-entr_*zmid) * integral + E_tot_(1) * exp(-entr_ * zmid)
    end subroutine

    subroutine compute_hd()
        implicit none
        integer :: i
        real(kind = 8) :: B(nz), buoyancy_drag(nz), w_kin_init

        B = g * (theta_u_ - theta_) / (theta_ * 3)
        buoyancy_drag = exp(-4*entr_*z) * cumtrapz(z, B*exp(4*entr_*z))
        w_kin_init = 0.5 * (0.5 * w_star_surf_)**2
        w_kin_ = buoyancy_drag + w_kin_init * exp(-4*entr_*z)

        do i = 1, nz
            if (w_kin_(i) < 0) then
                hd_ = max(((z(i-1)*w_kin_(i) - z(i)*w_kin_(i-1)) / (w_kin_(i)-w_kin_(i-1))), 100.1D0)
                !print *, hd_
                exit
            end if
        end do

    end subroutine

    function cumtrapz(x, y)
        implicit none
        real(kind = 8), intent(in) :: x(:), y(:)
        real(kind = 8) :: cumtrapz(size(x))
        integer :: i

        cumtrapz(1) = 0
        do i = 2, size(x)
            cumtrapz(i) = cumtrapz(i-1) + 0.5 * (x(i) - x(i-1)) * (y(i) + y(i-1))
        end do

    end function

    function interp1_scalar(x, Y, xq)
        implicit none
        real(kind = 8), intent(in) :: x(:), Y(:), xq
        real(kind = 8) :: interp1_scalar
        integer :: i,n

        n = 1
        do i = 2,size(x)
            if (x(i) >= xq) then
                n = i-1
                exit
            end if
        end do

        interp1_scalar = (Y(n+1) - Y(n)) / (x(n+1) - x(n)) * (xq - x(n)) + Y(n)

    end function

    subroutine compute_mass_flux()
        implicit none
        real(kind = 8) :: detrainment_integral(nz), flux_first_level

        detrainment_integral = cumtrapz(z, detrainment_)
        flux_first_level = 0.03 * w_star_surf_
        mass_flux_ = flux_first_level * exp(entr_*z - detrainment_integral)
        where (z > hd_)
            mass_flux_ = 0
        end where

    end subroutine

    ! subroutine compute_dynamics(progn, stepType)
    ! PURPOSE: Compute dynamics (tendencies from BL turbulence).
    ! INPUT:
    !   (prognostics_type) : progn [contains the full atmospheric state]
    subroutine compute_dynamics(progn, stepType)
        implicit none
        type(prognostics_type), intent(inout) :: progn
        integer, intent(in) :: stepType
        CHARACTER(255) :: outfmt

        w_kin_ = 0
        mass_flux_ = 0

        ! Save copies to member variables of dynamics_mod.
        u_ = progn%ua
        v_ = progn%va
        q_ = progn%q
        theta_ = progn%theta
        E_tot_ = progn%E_tot
        hd_ = progn%hd

        ! Compute derivatives and vertical shear: d(sqrt(u**2 + v**2))/dz
        dUDz_ = zDeriv(u_)
        DVDz_ = zDeriv(v_)
        dqdz_ = zDeriv(q_)
        dThetaDz_ = zDeriv(theta_)
        dThetaDz_(nz-1) = 3D-3
        verticalWindShear_ = compute_verticalWindShear()

        ! Richardson number
        call compute_richardsonNum_and_N()
        !print *, 'Passed RI_N'

        if (scheme == 2) then
            call compute_normalized_stresses()
            !print *, 'Passed norm. stress'
            call compute_Ekin_Epot()
            !print *, 'Passed Es'
            call compute_mixingLength_TEMF()
            !print *, 'Passed mixLen'
            call compute_near_surface_turbulence()
            DEDz_ = zDeriv_E_tot(E_tot_)
            !print *, 'Passed near_surf_turb'
        end if


        ! Turbulent coefficients
        if (scheme == 1) then
            Km_ =  compute_Km()
            Kh_ =  compute_Kh()
            flux_uw_ = -Km_ * dUDz_
            flux_vw_ = -Km_ * dVDz_
            flux_ThetaW = -Kh_ * dThetaDz_
        else if (scheme == 2) then
            if (flux_ThetaW(1) > 0) then
                theta_u_ = compute_updraft(theta_)
                call compute_hd()
                entr_ = C_r / max(100.0D0, hd_)
                u_u_ = compute_updraft(u_)
                v_u_ = compute_updraft(v_)
                q_u_ = compute_updraft(q_)
                call compute_E_u()
                detrainment_ = max(0D0, entr_ + 0.05 * (1 / (hd_ - z))) ! Risk of division by zero
                call compute_mass_flux()
            end if

            Km_ =  compute_Km_TEMF()
            Kh_ =  compute_Kh_TEMF()
            !print *, 'Passed KmKh'
            !theta_u_ = 0 !!!! TESTAUS
            call compute_fluxes(flux_uw_, Km_, dUDz_, u_u_, u_)
            call compute_fluxes(flux_vw_, Km_, dVDz_, v_u_, v_)
            call compute_E_flux()
            !dqdz_ = 0
            flux_qw_k(2:nz-1) = -Kh_(2:nz-1) * dqdz_(2:nz-1)
            flux_qw_k(1) = flux_qw_(1)
            if (time < 4.5*86400) then
                call compute_fluxes(flux_qw_, Kh_, dqdz_, q_u_, q_)
            else
                call compute_fluxes(flux_qw_, Kh_, dqdz_, q_, q_)
            end if
            call compute_fluxes(flux_ThetaW, Kh_, dThetaDz_, theta_u_, theta_)
        end if

        progn%DEDt = compute_E_tot_Tendency()
        !mixingLength_ = (vonk_ * z(2:size(z))) / (1 + vonk_ * z(2:size(z)) / lambda_)
        !Km_ =  compute_Km()
        !Kh_ =  compute_Kh()
        !flux_uw_ = -Km_ * dUDz_
        !flux_vw_ = -Km_ * DVDz_
        !flux_ThetaW = -Kh_ * dThetaDz_
        ! Compute tendencies for meteorology prognostics

        progn%dudt = compute_uTendency() ! Equation (1)
        progn%dvdt = compute_vTendency() ! Equation (2)
        progn%dThetaDt = compute_thetaTendency() ! Equation (3)
        progn%dqdt = compute_qTendency()
        !progn%dqdt(2:nz-2-1) = (progn%dqdt(1:nz-2-2) + progn%dqdt(2:nz-2-1) + progn%dqdt(3:nz-2)) / 3

        progn%E_tot(1) = E_tot_(1)
        progn%hd = hd_

        !print *, progn%dudt(1), progn%dvdt(1), progn%dThetaDt(1), progn%DEDt(2)
        !
        !
        !    ! Compute chemistry dynamical tendencies every timestep, but the parameterizations (in parameterizations_mod) only every dt_chem
        !    call progn % O3     % compute_dynamical_tendency(Kh_)
        !    call progn % O1D     % compute_dynamical_tendency(Kh_)
        !    call progn % OH     % compute_dynamical_tendency(Kh_)
        !    call progn % REST     % compute_dynamical_tendency(Kh_)
        !    call progn % NO2     % compute_dynamical_tendency(Kh_)
        !    call progn % NO     % compute_dynamical_tendency(Kh_)
        !    call progn % CH2O     % compute_dynamical_tendency(Kh_)
        !    call progn % HO2     % compute_dynamical_tendency(Kh_)
        !    call progn % CO     % compute_dynamical_tendency(Kh_)
        !    call progn % CO2     % compute_dynamical_tendency(Kh_)
        !    call progn % CH4     % compute_dynamical_tendency(Kh_)
        !    call progn % CH3O2     % compute_dynamical_tendency(Kh_)
        !    call progn % isoprene     % compute_dynamical_tendency(Kh_)
        !    call progn % RO2     % compute_dynamical_tendency(Kh_)
        !    call progn % MVK     % compute_dynamical_tendency(Kh_)
        !    call progn % H2O2     % compute_dynamical_tendency(Kh_)
        !    call progn % HNO3     % compute_dynamical_tendency(Kh_)
        !    call progn % NO3     % compute_dynamical_tendency(Kh_)
        !    call progn % N2O5     % compute_dynamical_tendency(Kh_)
        !    call progn % SO2     % compute_dynamical_tendency(Kh_)
        !    call progn % H2SO4     % compute_dynamical_tendency(Kh_)
        !    call progn % H2SO4_P     % compute_dynamical_tendency(Kh_)
        !    call progn % alpha_pinene     % compute_dynamical_tendency(Kh_)
        !    call progn % HNO3_P     % compute_dynamical_tendency(Kh_)
        !    call progn % ELVOC     % compute_dynamical_tendency(Kh_)
        !

        ! Output Km, Kh, Ri to files for diagnostics.
        IF ( time >= time_start_output .and. MOD( NINT((time - time_start)*10.0), NINT(dt_output*10.0)) == 0 ) THEN
            WRITE(outfmt, '(a, i3, a)') '(', nz-1, 'es25.16)'
            !        WRITE(*, '(a8, f8.3, a6)') 'time = ', time/one_hour, '  hours'
            WRITE(16, outfmt) Km_
            WRITE(17, outfmt) Kh_
            WRITE(18, outfmt) richardsonNum_
            WRITE(23, outfmt) E_tot_
            WRITE(24, outfmt) dissipation_
            WRITE(25, outfmt) transport_
            WRITE(26, outfmt) shear_prod_
            WRITE(27, outfmt) mixingLength_
            WRITE(28, outfmt) flux_uw_
            WRITE(29, outfmt) flux_vw_
            WRITE(31, outfmt) flux_ThetaW
            WRITE(32, outfmt) flux_Etot_w_
            WRITE(33, outfmt) E_kin_
            WRITE(34, outfmt) E_pot_
            WRITE(35, outfmt) buoyancy_
            WRITE(36, outfmt) theta_u_(1:nz-1)
            WRITE(38, outfmt) w_kin_(1:nz-1)
            WRITE(39, *) hd_
            WRITE(41, outfmt) flux_qw_
            WRITE(45, outfmt) flux_qw_k
            WRITE(outfmt, '(a, i3, a)') '(', nz, 'es25.16)'
            WRITE(46, outfmt) mass_flux_
        END IF

    end subroutine


    !-----------------------------------------------------------------------------------------
    ! Interface
    !
    SUBROUTINE dynamics_init()
        implicit none
        CHARACTER(255), PARAMETER :: outdir = 'output'
        mixingLength_ = compute_mixingLength() ! SAILYTA
    END SUBROUTINE dynamics_init




END MODULE dynamics_mod

