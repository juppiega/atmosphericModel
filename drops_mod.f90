module drops_mod
    use parameters_mod
    use time_mod
    implicit none
    private
    public drops_init, compute_drops, drop_fall_velocity

    REAL(DP), DIMENSION(n_drop_bins) :: diameter  ,&    ! Diameter of each size bin
        drop_area                      ,&       ! area concentration in each size bin
        drop_volume                      ,&    ! volume concentration in each size bin
        Re ,&
        coag_loss                            ,&    ! coagulation loss rate of particles in each size bin
        v_dep  ,&                                    ! Dry deposition velocity of particles
        P, dm, midpoints, drop_volume_mid, drop_fall_velocity
    logical :: rain(n_drop_bins)
    real(dp), allocatable :: collision_kernel(:,:), Q(:,:)

    real(dp), parameter :: M_w = 18.01, M_s = 58.44, rho_w = 1D3, rho_s = 2160, surf_tens = 76D-3, R_v = 461, R_a = 287
    real(dp), parameter :: D_w = 0.22D-4, K_a = 25D-3, L = 2.5D6, visc_a = 1.3D-5, alpha = 0.2 ! Condensation coeff Raatikainen et al. (2013)

contains


subroutine drops_init(drop_conc)
    implicit none
    real(dp), intent(inout) :: drop_conc(:)
    real(dp) :: max_diameter
    integer :: i, j

    ! Particle diameters between 2D-9 and max_diameter m:
    max_diameter = 3250D-6
    diameter(1)=4D-6
    DO i=2,n_drop_bins
        diameter(i)=diameter(i-1)*(max_diameter/diameter(1))**(1D0/(n_drop_bins-1))
    END DO

    drop_area = pi * diameter**2
    drop_volume = 1D0/6D0 * pi * diameter**3
    drop_conc = 0

    if (.not. allocated(collision_kernel)) allocate(collision_kernel(n_drop_bins, n_drop_bins), Q(n_drop_bins, n_drop_bins))
    collision_kernel = 0
    do i = 2,size(diameter)
        do j = 1,i-1
            collision_kernel(i,j) = coal_eff(diameter(i),diameter(j)) * pi/4 * (diameter(i) + diameter(j))**2 &
                                    * (terminal_velocity(diameter(i)) - terminal_velocity(diameter(j)));
        end do
    end do
    do i = 1,size(diameter)-1
        do j = i+1,size(diameter)
            collision_kernel(i,j) = collision_kernel(j,i);
        end do
    end do

    do i = 1, size(diameter)
        Re(i) = diameter(i) * terminal_velocity(diameter(i)) / visc_a
    end do

    Q = 0
    do j = 1,size(diameter)-1
        do i = j+1,size(diameter)
            Q(i,j) = 145.37/(rho_w*drop_volume(j)) * diameter(j)/diameter(i) * exp(-7*diameter(j)/diameter(i))
        end do
    end do

    P = 2.94E-7 * exp(34*diameter/2)

    midpoints(2:) = exp(0.5*(log(diameter(1:n_drop_bins-1)) + log(diameter(2:))))
    midpoints(1) = midpoints(2) - (midpoints(3) - midpoints(2))
    dm(1:n_drop_bins-1) = pi/6 * (midpoints(2:)**3 - midpoints(1:n_drop_bins-1)**3)
    dm(n_drop_bins) = 0 ! Never used
    drop_volume_mid = pi/6 * midpoints**3

    rain = diameter > 100E-6

    do j = 1, n_drop_bins
        drop_fall_velocity(j) = terminal_velocity(diameter(j))
    end do

end subroutine drops_init

subroutine compute_drops(drop_conc, T, p, saturation, swelled_diameter, N_new_drops, dry_diam, &
                         a, q, N_tot, LWC, r_eff, rain_rate, total_area)
    implicit none
    real(dp), intent(inout) :: drop_conc(:), saturation, q, N_tot, LWC, r_eff, rain_rate, total_area
    real(dp), intent(in) :: T, p, swelled_diameter(:), N_new_drops(:), dry_diam(:), a
    real(dp) :: b(n_drop_bins), saturation_change, dt_cond, dt_coal, dt_breakup, time_cond, q_change, time_tot


    call compute_nucleation(drop_conc, swelled_diameter, N_new_drops, dry_diam, b)

    !print *, drop_conc
    dt_cond = 0.01
    dt_coal = min(10D0, dt_micro)
    dt_breakup = min(10D0, dt_micro)
    !print *, saturation, sum(drop_conc), drop_conc
    !print *, saturation
    !print *, q

    time_tot = 0
    do while (time_tot <= dt_micro + 1E-6)
        time_cond = 0
        do while(time_cond <= dt_breakup + 1E-6)
            call compute_condensation(drop_conc, a, b, T, p, saturation, dt_cond, saturation_change, q_change)
            saturation = saturation + saturation_change
            q = q + q_change
            time_cond = time_cond + dt_cond

        end do

        call compute_coalescence(drop_conc, dt_coal)
        call compute_breakup(drop_conc, dt_breakup)
        time_tot = time_tot + dt_breakup

    !print *, drop_conc
    !print *, sum(drop_conc*drop_volume)

    !print *, q_change


    !print *, drop_conc
    !print *, sum(drop_conc*drop_volume)


    !print *, drop_conc
    !print *, sum(drop_conc*drop_volume)

    end do

    N_tot = sum(drop_conc)
    LWC = sum(drop_volume*drop_conc)*rho_w
    r_eff = sum((diameter/2)**3 * drop_conc) / sum((diameter/2)**2 * drop_conc)
    rain_rate = sum(pack(drop_volume*drop_conc*rho_w*drop_fall_velocity, rain))
    total_area = sum(pi*(diameter/2)**2 * drop_conc)

end subroutine compute_drops

subroutine compute_nucleation(drop_conc, swelled_diameter, N_new_drops, dry_diam, b)
    implicit none
    real(dp), intent(inout) :: drop_conc(:)
    real(dp), intent(in) :: swelled_diameter(:), N_new_drops(:), dry_diam(:)
    real(dp) :: r3_mean, new_drops_this
    real(dp), intent(out) :: b(n_drop_bins)
    integer :: i
    logical :: ind(n_aer_bins)

    do i = 1, n_drop_bins
        if (i == 1) then
            ind = swelled_diameter < diameter(1)
        else
            ind = diameter(i-1)<=swelled_diameter .and. swelled_diameter<diameter(i)
        end if
        new_drops_this = sum(pack(N_new_drops, ind))
        drop_conc(i) = drop_conc(i) + new_drops_this
        r3_mean = sum(pack(N_new_drops*(dry_diam/2)**3, ind))/(new_drops_this + 1D-5)
        b(i) = 2*M_w*rho_s*r3_mean / (rho_w*M_s)
    end do

end subroutine

subroutine compute_condensation(drop_conc, a, b, T, p, saturation, dt_drops, saturation_change, q_change)
    implicit none
    real(dp), intent(inout) :: drop_conc(:)
    real(dp), intent(in) :: a, b(:), T, p, saturation, dt_drops
    real(dp), intent(out) :: saturation_change, q_change
    real(dp) :: es, r_corr, c1, c2(n_drop_bins), S_corr(n_drop_bins), dt_V(n_drop_bins), &
                volume_new(n_drop_bins), conc_half(n_drop_bins), conc_full(n_drop_bins),&
                 x1, x2, f(n_drop_bins), LWC_before, LWC_after, evap_smallest, Flux_left, Flux_right

    integer :: j

    LWC_before = sum(drop_volume*drop_conc)*rho_w

    es = 288513966.0*exp(-4302.645/(T-29.65)) * 1D2 ! [Pa]
    where (Re < 2.5)
        f = 1 + 0.09*Re
    elsewhere
        f = 0.78+0.28*sqrt(Re)
    end where
    c1 = D_w*L*L*es/(K_a*R_v*R_v*T**3)
    c2 = D_w*es*f/(rho_w*R_v*T)
    r_corr = D_w*sqrt(2*PI/(R_v*T))/alpha/(1+c1)
    S_corr = saturation - a/(midpoints/2) + b/(midpoints/2)**3


    dt_V = 0.5*PI*midpoints**2 * f*c2*S_corr/((midpoints/2 + r_corr)*(1+c1))
    dt_V(n_drop_bins) = 0

    volume_new = drop_volume_mid + dt_drops*dt_V

    conc_half = 0
    conc_half(n_drop_bins) = drop_conc(n_drop_bins)

    Flux_left = 0
    do j = 1,n_drop_bins-1
        !x1 = (drop_volume(j+1) - volume_new(j)) / (drop_volume(j+1) - drop_volume(j))
        !x2 = 1 - x1
        if (dt_V(j+1) > 0) then
            Flux_right = (volume_new(j+1)-drop_volume_mid(j+1))/(drop_volume_mid(j+1)-drop_volume_mid(j))*drop_conc(j)
        else if (j < n_drop_bins-1) then
            Flux_right = (volume_new(j+1)-drop_volume_mid(j+1))/(drop_volume_mid(j+2)-drop_volume_mid(j+1))*drop_conc(j+1)
        else
            Flux_right = 0
        end if

        conc_half(j) = drop_conc(j) - (Flux_right - Flux_left)
        Flux_left = Flux_right
    end do

!    do j = 1,n_drop_bins-1
!        x1 = (drop_volume(j+1) - volume_new(j)) / (drop_volume(j+1) - drop_volume(j))
!        x2 = 1 - x1
!        conc_full(j) = conc_full(j) + x1*conc_half(j)
!        conc_full(j+1) = conc_full(j+1) + x2*conc_half(j)
!    end do
!
!    !drop_conc = drop_conc + (conc_full - conc_half)

    if (conc_half(1) > drop_conc(1)) then
        evap_smallest = (conc_half(1)-drop_conc(1))*rho_w*drop_volume(1)
    else
        evap_smallest = 0
        drop_conc(1) = conc_half(1)
    end if
    drop_conc(2:) = conc_half(2:)

    LWC_after = sum(drop_volume*drop_conc)*rho_w
    q_change = R_a*T * (LWC_before - LWC_after + evap_smallest) / p
    saturation_change = q_change * p / (0.622*es)

    !print *, saturation, saturation_change, drop_conc

    if (any(drop_conc < 0)) then
        print *, drop_conc
        stop 'Droplet concentration negative after condensation'
    end if

end subroutine

subroutine compute_coalescence(drop_conc, dt_drops)
    implicit none
    real(dp), intent(inout) :: drop_conc(:)
    real(dp), intent(in) :: dt_drops
    real(dp), dimension(n_drop_bins) :: dt_V, volume_new, coag_loss, conc_new
    real(dp) :: x1, x2
    integer :: i, j

    !print *, drop_conc
    !print *, sum(drop_volume*drop_conc)

    do i = 1, n_drop_bins
        dt_V(i) = sum(drop_volume(1:i) * drop_conc(1:i) * collision_kernel(1:i,i))
        coag_loss(i) = drop_conc(i) * sum(collision_kernel(i+1:,i) * drop_conc(i+1:)) ! Assumes no self-coagulation
    end do

    volume_new = drop_volume + dt_drops*dt_V
    drop_conc = drop_conc - dt_drops*coag_loss

    conc_new = 0
    conc_new(n_drop_bins) = drop_conc(n_drop_bins)

    do j = 1,n_drop_bins-1
        x1 = (drop_volume(j+1) - volume_new(j)) / (drop_volume(j+1) - drop_volume(j))
        x2 = 1 - x1
        conc_new(j) = conc_new(j) + x1*drop_conc(j)
        conc_new(j+1) = conc_new(j+1) + x2*drop_conc(j)
    end do

    !print *, sum(drop_conc)
    !print *, sum(conc_new)
    !print *, sum(drop_volume*conc_new)
    !print *, conc_new
    !print *, coag_loss
    !print *, collision_kernel(:,n_drop_bins-10)
    !print *, diameter
    !stop 'f'

    drop_conc = conc_new

end subroutine

subroutine compute_breakup(drop_conc, dt_drops)
    implicit none
    real(dp), intent(inout) :: drop_conc(:)
    real(dp), intent(in) :: dt_drops
    real(dp) :: dt_N(n_drop_bins)
    integer :: i

    do i = 1, n_drop_bins-1
        dt_N(i) = -P(i)*drop_conc(i) + dm(i)*sum(P(i+1:)*Q(i+1:,i)*drop_conc(i+1:))
    end do
    dt_N(n_drop_bins) = -P(n_drop_bins)*drop_conc(n_drop_bins)

    drop_conc = drop_conc + dt_N * dt_drops

end subroutine

function terminal_velocity(D_p)
    implicit none
    real(dp), intent(in) :: D_p ! Drop diameter in meters
    real(dp) :: terminal_velocity

    if (D_p < 60E-6) then
        terminal_velocity = 2.98D7*D_p**2;
    else if (D_p < 1.2E-3) then
        terminal_velocity = 26462.62*D_p**1.2772;
    else
        terminal_velocity = 142.1*sqrt(D_p);
    end if

end function

function coal_eff(D_big, d_small) result(e)
    implicit none
    real(dp), intent(in) :: D_big, d_small
    real(dp) :: e, R_big, r_small, y, CKE, S_C, E_T, e_small, e_big, Dmax, Dmin

    R_big = 100*D_big/2; r_small = 100*d_small/2;

    e_small = min(max(4.5E4*R_big**2 * (1 - 3E-4/r_small),0D0),1D0);

    R_big = R_big/100; r_small = r_small / 100;
    y = R_big / r_small;
    CKE = 2D0/3D0*pi*rho_w*(R_big)**3 * (terminal_velocity(2*R_big)-terminal_velocity(2*r_small))**2 / (1+y**3);
    S_C = 4*pi*surf_tens*(r_small**3 + R_big**3)**(2D0/3D0);
    E_T = CKE + 4*pi*surf_tens*(r_small**2 + R_big**2) - S_C;
    if (E_T < 5E-6) then
        e_big = 0.778*(y/(y+1))**2 * exp(-2.61E6*surf_tens*E_T**2/S_C);
    else
        e_big = 0;
    end if

    Dmax = 500E-6; Dmin = 100e-6;
    if (D_big < Dmin) then
        e = e_small;
    else if (D_big > Dmax) then
        e = e_big;
    else
        e = e_small*cos(pi/2 * (D_big-Dmin)/(Dmax-Dmin))**2 + e_big*sin(pi/2 * (D_big-Dmin)/(Dmax-Dmin))**2;
    end if

end function

end module drops_mod
