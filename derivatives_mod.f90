module derivatives_mod
    use grid_mod, only: nz, z, zmid
    use fft_mod
    implicit none

    contains

    ! function zDeriv(vector)
    ! PURPOSE: Compute height derivative of a vector
    ! INPUT:
    !   real*8    vector(nz) : Vector to differentiate.
    ! OUTPUT:
    !   real*8    zDeriv(nz-1) : Height derivative of the vector (One value shorter than vector)
    function zDeriv(vector)
        implicit none
        real(kind = 8), intent(in) :: vector(:)
        real(kind = 8), allocatable :: zDeriv(:)

        zDeriv = (vector(2:nz) - vector(1:nz-1)) / (z(2:nz) - z(1:nz-1))
    end function zDeriv

    ! function zDerivMidlevel(vector)
    ! PURPOSE: Compute height derivative on midlevel points (mainly for Km, Kh)
    ! INPUT:
    !   real*8    vector(nz-1) : Vector to differentiate.
    ! OUTPUT:
    !   real*8    zDerivMidlevel(nz-2) : Height derivative of the vector defined on midlevels (one shorter than vector)
    function zDerivMidlevel(vector)
        implicit none
        real(kind = 8), intent(in) :: vector(:)
        real(kind = 8), allocatable :: zDerivMidlevel(:)
        integer :: N = nz-1 ! For clarity

        zDerivMidlevel = (vector(2:N) - vector(1:N-1)) / (z(3:N+1) - z(1:N-1))
    end function

    function zDeriv_E_tot(vector)
        implicit none
        real(kind = 8), intent(in) :: vector(:)
        real(kind = 8) :: zDeriv_E_tot(nz-1), diff(nz-2), dz
        integer :: N = nz-1 ! For clarity

        !zDeriv_E_tot = deriv_fft(vector, zmid(2) - zmid(1))

        diff = (vector(2:N) - vector(1:N-1)) / (zmid(2:N) - zmid(1:N-1)) ! midlevels 1 to N-1
        dz = zmid(2) - zmid(1)
        zDeriv_E_tot(1) = (-1.5*vector(1) + 2*vector(2) - 0.5*vector(3)) / (zmid(2) - zmid(1))
        zDeriv_E_tot(2) = 0.5 * (diff(1) + diff(2))

        zDeriv_E_tot(3:N-2) = (vector(1:N-4) - 8*vector(2:N-3)&
                           + 8*vector(4:N-1) - vector(5:N)) / (12*dz)

        zDeriv_E_tot(N-1) = 0.5 * (diff(N-2) + diff(N-1))
        zDeriv_E_tot(N) = (0.5*vector(N-2) - 2*vector(N-1) + 1.5*vector(N)) / (zmid(N) - zmid(N-1))
       !zDeriv_E_tot(2:N-1) = 0.5 * (diff(1:N-2) + diff(2:N-1))
!        zDeriv_E_tot(1) = diff(1)
!        zDeriv_E_tot(N) = 0
!        zDeriv_E_tot(2) = 0.5 * (diff(1) + diff(2))
!        zDeriv_E_tot(N-1) = 0.5 * (diff(N-2) + diff(N-1))

    end



end module derivatives_mod
