module derivatives_mod
    use grid_mod, only: nz, z
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

end module derivatives_mod
