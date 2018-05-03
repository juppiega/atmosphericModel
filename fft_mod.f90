
module fft_mod
  use parameters_mod
  implicit none
  !integer,       parameter, private :: dp=selected_real_kind(15,300)
  !real(kind=dp), parameter :: pi=3.141592653589793238460_dp
contains

  ! In place Cooley-Tukey FFT
  recursive subroutine fft(x)
    complex(kind=dp), dimension(:), intent(inout)  :: x
    complex(kind=dp)                               :: t
    integer                                        :: N
    integer                                        :: i
    complex(kind=dp), dimension(:), allocatable    :: even, odd

    N=size(x)

    if(N .le. 1) return

    allocate(odd((N+1)/2))
    allocate(even(N/2))

    ! divide
    odd =x(1:N:2)
    even=x(2:N:2)

    ! conquer
    call fft(odd)
    call fft(even)

    ! combine
    do i=1,N/2
       t=exp(cmplx(0.0_dp,-2.0_dp*pi*real(i-1,dp)/real(N,dp),kind=dp))*even(i)
       x(i)     = odd(i) + t
       x(i+N/2) = odd(i) - t
    end do

    deallocate(odd)
    deallocate(even)

  end subroutine fft

  subroutine ifft(x)
     implicit none
     complex(kind=dp), dimension(:), intent(inout)  :: x

     x = conjg(x)
     call fft(x)
     x = conjg(x)
     x = x / size(x)

  end subroutine

  function deriv_fft(x, dx)
    implicit none
    real(kind = 8), intent(in) :: x(:), dx
    real(kind = 8) :: deriv_fft(size(x))
    complex(kind = 8), allocatable :: x_c(:), k(:)
    integer :: N_orig, N_buff, N_pow2, i, N
    REAL(dp), PARAMETER :: PI   = 2*ASIN(1.0_dp)   ! [-], the constant pi

    N_orig = size(x)
    if (mod(N_orig,2) /= 0) stop 'Number of levels has to be even'
    N_pow2 = 2*2**ceiling(log(dble(N_orig))/log(2D0))
    N_buff = (N_pow2 - N_orig) / 2

    allocate(x_c(2*N_pow2), k(2*N_pow2))
    x_c(N_buff+1 : N_buff+N_orig) = x
    x_c(N_buff:1:-1) = x(1) - [(i, i = 1, N_buff)]*(x(2)-x(1))
    x_c(N_buff+N_orig+1 : N_pow2) = x(N_orig) + [(i, i = 1, N_buff)]*(x(N_orig) - x(N_orig-1))

    x_c(N_pow2+1:size(x_c)) = x_c(N_pow2:1:-1)
    !print *, x_c
    N = size(x_c)

    call fft(x_c)
    !print *, x_c
    k = [ [(i, i = 0, N/2-1)], 0, [(i, i = -N/2+1, -1)] ]
    k = 2*pi*k/(dx*N)
    x_c = cmplx(0D0, 1D0, kind = 8) * cmplx(k, kind = 8) * x_c

    call ifft(x_c)

    deriv_fft = real(x_c(N_buff+1 : N_buff+N_orig))


  end function

end module fft_mod
