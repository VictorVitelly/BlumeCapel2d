module arrays
    use iso_fortran_env, only : dp => real64, i32 => int32
    use parameters, only : L,D
    implicit none

    integer(i32), allocatable :: spin(:,:)
    integer(i32), allocatable :: ip(:), im(:)
    real(dp) :: p_p1(9), p_0(9), p_m1(9)

contains

  subroutine init_vecs()
    integer(i32) :: i
    allocate(spin(L,L))
    allocate(ip(L),im(L))
    do i=1,L-1
      ip(i)=i+1
    end do
    ip(L)=1
    do i=2,L
      im(i)=i-1
    end do
    im(1)=L
  end subroutine init_vecs

  subroutine init_probabilities(T0,D0)
  real(dp), intent(in) :: T0,D0
  integer(i32) :: i
  real(dp) :: A,sigma
  do i=1,9
    sigma=real(-5+i,dp)
    A=1._dp+2._dp*exp(-D0/T0)*cosh(sigma/T0)
    p_p1(i)=exp(-(-sigma+D0)/T0)/A
    p_0(i)=1._dp/A
    p_m1(i)=exp(-(sigma+D0)/T0)/A
  end do

  end subroutine init_probabilities


end module arrays
