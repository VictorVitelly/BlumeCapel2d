module arrays
    use iso_fortran_env, only : dp => real64, i32 => int32
    use parameters, only : L,D
    implicit none

    integer(i32), allocatable :: spin(:,:)
    integer(i32), allocatable :: ip(:), im(:)
    real(dp) :: p_p1(9), p_0(9), p_m1(9), ptot(3,9)
    integer(i32) :: SpinVal(3)=(/1,0,-1/)

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

  subroutine init_probabilities(T)
  real(dp) :: T
  integer(i32) :: i,sigma
  real(dp) :: A
  do i=1,9
    sigma=-5+i
    A=1._dp+2._dp*exp(-D/T)*cosh(real(sigma,dp)/T)
    p_p1(i)=exp(-(-real(sigma,dp)+D)/T)/A
    p_0(i)=1._dp/A
    p_m1(i)=exp(-(real(sigma,dp)+D)/T)/A
    ptot(1,i)=p_p1(i)
    ptot(2,i)=p_0(i)
    ptot(3,i)=p_m1(i)
  end do

  end subroutine init_probabilities


end module arrays
