module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: L=8,thermalization=5000,eachsweep=150,Nmsrs=120
    integer(i4),parameter :: Nmsrs2=120,Mbin(5)=(/4,5,10,15,20/)
    real :: starting,ending
    real(dp) :: D,T
    real(dp), parameter :: PI=4._dp*Atan(1.0_dp)

end module parameters
