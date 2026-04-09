module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: L=4,thermalization=5000,eachsweep=20,Nmsrs=120
    integer(i4),parameter :: Nmsrs2=120,Mbin(5)=(/4,5,10,15,20/)
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    real :: starting,ending
    real(dp) :: D=1._dp
    real(dp), parameter :: PI=4._dp*Atan(1.0_dp)

end module parameters
