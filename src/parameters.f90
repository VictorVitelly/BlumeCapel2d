module parameters

    use iso_fortran_env, only : dp => real64, i32 => int32
    implicit none

    integer(i32), parameter :: L=64,thermalization=10000,eachsweep=25,Nmsrs=120
    integer(i32),parameter :: Nmsrs2=20,Mbin(5)=(/4,5,10,15,20/)
    real :: starting,ending
    real(dp), parameter :: PI=4._dp*Atan(1.0_dp)

end module parameters
