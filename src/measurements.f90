module measurements

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine thermalize(T0,D0)
  real(dp), intent(in) :: T0,D0
  integer(i32) :: i
  open(1, file = 'data/therm.dat', status = 'replace')
    call init_probabilities(T0,D0)
    !call hot_start(spin)
    spin=1
    do i=1,thermalization
      if(i==1 .or. mod(i,2)==0 ) then
        write(1,*) i, Hamilt(spin,D0)/(real(L**2,dp) ), Magnet(spin)/real(L**2,dp)
      end if
      call heat_bath(spin)  
    end do
  close(1)
  end subroutine thermalize
  
  subroutine vary_T(Ti,Tf,D0,Nts)
  real(dp), intent(in) :: Ti,Tf,D0
  integer(i32),intent(in) :: Nts
  integer(i32) :: i1,i2,i3,i4
  real(dp) T0,norm,vol,EE,MM,E_ave,E_err,M_ave,M_err,E2,M2,M4
  real(dp), dimension(Nmsrs2) :: E,M,suscep,heat,U4
  real(dp) :: suscep_ave,suscep_err,heat_ave,heat_err,U4_ave,U4_err
  open(10, file = 'data/energy.dat', status = 'replace')
  open(20, file = 'data/magnetization.dat', status = 'replace')
  open(30, file = 'data/susceptibility.dat', status = 'replace')
  open(40, file = 'data/heat.dat', status = 'replace')
  open(50, file = 'data/binder.dat', status = 'replace')
  norm=real(Nmsrs,dp)
  vol=real(L**2,dp)
  do i1=1,Nts
    T0=Tf+(Ti-Tf)*real(i1-1,dp)/real(Nts-1,dp)
    E(:)=0._dp
    M(:)=0._dp
    write(*,*) i1,T0
    call thermalize(T0,D0)
    do i2=1,Nmsrs2
      E2=0._dp
      M2=0._dp
      M4=0._dp    
      do i3=1,Nmsrs
        do i4=1,eachsweep
          call heat_bath(spin)
        end do
        MM=Magnet(spin)
        EE=Hamilt(spin,D0)
        E(i2)=E(i2)+EE
        M(i2)=M(i2)+abs(MM)
        E2=E2+EE**2
        M2=M2+MM**2
        M4=M4+MM**4
      end do
      E(i2)=E(i2)/norm
      M(i2)=M(i2)/norm
      E2=E2/norm
      M2=M2/norm
      M4=M4/norm
      suscep(i2)=M2-M(i2)**2
      heat(i2)=E2-E(i2)**2
      U4(i2)=1._dp-M4/(3._dp*M2**2)
    end do
    call mean_scalar(E,E_ave,E_err)
    call mean_scalar(M,M_ave,M_err)
    call mean_scalar(suscep,suscep_ave,suscep_err)
    call mean_scalar(heat,heat_ave,heat_err)
    call mean_scalar(U4,U4_ave,U4_err)
    write(10,*) T0, E_ave/vol, E_err/vol
    write(20,*) T0, M_ave/vol, M_err/vol
    write(30,*) T0, suscep_ave/(vol*T0), suscep_err/(vol*T0)
    write(40,*) T0, heat_ave/(vol*T0*T0), heat_err/(vol*T0*T0)
    write(50,*) T0, U4_ave, U4_err
  end do
  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  
  end subroutine vary_T
  
end module measurements
