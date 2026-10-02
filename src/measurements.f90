module measurements

  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine thermalize(T0,D0)
  real(dp), intent(in) :: T0,D0
  integer(i32) :: i
  !open(1, file = 'data/therm.dat', status = 'replace')
  !open(2, file = 'data/configs.dat', status = 'replace')
    call init_probabilities(T0,D0)
    !call hot_start(spin)
    !spin=0
    !call half_start(spin)
    do i=1,thermalization
      !if(i==1 .or. mod(i,2)==0 ) then
        !write(1,*) i, Hamilt(spin,D0)/(real(L**2,dp) ), Magnet(spin)/real(L**2,dp
      !end if
      !call cycles(spin,T0,500)
      call heat_bath(spin)
    end do
  !close(1)
  !close(2)
  end subroutine thermalize
  
  subroutine vary_T(Ti,Tf,D0,Nts)
  real(dp), intent(in) :: Ti,Tf,D0
  integer(i32),intent(in) :: Nts
  integer(i32) :: i1,i2,i3,i4
  real(dp) T0,norm,vol,EE,MM,E_ave,E_err,M_ave,M_err,E2,M2,M4
  real(dp), dimension(Nmsrs2) :: E,M,suscep,heat,U4,rank0,rank1
  real(dp) :: suscep_ave,suscep_err,heat_ave,heat_err,U4_ave,U4_err
  real(dp) :: rank0_tot,rank1_tot,rank0_ave,rank0_err,rank1_ave,rank1_err
  real(dp) :: x1,x2
  open(10, file = 'data/energy.dat', status = 'replace')
  open(20, file = 'data/magnetization.dat', status = 'replace')
  open(30, file = 'data/susceptibility.dat', status = 'replace')
  open(40, file = 'data/heat.dat', status = 'replace')
  open(50, file = 'data/binder.dat', status = 'replace')
  open(60, file = 'data/isles.dat', status = 'replace')
  
  norm=real(Nmsrs,dp)
  vol=real(L**2,dp)
  do i1=1,Nts
    T0=Tf+(Ti-Tf)*real(i1-1,dp)/real(Nts-1,dp)
    !call init_probabilities(T0,D0)
    E(:)=0._dp
    M(:)=0._dp
    rank0(:)=0._dp
    rank1(:)=0._dp
    write(*,*) i1,T0
    call thermalize(T0,D0)
    do i2=1,Nmsrs2
      E2=0._dp
      M2=0._dp
      M4=0._dp    
      do i3=1,Nmsrs
        do i4=1,eachsweep
          call cycles(spin,T0,4)
        end do
        MM=Magnet(spin)
        EE=Hamilt(spin,D0)
        E(i2)=E(i2)+EE
        M(i2)=M(i2)+abs(MM)
        E2=E2+EE**2
        M2=M2+MM**2
        M4=M4+MM**4
        call cluster2(spin,rank0_tot,rank1_tot,x1,x2) 
        rank0(i2)=rank0(i2)+rank0_tot
        rank1(i2)=rank1(i2)+rank1_tot
      end do
      E(i2)=E(i2)/norm
      M(i2)=M(i2)/norm
      E2=E2/norm
      M2=M2/norm
      M4=M4/norm
      suscep(i2)=M2-M(i2)**2
      heat(i2)=E2-E(i2)**2
      U4(i2)=1._dp-M4/(3._dp*M2**2)
      rank0(i2)=rank0(i2)/(norm*vol)
      rank1(i2)=rank1(i2)/(norm*vol)
    end do
    call mean_scalar(E,E_ave,E_err)
    call mean_scalar(M,M_ave,M_err)
    call mean_scalar(suscep,suscep_ave,suscep_err)
    call mean_scalar(heat,heat_ave,heat_err)
    call mean_scalar(U4,U4_ave,U4_err)
    call mean_scalar(rank0,rank0_ave,rank0_err)
    call mean_scalar(rank1,rank1_ave,rank1_err)
    write(10,*) T0, E_ave/vol, E_err/vol
    write(20,*) T0, M_ave/vol, M_err/vol
    write(30,*) T0, suscep_ave/(vol*T0), suscep_err/(vol*T0)
    write(40,*) T0, heat_ave/(vol*T0*T0), heat_err/(vol*T0*T0)
    write(50,*) T0, U4_ave, U4_err
    write(60,*) T0, rank0_ave, rank0_err, rank1_ave, rank1_err
    !do i4=1,L
    !  write(*,*) spin(i4,:)
    !end do  
    !write(*,*) '_________________________________________________'
  end do
  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  close(60)
  end subroutine vary_T
  
  subroutine correlate(Ti,Tf,D0,NTs)
  real(dp), intent(in) :: Ti,Tf,D0
  integer(i32), intent(in) :: NTs
  integer(i32) :: i,j,k,i2
  real(dp), allocatable :: corr1(:)
  real(dp), allocatable :: corr2(:,:)
  real(dp), allocatable, dimension(:,:) :: CF,CF_ave,CF_err
  !real(dp) :: M(Nmsrs2),MM(Nmsrs2),M_ave,M_err
  real(dp) :: T0,vol,norm,xi2_ave,xi2_err
  open(70, file = 'data/corrfunc.dat', status = 'replace')
  open(80, file = 'data/xi2.dat', status = 'replace')
    vol=real(L**2,dp)
    norm=real(Nmsrs,dp)
    allocate(corr1(L))
    allocate(corr2(L,L))
    allocate(CF(L,Nmsrs2))
    allocate(CF_ave(L,Nts))
    allocate(CF_err(L,Nts))
    do k=1,Nts
      T0=Ti+(Tf-Ti)*real(k-1,dp)/real(Nts-1)
      call init_probabilities(T0,D0)
      write(*,*) k, 'de', Nts 
      !M(:)=0._dp
      call hot_start(spin)
      do j=1,thermalization
        call cycles(spin,T0,4)
      end do
      do j=1,Nmsrs2
        !call hot_start(spin)
        !do i=1,thermalization
        !  call cycles(spin,T0,4)
        !end do
        call initialize(corr1,corr2)
        corr1(:)=0._dp
        do i=1,Nmsrs
          do i2=1,eachsweep       
            call cycles(spin,T0,4)
          end do
          !M(j)=M(j)+abs(magnet(spin))/vol
          call correlationoptim(spin,corr1)
        end do
        !MM(j)=(M(j)**2)/(norm**2)
        corr1(:)=corr1(:)/norm
        CF(:,j)=corr1(:)
      end do
      !call mean_scalar(MM,M_ave,M_err)
      do j=1,L
        call mean_scalar(CF(j,:),CF_ave(j,k) ,CF_err(j,k))
      end do
      call secondmomentum(CF,xi2_ave,xi2_err)
      write(80,*) T0,xi2_ave,xi2_err
      do j=1,L
        call mean_scalar(CF(j,:),CF_ave(j,k) ,CF_err(j,k))
      end do
    end do
    
    do k=1,L+1
      write(70,*) abs(k-1), CF_ave(iv(k),:), CF_err(iv(k),:)
    end do
    close(70)
    close(80)
    deallocate(spin)
    deallocate(corr1,corr2,CF,CF_ave,CF_err)
  end subroutine correlate
  
  subroutine cool_T(Tc,tau,D0)
  real(dp), intent(in) :: Tc,D0
  integer(i32),intent(in) :: tau
  integer(i32) :: i1,i2,npoints,measurements,n_zeros
  real(dp) norm,vol,EE,MM,E_ave,E_err,M_ave,M_err !,E2,M2,M4
  real(dp) :: rank0_tot,rank1_tot,rank0_ave,rank0_err,rank1_ave,rank1_err
  real(dp), allocatable :: E(:,:),M(:,:),T(:),rank0(:,:),rank1(:,:)
  real(dp), allocatable :: rank0n(:,:),rank1n(:,:),n_0(:,:)
  real(dp) :: rank0n_tot,rank1n_tot,rank0n_ave,rank0n_err,rank1n_ave,rank1n_err
  real(dp) :: n_zeros_ave,n_zeros_err
  integer(i32) ::spin(L,L)
  open(10, file = 'data/ene.dat', status = 'replace')
  open(20, file = 'data/mag.dat', status = 'replace')
  open(30, file = 'data/rnk.dat', status = 'replace')
  open(40, file = 'data/ncl.dat', status = 'replace')
  open(50, file = 'data/n0s.dat', status = 'replace')
  norm=real(Nmsrs,dp)
  vol=real(L**2,dp)
  npoints=2*tau+1
  !npoints=2*tau
  measurements=Nmsrs2*Nmsrs 
  
  allocate(T(npoints),E(measurements,npoints),M(measurements,npoints))
  allocate(rank0(measurements,npoints),rank1(measurements,npoints))  
  allocate(rank0n(measurements,npoints),rank1n(measurements,npoints))
  allocate(n_0(measurements,npoints))
  do i1=1,npoints
    T(i1)=Tc*(1._dp+real(tau+1-i1,dp)/real(tau,dp) )
    !T(i1)=2._dp*Tc*(1._dp-real(i1-1,dp)/real(npoints-1,dp) )
    !write(*,*) T(i1)
  end do
  call init_probabilities(T(1),D0)
  call set_cool_probabilities(T,D0)
  E=0._dp
  M=0._dp
  rank0=0._dp
  rank1=0._dp
  rank0n=0._dp
  rank1n=0._dp
  n_0=0._dp
  
  !spin=0
  !call thermalize(T(1),D0)
  do i1=1,measurements
    !write(*,*) i1
    !spin=0
    call hot_start(spin)
    call thermalize(T(1),D0)
    do i2=1,npoints
      call cool_heat_bath(spin,i2)
      !call brute_heat_bath(spin,T(i2),D0)
      MM=Magnet(spin)
      EE=Hamilt(spin,D0)
      E(i1,i2)=EE
      M(i1,i2)=abs(MM)
      call cluster2(spin,rank0_tot,rank1_tot,rank0n_tot,rank1n_tot) 
      rank0(i1,i2)=rank0_tot
      rank1(i1,i2)=rank1_tot
      rank0n(i1,i2)=rank0n_tot
      rank1n(i1,i2)=rank1n_tot
      call count_zeros(spin, n_zeros)
      n_0(i1,i2)=real(n_zeros,dp)
    end do
  end do
  
  do i1=1,npoints
    call mean_scalar(E(:,i1),E_ave,E_err)
    call mean_scalar(M(:,i1),M_ave,M_err)
    call mean_scalar(rank0(:,i1),rank0_ave,rank0_err)
    call mean_scalar(rank1(:,i1),rank1_ave,rank1_err)
    call mean_scalar(rank0n(:,i1),rank0n_ave,rank0n_err)
    call mean_scalar(rank1n(:,i1),rank1n_ave,rank1n_err)    
    call mean_scalar(n_0(:,i1),n_zeros_ave,n_zeros_err)
    write(10,*) i1,T(i1),E_ave/vol,E_err/vol
    write(20,*) i1,T(i1),M_ave/vol,M_err/vol
    write(30,*) i1,T(i1),rank1_ave/vol,rank1_err/vol,rank0_ave/vol,rank0_err/vol 
    write(40,*) i1,T(i1),rank1n_ave/vol,rank1n_err/vol,rank0n_ave/vol,rank0n_err/vol 
    write(50,*) i1,T(i1),n_zeros_ave/vol,n_zeros_err/vol
  end do

  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  deallocate(T,E,M,rank0,rank1,rank0n,rank1n)

  end subroutine cool_T
  
  
end module measurements
