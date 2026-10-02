module arrays
    use iso_fortran_env, only : dp => real64, i32 => int32
    use parameters, only : L
    implicit none

    integer(i32), allocatable :: spin(:,:)
    integer(i32), allocatable :: ip(:), im(:)
    real(dp) :: p_p1(9), p_0(9), p_m1(9)
    real(dp), allocatable :: p_p1_T(:,:), p_0_T(:,:), p_m1_T(:,:)

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
    
  subroutine set_cool_probabilities(T,D0)
  real(dp), intent(in) :: T(:),D0
  integer(i32) :: i,j,Narr,ntied
  real(dp) :: A,sigma,Ep1,E0,Em1,Emin
  real(dp), parameter :: T_thresh = 6.0e-3_dp, tol = 1.0e-10_dp

  Narr=size(T)
  allocate(p_p1_T(9,Narr),p_0_T(9,Narr), p_m1_T(9,Narr))

  do j=1,Narr
    if (T(j) > T_thresh) then
      do i=1,9
        sigma=real(-5+i,dp)
        A=1._dp+2._dp*exp(-D0/T(j))*cosh(sigma/T(j))
        p_p1_T(i,j)=exp(-(-sigma+D0)/T(j))/A
        p_0_T(i,j)=1._dp/A
        p_m1_T(i,j)=exp(-(sigma+D0)/T(j))/A
      end do
    else
      do i=1,9
        sigma=real(-5+i,dp)
        Ep1 = D0 - sigma
        E0  = 0._dp
        Em1 = D0 + sigma
        Emin = min(Ep1, E0, Em1)

        ntied = 0
        if (abs(Ep1-Emin) < tol) ntied = ntied + 1
        if (abs(E0 -Emin) < tol) ntied = ntied + 1
        if (abs(Em1-Emin) < tol) ntied = ntied + 1

        p_p1_T(i,j) = 0._dp
        p_0_T(i,j)  = 0._dp
        p_m1_T(i,j) = 0._dp
        if (abs(Ep1-Emin) < tol) p_p1_T(i,j) = 1._dp/real(ntied,dp)
        if (abs(E0 -Emin) < tol) p_0_T(i,j)  = 1._dp/real(ntied,dp)
        if (abs(Em1-Emin) < tol) p_m1_T(i,j) = 1._dp/real(ntied,dp)
      end do
    end if
  end do
  
  !do j=1,Narr
  !  do i=1,9
  !    write(*,*) T(j), real(-5+i,dp), p_p1_T(i,j), p_0_T(i,j), p_m1_T(i,j)
  ! end do
  !end do
  end subroutine set_cool_probabilities
  
  subroutine set_cool_probabilities_test(T,D0)
  real(dp), intent(in) :: T(:),D0
  integer(i32) :: i,j,Narr,x(9)
  real(dp) :: A,sigma(9),dHs(3,9)
  Narr=size(T)
  allocate(p_p1_T(9,Narr),p_0_T(9,Narr), p_m1_T(9,Narr))
  do i=1,9
    sigma(i)=real(-5+i,dp)
  end do
  do j=1,Narr-1
    do i=1,9
      A=1._dp+2._dp*exp(-D0/T(j))*cosh(sigma(i)/T(j))
      p_p1_T(i,j)=exp(-(-sigma(i)+D0)/T(j))/A
      p_0_T(i,j)=1._dp/A
      p_m1_T(i,j)=exp(-(sigma(i)+D0)/T(j))/A
    end do
  end do

  p_p1_T(:,Narr)=0._dp
  p_0_T(:,Narr)=0._dp
  p_m1_T(:,Narr)=0._dp 
  do i=1,9
    dHs(1,i)=-sigma(i)+D0
    dHs(2,i)=0
    dHs(3,i)=sigma(i)+D0
  end do
  do i=1,9
    x(i)=minloc(dHs(:,i),dim=1)
  end do
  do i=1,9
    if(x(i)==1 ) then 
      p_p1_T(i,Narr)=1._dp
      !write(*,*) 'dH=',dHs(1,i), 'p(+1)=1'
    else if (x(i)==2 ) then 
      p_0_T(i,Narr)=1._dp
      !write(*,*) 'dH=',dHs(2,i), 'p(0)=1'
    else  
      p_m1_T(i,Narr)=1._dp
      !write(*,*) 'dH=',dHs(3,i), 'p(-1)=1'
    end if
  end do
  
  if(D0 .le. 0.99999_dp) then 
    p_p1_T(4,Narr)=0._dp
    p_0_T(4,Narr)=0._dp
    p_m1_T(4,Narr)=1._dp
  else  if(D0 .ge. 1.00001_dp) then 
    p_p1_T(4,Narr)=0._dp
    p_0_T(4,Narr)=1._dp
    p_m1_T(4,Narr)=0._dp
  else 
    p_p1_T(4,Narr)=0.0_dp
    p_0_T(4,Narr)=0.5_dp
    p_m1_T(4,Narr)=0.5_dp
  end if
  
  if(D0<0.99999_dp) then 
    p_p1_T(6,Narr)=1._dp
    p_0_T(6,Narr)=0._dp
    p_m1_T(6,Narr)=0._dp
  else  if(D0>1.00001_dp) then 
    p_p1_T(6,Narr)=0._dp
    p_0_T(6,Narr)=1._dp
    p_m1_T(6,Narr)=0._dp
  else 
    p_p1_T(6,Narr)=0.5_dp
    p_0_T(6,Narr)=0.5_dp
    p_m1_T(6,Narr)=0._dp
  end if
  
  if(D0<0.00001_dp) then 
    p_p1_T(5,Narr)=1._dp/3._dp
    p_0_T(5,Narr)=1._dp/3._dp
    p_m1_T(5,Narr)=1._dp/3._dp
  end if
  
  !consistency tests
  !do j=Narr-1,Narr
  !  do i=1,9
  !    write(*,*) p_p1_T(i,j), p_0_T(i,j), p_m1_T(i,j)
  !  end do
  !end do
  
  !do i=1,9
  !  select case(x(i))
  !    case(1); p_p1_T(i,Narr) = 1._dp   
  !    case(2); p_0_T(i,Narr)  = 1._dp  
  !    case(3); p_m1_T(i,Narr) = 1._dp   
  !  end select
  !  write(*,*) 'sigma=',sigma(i),'dH=',dHs(x(i),i),' x=',x(i)
  !end do
  
  end subroutine set_cool_probabilities_test



end module arrays
