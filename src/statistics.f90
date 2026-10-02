module statistics
  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use arrays
  use functions
  implicit none

contains

! 
  subroutine hot_start(spin)
    integer(i32), dimension(L,L), intent(out) :: spin
    integer(i32) :: i1,i2
    real(dp) :: r1
    do i1=1,L
      do i2=1,L
        call random_number(r1)
        if(r1 .le. (1._dp/3._dp) ) then
          spin(i1,i2)=1
        elseif(r1 .le. (2._dp/3._dp)) then
          spin(i1,i2)=0
        else
          spin(i1,i2)=-1
        end if
      end do
    end do
  end subroutine hot_start
  
  subroutine half_start(spin)
    integer(i32), dimension(L,L), intent(out) :: spin
    integer(i32) :: i1,i2
    spin(:,:)=1
    do i1=1,L
      do i2=1,L/2
        spin(i1,i2)=0
      end do
    end do
  end subroutine half_start
  
  subroutine flip_sign(spin)
    integer(i32), dimension(L,L), intent(inout) :: spin
    spin=-spin
  end subroutine flip_sign
  
  subroutine heat_bath(spin)
    integer(i32), dimension(L,L), intent(inout) :: spin
    integer(i32) :: i1,i2,x
    real(dp) :: r
    do i1=1,L
      do i2=1,L
        x=5+neighbors(spin,i1,i2)
        call random_number(r)
        if(r<p_p1(x)) then 
          spin(i1,i2)=1
        else if(r<p_p1(x)+p_0(x)) then
          spin(i1,i2)=0
        else 
          spin(i1,i2)=-1
        end if
      end do
    end do
  end subroutine heat_bath
  
  subroutine cool_heat_bath(spin,i)
    integer(i32), dimension(L,L), intent(inout) :: spin
    integer(i32), intent(in) :: i
    integer(i32) :: i1,i2,x,n
    real(dp) :: r,r1,r2
    !Lexicographic
    !do i1=1,L
    !  do i2=1,L
    !    x=5+neighbors(spin,i1,i2)
    !    call random_number(r)
    !    if(r<p_p1_T(x,i)) then 
    !      spin(i1,i2)=1
    !    else if(r< (p_p1_T(x,i)+p_0_T(x,i)) ) then
    !      spin(i1,i2)=0
    !    else 
    !      spin(i1,i2)=-1
    !    end if
    !  end do
    !end do
    !Random sites 
    do n=1,L*L
      call random_number(r1); i1 = int(r1*L) +1
      call random_number(r2); i2 = int(r2*L) +1
      x=5+neighbors(spin,i1,i2)
      call random_number(r)
      if(r<p_p1_T(x,i)) then 
        spin(i1,i2)=1
      else if(r< (p_p1_T(x,i)+p_0_T(x,i)) ) then
        spin(i1,i2)=0
      else 
        spin(i1,i2)=-1
      end if  
    end do
  end subroutine cool_heat_bath
  
  subroutine brute_heat_bath(spin,T0,D0)
    integer(i32), dimension(L,L), intent(inout) :: spin
    real(dp), intent(in) :: T0,D0
    integer(i32) :: i1,i2
    real(dp) :: r,p1,p2,p3,AA,x
    do i1=1,L
      do i2=1,L
        x=real(neighbors(spin,i1,i2),dp)
        call random_number(r)
        AA=1._dp+exp(-(-x+D0)/T0)+exp(-(x+D0)/T0)
        p1=exp(-(-x+D0)/T0)/AA
        p2=1._dp/AA
        p3=exp(-(x+D0)/T0)/AA
        if(r<p1) then 
          spin(i1,i2)=1
        else if(r<p1+p2) then
          spin(i1,i2)=0
        else 
          spin(i1,i2)=-1
        end if
      end do
    end do
  end subroutine brute_heat_bath

  subroutine cluster(spin,T)
    real(dp),intent(in) :: T
    integer(i32), dimension(L,L),intent(inout) :: spin
    logical, dimension(L,L) :: bond_x,bond_y
    integer(i32) :: i,j,label(L,L),parent(L*L),next_label,left_label,up_label
    logical, allocatable :: flip_cluster(:)
    real(dp) :: r,p

    do i=1,L
      do j=1,L
        if(spin(i,j)==spin(mod(i,L)+1,j) .and. (spin(i,j).ne.0) ) then
          p=1._dp-exp(-2._dp/T )
          call random_number(r)
          bond_x(i,j)=(r<p)
        else
          bond_x(i,j)=.false.
        end if
        if(spin(i,j)==spin(i,mod(j,L)+1) .and. (spin(i,j).ne.0) ) then
          p=1._dp-exp(-2._dp/T )
          call random_number(r)
          bond_y(i,j)=(r<p)
        else
          bond_y(i,j)=.false.
        end if
      end do
    end do

    label(:,:)=0
    do i=1,L*L
      parent(i)=i
    end do
    next_label=1
    left_label=0
    up_label=0

    do i=1,L
      do j=1,L
        left_label=0
        up_label=0
        if(i>1 .and. bond_x(im(i),j) ) then
          left_label=label(im(i),j)
        end if
        if(j>1 .and. bond_y(i,im(j)) ) then
          up_label=label(i,im(j))
        end if
        if(left_label==0 .and. up_label==0) then
          label(i,j)=next_label
          next_label=next_label+1
        else if(left_label /= 0 .and. up_label==0) then
          label(i,j)=left_label
        else if(left_label== 0 .and. up_label/=0) then
          label(i,j)=up_label
        else
          label(i,j)=min(left_label,up_label)
          call union(left_label,up_label,parent)
        end if
      end do
    end do

    do j=1,L
      if(bond_x(L,j) ) then
        call union(label(1,j),label(L,j),parent )
      end if
    end do

    do i=1,L
      if(bond_y(i,L) ) then
        call union(label(i,1),label(i,L),parent )
      end if
    end do

    do i=1,L
      do j=1,L
        label(i,j)=find(label(i,j),parent)
      end do
    end do

    allocate(flip_cluster(next_label) )
    flip_cluster(:)=.false.

    do i=1,next_label-1
      call random_number(r)
      flip_cluster(i)=(r<0.5_dp)
    end do
    do i=1,L
      do j=1,L
        if(flip_cluster(label(i,j))) then
          spin(i,j)=-spin(i,j)
        end if
      end do
    end do

  end subroutine cluster
  
  subroutine cycles(spin,T,metros)
    real(dp),intent(in) :: T
    integer(i32), intent(in) :: metros
    integer(i32), dimension(L,L),intent(inout) :: spin
    integer(i32) :: i
    do i=1,metros
      call heat_bath(spin)
    end do
    call cluster(spin,T)
  end subroutine cycles
  
  subroutine cluster2(spin,rank0_tot,rank1_tot,rank0n_tot,rank1n_tot) 
    integer(i32), dimension(L,L),intent(in) :: spin
    real(dp), intent(out) :: rank0_tot,rank1_tot,rank0n_tot,rank1n_tot
    logical, dimension(L,L) :: bond_x,bond_y
    integer(i32) :: i,j,label(L,L),parent(L*L),next_label,left_label,up_label
    integer(i32),allocatable :: rank0(:),rank1(:)
    integer(i32) :: rank0_ave,rank0_n,rank1_ave,rank1_n
    
    do i=1,L
      do j=1,L
        if(spin(i,j)==spin(mod(i,L)+1,j) ) then
          bond_x(i,j)=.true.
        else
          bond_x(i,j)=.false.
        end if
        if(spin(i,j)==spin(i,mod(j,L)+1) ) then
          bond_y(i,j)=.true.
        else
          bond_y(i,j)=.false.
        end if
      end do
    end do

    label(:,:)=0
    do i=1,L*L
      parent(i)=i
    end do
    next_label=1
    left_label=0
    up_label=0

    do i=1,L
      do j=1,L
        left_label=0
        up_label=0
        if(i>1 .and. bond_x(im(i),j) ) then
          left_label=label(im(i),j)
        end if
        if(j>1 .and. bond_y(i,im(j)) ) then
          up_label=label(i,im(j))
        end if
        if(left_label==0 .and. up_label==0) then
          label(i,j)=next_label
          next_label=next_label+1  
        else if(left_label /= 0 .and. up_label==0) then
          label(i,j)=left_label
        else if(left_label== 0 .and. up_label/=0) then
          label(i,j)=up_label
        else
          label(i,j)=min(left_label,up_label)
          call union(left_label,up_label,parent)
        end if
      end do
    end do
    
    do j=1,L
      if(bond_x(L,j) ) then
        call union(label(1,j),label(L,j),parent )
      end if
    end do
    
    do i=1,L
      if(bond_y(i,L) ) then
        call union(label(i,1),label(i,L),parent )
      end if 
    end do
    
    do i=1,L
      do j=1,L
        label(i,j)=find(label(i,j),parent)
      end do
    end do
    
    allocate(rank0(next_label) )
    allocate(rank1(next_label) )
    rank0(:)=0 
    rank1(:)=0

    do i=1,L
      do j=1,L
        if(spin(i,j)==0) then
          rank0(label(i,j))=rank0(label(i,j))+1 !cluster size
        else 
          rank1(label(i,j))=rank1(label(i,j))+1 
        end if
      end do
    end do
    
    rank0_ave=0 
    rank0_n=0 
    rank1_ave=0 
    rank1_n=0 
    
    rank0_tot=0._dp
    rank1_tot=0._dp
    do i=1,next_label
      if(rank0(i).ne.0) then 
        rank0_ave=rank0_ave+rank0(i) 
        rank0_n=rank0_n+1 
      end if 
      if(rank1(i).ne.0) then 
        rank1_ave=rank1_ave+rank1(i) 
        rank1_n=rank1_n+1 
      end if 
    end do 
    
    rank0_tot=0._dp
    rank1_tot=0._dp
    if(rank0_ave .ne. 0) then 
      rank0_tot=real(rank0_ave,dp)/real(rank0_n,dp) 
    end if
    if(rank1_ave .ne. 0) then 
      rank1_tot=real(rank1_ave,dp)/real(rank1_n,dp)
    end if
    rank0n_tot=real(rank0_n,dp)
    rank1n_tot=real(rank1_n,dp)
    deallocate(rank0,rank1) !cluster size
    
  end subroutine cluster2


  !ERROR STATISTICS
  subroutine standard_error(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: variance
    integer(i32) :: k,Narr
    Narr=size(x)
    deltay=0._dp
    variance=0._dp
    do k=1,Narr
      variance=variance+(x(k) -y)**2
    end do
    variance=variance/real(Narr-1,dp)
    deltay=Sqrt(variance/real(Narr,dp))
  end subroutine standard_error

  subroutine jackknife(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: jackk
    real(dp), allocatable :: xmean(:), delta_y(:)
    integer(i32) :: k,Narr,i,j
      Narr=size(x)
      allocate(delta_y(size(Mbin)))
      do j=1,size(Mbin)
        allocate(xmean(Mbin(j)))
        jackk=0._dp
        xmean=0._dp
        do i=1,Mbin(j)
          do k=1,Narr
            if(k .le. (i-1)*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            else if(k > i*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            end if
          end do
          xmean(i)=xmean(i)/(real(Narr,dp) -real(Narr/Mbin(j),dp))
        end do
        do k=1,Mbin(j)
          jackk=jackk+(xmean(k)-y )**2
        end do
        delta_y(j)=Sqrt(real(Mbin(j)-1,dp)*jackk/real(Mbin(j),dp))
        deallocate(xmean)
      end do
      deltay=maxval(delta_y)
  end subroutine jackknife

  subroutine mean_0(x,y)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y
    integer(i32) :: k,Narr
    Narr=size(x)
    y=0._dp
    do k=1,Narr
      y=y+x(k)
    end do
    y=y/real(Narr,dp)
  end subroutine mean_0

  subroutine mean_scalar(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y,deltay
    call mean_0(x,y)
    !call standard_error(x,y,deltay)
    call jackknife(x,y,deltay)
  end subroutine mean_scalar
  
  !!Correlation function and length
  
  subroutine initialize(corr1,corr2)
    real(dp), dimension(L), intent(inout) :: corr1
    real(dp), dimension(L,L), intent(inout) :: corr2
      corr1=0._dp
      corr2=0._dp
  end subroutine initialize

  subroutine correlation(spin,corr1,corr2)
    integer(i32), dimension(L,L), intent(in) :: spin
    real(dp), dimension(L), intent(inout) :: corr1
    real(dp), dimension(L,L), intent(inout) :: corr2
    real(dp):: spinvec(L)
    integer(i32) :: i1,i2
    spinvec=0._dp
    !do i1=1,L
    !  do i2=1,L
    !    spinvec(i1)=spinvec(i1)+real(spin(i1,i2),dp)
    !  end do
    !end do
    do i1=1,L
      !corr1(i1)=abs(spinvec(i1))
      do i2=1,L
        corr2(i1,i2)=corr2(i1,i2)+spin(i1,i2)*spin(1,1)
        !corr2(i1,i2)=corr2(i1,i2)+spinvec(i1,i2)*spinvec(1,1)
      end do
    end do
  end subroutine correlation
  
  subroutine correlationoptim(spin,corr1)
    integer(i32), dimension(L,L), intent(in) :: spin
    real(dp), dimension(L), intent(inout) :: corr1
    real(dp):: spinvec(L)
    integer(i32) :: i1,i2
    spinvec=0._dp
    do i1=1,L
      do i2=1,L
        spinvec(i1)=spinvec(i1)+real(spin(i1,i2),dp)
      end do
    end do
    do i1=1,L
        corr1(i1)=corr1(i1)+spinvec(i1)*spinvec(1)
        !corr1(i1)=corr1(i1)+spin(1,1)*spin(i1,1)
    end do
  end subroutine correlationoptim

  subroutine correlation_function(corr1,corr2,CF)
    real(dp), dimension(L), intent(in) :: corr1
    real(dp), dimension(L,L), intent(in) :: corr2
    real(dp), dimension(L), intent(out) :: CF
    integer(i32) :: i1
    do i1=1,L
      CF(i1)=corr2(i1,1)!-corr1(1)**2
    end do
    !CF(:)=CF(:)/real(L**2,dp)
  end subroutine correlation_function
  
  subroutine secondmomentum(CF,xi2_ave,xi2_err)
  real(dp),dimension(L,Nmsrs2),intent(in) :: CF
  real(dp),intent(out) :: xi2_ave,xi2_err 
  integer(i32) :: i1,i2
  real(dp) :: F1(Nmsrs2),F2(Nmsrs2),F12(Nmsrs2),F12_ave,F12_err
  F1(:)=0._dp
  F2(:)=0._dp
  do i1=1,Nmsrs2
    do i2=1,L
      F1(i1)=F1(i1)+CF(i2,i1)
      F2(i1)=F2(i1)+CF(i2,i1)*COS(real(i2-1,dp)*2._dp*PI/real(L,dp))
    end do
  end do
  do i1=1,Nmsrs2
    F12(i1)=F1(i1)/F2(i1)
  end do
  call mean_scalar(F12,F12_ave,F12_err)
  xi2_ave=sqrt( (F12_ave -1._dp))/(2._dp*abs(SIN(PI/real(L,dp))) ) 
  xi2_err=F12_err/(4._dp*sqrt(F12_ave-1._dp)*abs(SIN(PI/real(L,dp))) )
  !write(*,*) xi2_ave,xi2_err
  end subroutine secondmomentum
  
  subroutine secondmomentum2(CF,xi2_ave,xi2_err)
  real(dp),dimension(L,L,Nmsrs2),intent(in) :: CF
  real(dp),intent(out) :: xi2_ave,xi2_err 
  integer(i32) :: i1,i2,i3
  real(dp) :: F1(Nmsrs2),F2(Nmsrs2),F12(Nmsrs2),F12_ave,F12_err
  F1(:)=0._dp
  F2(:)=0._dp
  do i1=1,Nmsrs2
    do i2=1,L
      do i3=1,L
        F1(i1)=F1(i1)+CF(i3,i2,i1)
        F2(i1)=F2(i1)+CF(i3,i2,i1)*COS(real(i2-1,dp)*2._dp*PI/real(L,dp))
      end do
    end do
  end do
  do i1=1,Nmsrs2
    F12(i1)=F1(i1)/F2(i1)
  end do
  call mean_scalar(F12,F12_ave,F12_err)
  xi2_ave=sqrt( (F12_ave -1._dp))/(2._dp*abs(SIN(PI/real(L,dp))) ) 
  xi2_err=F12_err/(4._dp*sqrt(F12_ave-1._dp)*abs(SIN(PI/real(L,dp))) )
  !write(*,*) xi2_ave,xi2_err
  end subroutine secondmomentum2

end module statistics
