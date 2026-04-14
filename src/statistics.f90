module statistics
  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use arrays
  use functions
  implicit none

contains


  subroutine hot_start(spin)
    integer(i32), dimension(L,L), intent(out) :: spin
    integer(i32) :: i1,i2
    real(dp) :: r1
    do i1=1,L
      do i2=1,L
        call random_number(r1)
        if(r1 .le. (1._dp/3._dp) ) then
          spin(i1,i2)=1
        elseif(r1 .ge. (2._dp/3._dp)) then
          spin(i1,i2)=0
        else
          spin(i1,i2)=-1
        end if
      end do
    end do
  end subroutine hot_start
  
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

  subroutine cluster(spin,T,D)
    real(dp),intent(in) :: T,D
    integer(i4), dimension(L,L),intent(inout) :: spin
    logical, dimension(L,L) :: bond_x,bond_y
    integer(i4) :: i,j,label(L,L),parent(L*L),next_label,left_label,up_label
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
        if(i>1 .and. bond_x(i-1,j) ) then
          left_label=label(i-1,j)
        end if
        if(j>1 .and. bond_y(i,j-1) ) then
          up_label=label(i,j-1)
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

end module statistics
