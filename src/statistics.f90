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
