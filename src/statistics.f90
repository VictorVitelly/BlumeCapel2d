module statistics
  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use arrays
  use functions
  implicit none

contains


  subroutine hot_start(spin)
    integer(i4), dimension(L,L), intent(out) :: spin
    integer(i4) :: i1,i2
    real(dp) :: r1
    do i1=1,L
      do i2=1,L
        call random_number(r1)
        if(r1 .le. (1._dp/3._dp) ) then
          spin(i1,i2)=1
        elseif(r1 .ge. (2._dp/3._dp)) then
          spin(i1,i2)=-1
        else
          spin(i1,i2)=0
        end if
      end do
    end do
  end subroutine hot_start

  subroutine change_spin(spin,x)
    integer(i32), intent(out) :: spin, x
    spin=SpinVal(x)
  end subroutine change_spin

  subroutine heat_bath(spin,i1,i2,y)
    integer(i32), dimension(L,L), intent(inout) :: spin
    integer(i32), intent(in) :: i1,i2,y
    if(y==3) then
      call heat_bath3(spin,i1,i2)
    else if(y==2) then
      call heat_bath2()
    else
      call heat_bath1(spin,i1,i2)
    end if
  end subroutine heat_bath

  subroutine heat_bath1(spin,i1,i2)
    integer(i32), dimension(L,L), intent(inout) :: spin
    integer(i32), intent(in) :: i1,i2
    real(dp) :: r1
    integer(i32) x,x1,x2,x3
      write(*,*) 'Sitio',i1,i2
      x=5+neighbors(spin,i1,i2)
      x1=maxloc(ptot(:,x),dim=1)
      x2=minloc(ptot(:,x),dim=1)
      x3=6-x1-x2
      call random_number(r1)
      write(*,*) 'Neighbors',neighbors(spin,i1,i2)
      write(*,*) 'Las probabilidades son', ptot(:,x)
      write(*,*) 'r1=',r1
      if(r1>maxval(ptot(:,x))) then
        write(*,*) 'Como r1>',maxval(ptot(:,x)), 'en', x1
        call change_spin(spin(i1,i2),x1)
        write(*,*) 'Cambié el espín a', spin(i1,i2)
      elseif(r1<minval(ptot(:,x))) then
        write(*,*) 'Como r1<',minval(ptot(:,x)), 'en', x2
        call change_spin(spin(i1,i2),x2)
        write(*,*) 'Cambié el espín a', spin(i1,i2)
      else
        write(*,*) 'Quedó en medio', x3
        call change_spin(spin(i1,i2),x3)
        write(*,*) 'Cambié el espín a', spin(i1,i2)
      end if
  end subroutine heat_bath1

  subroutine heat_bath3(spin,i1,i2)
    integer(i32), dimension(L,L), intent(inout) :: spin
    integer(i32), intent(in) :: i1,i2
    real(dp) :: r1
    integer(i32) x,x1,x2,x3
      call random_number(r1)
      if(r1<(1._dp/3._dp) ) then
        call change_spin(spin(i1,i2),1)
      elseif(r1<2._dp/3._dp) then
        call change_spin(spin(i1,i2),2)
      else
        call change_spin(spin(i1,i2),3)
      end if
  end subroutine heat_bath3


end module statistics
