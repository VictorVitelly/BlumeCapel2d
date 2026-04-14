module functions
    use iso_fortran_env, only : dp => real64, i32 => int32
    use parameters
    use arrays
    implicit none

    contains

  function neighbors(spin,i1,i2)
    integer(i32), dimension(L,L), intent(in) :: spin
    integer(i32), intent(in) :: i1,i2
    integer(i32) :: neighbors
    neighbors=spin(ip(i1),i2)+spin(i1,ip(i2))+spin(im(i1),i2)+spin(i1,im(i2))
  end function neighbors

  function Hamilt(spin,D0)
    integer(i32), dimension(L,L), intent(in) :: spin
    real(dp), intent(in) :: D0
    integer(i32) :: i1,i2
    real(dp) :: Hamilt,neigh
    Hamilt=0._dp
    do i1=1,L
      do i2=1,L
        neigh=real(spin(ip(i1),i2)+spin(i1,ip(i2)),dp)
        Hamilt=Hamilt-real(spin(i1,i2),dp)*neigh+D0*real(spin(i1,i2)**2,dp)
      end do
    end do
  end function Hamilt
  
  function Magnet(spin)
    integer(i32), dimension(L,L), intent(in) :: spin
    integer(i32) :: i1,i2
    real(dp) :: Magnet
    Magnet=0._dp
    do i1=1,L
      do i2=1,L
        Magnet=Magnet+spin(i1,i2)
      end do
    end do  
  end function Magnet


end module functions
