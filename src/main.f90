program main

  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  call cpu_time(starting)
  call init_vecs()

  !call thermalize(0.66_dp,1.97_dp)
  !call vary_T(0.5_dp,0.8_dp,1.95_dp,10)
  !call correlate(1.69_dp,1.9_dp,0.0_dp,20)
  
  !cool_T(Tc,tauQ,D0)
  call cool_T(1.3_dp,12,1.25_dp)
  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

end program main
