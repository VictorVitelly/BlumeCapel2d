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

  !call thermalize(0.1_dp,1.9_dp)
  call vary_T(0.23_dp,0.8_dp,1.95_dp,20)
  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

end program main
