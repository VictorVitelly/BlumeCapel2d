program main

  use iso_fortran_env, only : dp => real64, i32 => int32
  use parameters
  use arrays
  use functions
  use statistics
  implicit none

  call cpu_time(starting)
  call init_vecs()
  call hot_start(spin)
  call init_probabilities(1._dp)
  write(*,*) spin
  call heat_bath(spin)
  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

end program main
