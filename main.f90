program main
   use mpmodule
   use bspline_gen
   use h2plus_sep
   implicit none

   integer :: d, n, n_remove, jz2
   type(mp_real) :: Z1, Z2, m, C, R, ximax, ximin, epsilon, eta_slp, zero, one
   logical :: save_step

   double precision :: tm_init, tm_fin

   zero = '0.d0'
   one = '1.d0'

   !-------------------------------------------------
   ! Define Important Constants
   d = 8 ! Order of the B-Spline (order Mathematica + 1)
   n = 30 ! Number of Usable B-spline
   n_remove = 0 ! Number of knots to remove from each end
   Z1 = '1.0d0' ! number of protons for the first atom
   Z2 = '1.0d0' ! number of protons for the second atom
   m = '1.0d0' ! mass of the electron
   C = '137.035999679d0' ! check CODATA 1986
   R = '1.0d0' ! distance between the two nuclei
   ximax = '30.0d0' ! maximum position of the B-spline on z-axis
   ximin = '1.0d0' ! minimum position of the B-spline on z-axis
   jz2 = one ! Quantum number will be divided by 2
   ! epsilon = '1.0d-5' ! machine epsilon
   epsilon = zero ! machine epsilon
   eta_slp = '4.0d-2' ! paramet for the generation of the knot vector on eta
   !-------------------------------------------------
   save_step = .false. ! Save matrices at each step
   !-------------------------------------------------

   ! Start the timer
   tm_init = second()
   call init_h2plus_sep(d, n, n_remove, Z1, Z2, m, C, R, ximax, ximin, jz2, epsilon, eta_slp, save_step)
   tm_fin = second()
   print *, "Elapsed time for global execution: ", tm_fin - tm_init, " seconds"
end program main

function epsilonn(alpha)
   !> @brief Calculate the machine epsilon
   !> @param alpha The value to calculate the machine epsilon
   USE mpmodule
   implicit type(mp_real) (a - h, o - z)

   ten = '10.d0'
   epsilonn = ten**(-mpipl)

   return
end function epsilonn
