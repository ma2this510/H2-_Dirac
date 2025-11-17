program main
   use mpmodule
   use bspline_gen
   use h2plus_sep
   implicit none

   integer :: d, n, n_remove, jz2, maxit
   character(len=20) :: Z1read, Z2read, mread, Cread, Rread, ximaxread, ximinread, jz2read, epsilonread, eta_slpread, xi_slpread, eig_read, folder_name
   type(mp_real) :: Z1, Z2, m, C, R, ximax, ximin, epsilon, eta_slp, xi_slp, zero, one, eig
   logical :: save_step, compute_wf, tot_diag

   double precision :: tm_init, tm_fin

   zero = '0.d0'
   one = '1.d0'

   !-------------------------------------------------
   ! Define Important Constants
   jz2 = one ! Quantum number will be divided by 2
   !-------------------------------------------------
   save_step = .false. ! Save matrices at each step
   !-------------------------------------------------

   !-------------------------------------------------
   ! Read Input Parameters
   !-------------------------------------------------

   read(*,*) d

   read(*,*) n

   read(*,*) n_remove

   read(*,*) Z1read
   Z1 = mpreal(Z1read)

   read(*,*) Z2read
   Z2 = mpreal(Z2read)

   read(*,*) mread
   m = mpreal(mread)

   read(*,*) Cread
   C = mpreal(Cread)

   read(*,*) Rread
   R = mpreal(Rread)

   read(*,*) ximaxread
   ximax = mpreal(ximaxread)

   read(*,*) ximinread
   ximin = mpreal(ximinread)

   read(*,*) epsilonread
   epsilon = mpreal(epsilonread)

   read(*,*) eta_slpread
   eta_slp = mpreal(eta_slpread)

   read(*,*) xi_slpread
   xi_slp = mpreal(xi_slpread)

   read(*,*) save_step

   read(*,*) folder_name

   read(*,*) tot_diag

   read(*,*) maxit

   read(*,*) eig_read
   eig = mpreal(eig_read)

   read(*,*) compute_wf

   ! Start the timer
   tm_init = second()
   call init_h2plus_sep(d, n, n_remove, Z1, Z2, m, C, R, ximax, ximin, jz2, epsilon, eta_slp, xi_slp, save_step, folder_name, tot_diag, maxit, eig, compute_wf)
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
