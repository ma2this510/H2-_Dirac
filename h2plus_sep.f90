module h2plus_sep
    use mpmodule
    use bspline_gen
    use tools_mp
    implicit none
    
    type(mp_real), save :: one, zero

    type(mp_real), dimension(:, :), allocatable, save :: pwr_knot_xi, pwr_knot_eta

   private

   public :: init_h2plus_sep

contains
    
    function get_pw_xi(i, k) result(result)
      !> @brief This function returns the power of the xi knot vector.
      !> @param i : integer : the index of the knot vector
      !> @param k : integer : the power of the knot vector
      !> @return result : real : the value of the power of the xi knot vector
      integer, intent(in) :: i, k
      type(mp_real) :: result

      if (allocated(pwr_knot_xi)) then
         if (k > size(pwr_knot_xi, 2)) then
            print *, 'Warning: k exceeds the size of pwr_knot_xi.'
            print *, 'k = ', k
            result = zero
         else
            result = pwr_knot_xi(i, k)
         end if
      else
         print *, 'Error: pwr_knot_xi is not allocated.'
         result = zero
      end if
   end function get_pw_xi

   function get_pw_eta(i, k) result(result)
      !> @brief This function returns the power of the eta knot vector.
      !> @param i : integer : the index of the knot vector
      !> @param k : integer : the power of the knot vector
      !> @return result : real : the value of the power of the eta knot vector
      integer, intent(in) :: i, k
      type(mp_real) :: result

      if (allocated(pwr_knot_eta)) then
         if (k > size(pwr_knot_eta, 2)) then
            print *, 'Warning: k exceeds the size of pwr_knot_eta.'
            print *, 'k = ', k
            result = zero
         else
            result = pwr_knot_eta(i, k)
         end if
      else
         print *, 'Error: pwr_knot_eta is not allocated.'
         result = zero
      end if
   end function get_pw_eta

   subroutine int_C11one(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C11one integral for the H2+ molecule.
      !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
      !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @return S11 : real : the value of the S11 integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out), dimension(:, :) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

      type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
      type(mp_real), dimension(:, :, :, :), allocatable :: prod_xi, prod_eta
      integer :: alpha, beta, i, j, k, l
      integer, dimension(2) :: i2, j2

      allocate (prod_xi(size(b_xi, 1), size(b_xi, 1), size(b_xi, 2), size(b_xi, 3)+size(b_xi, 3)-1), prod_eta(size(b_eta, 1), size(b_eta, 1), size(b_eta, 2), size(b_eta, 3)+size(b_eta, 3)-1))

      do i = 1, size(b_xi, 1) ! Loop over the first basis functions
         do j = 1, size(b_xi, 1) ! Loop over the second basis functions
            do k = 1, size(b_xi, 2) ! Loop over the polynomials
               ! Calculate the product of the B-spline coefficients
               prod_xi(i, j, k, :) = fusion_coef(b_xi(i, k, :), b_xi(j, k, :))
               prod_eta(i, j, k, :) = fusion_coef(b_eta(i, k, :), b_eta(j, k, :))
            end do
         end do
      end do

      allocate(xi_1(size(prod_xi, 1), size(prod_xi, 2)), xi_2(size(prod_xi, 1), size(prod_xi, 2)), eta_1(size(prod_eta, 1), size(prod_eta, 2)), eta_2(size(prod_eta, 1), size(prod_eta, 2)))

      xi_1 = zero
      xi_2 = zero
      eta_1 = zero
      eta_2 = zero

      do i = 1, size(b_xi, 1) ! Loop over first basis functions
         do j = 1, size(b_xi, 1) ! Loop over second basis functions
            do k = 1, size(b_xi, 2) - 1! Loop over the polynomials
                  do l = 1, size(b_xi, 3) ! loop over order
                     alpha = size(b_xi, 3) - l
                     xi_1(i, j) = xi_1(i, j) + prod_xi(i, j, k, l)*knotxi(k+1)**(2 + alpha)*(-((Z1 + Z2)/(2*R + R*alpha)) + (c**2*m*knotxi(k+1))/(3 + alpha)) ! Finale
                     xi_1(i, j) = xi_1(i, j) - prod_xi(i, j, k, l)*knotxi(k)**(2 + alpha)*(-((Z1 + Z2)/(2*R + R*alpha)) + (c**2*m*knotxi(k))/(3 + alpha)) ! Initial

                     xi_2(i, j) = xi_2(i, j) +  prod_xi(i, j, k, l)*knotxi(k+1)**(1 + alpha)/(1 + alpha) ! Finale
                     xi_2(i, j) = xi_2(i, j) -  prod_xi(i, j, k, l)*knotxi(k)**(1 + alpha)/(1 + alpha) ! Initial
                  end do
            end do

            do k = 1, size(b_eta, 2) - 1 ! Loop over the polynomials
                  do l = 1, size(b_eta, 3) ! loop over order
                     beta = size(b_eta, 3) - l
                     eta_1(i, j) = eta_1(i, j) + prod_eta(i, j, k, l)*knoteta(k + 1)**(1 + beta)/(1 + beta) ! Finale
                     eta_1(i, j) = eta_1(i, j) - prod_eta(i, j, k, l)*knoteta(k)**(1 + beta)/(1 + beta) ! Initial

                     eta_2(i, j) = eta_2(i, j) - prod_eta(i, j, k, l)*((knoteta(k + 1)**(2 + beta)*((-Z1 + Z2)/(2 + beta) + (c**2*m*R*knoteta(k + 1))/(3 + beta)))/R) ! Finale
                     eta_2(i, j) = eta_2(i, j) + prod_eta(i, j, k, l)*((knoteta(k)**(2 + beta)*((-Z1 + Z2)/(2 + beta) + (c**2*m*R*knoteta(k))/(3 + beta)))/R) ! Initial
                  end do
            end do
         end do
      end do

      open(10, file="test.csv", status='replace')
      do i = 1, size(xi_1, 1) 
         call write_csv(xi_1(i, :), 10, 30, 10)
      end do
      close(10)

      result = zero

      do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
         do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
            i2 = indexToPair(i, size(b_xi, 1))
            j2 = indexToPair(j, size(b_eta, 1))

            result(i, j) = 2 * mppi() * (R**3) * ( xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) + xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)) )
         end do
      end do
            

    end subroutine int_C11one



    subroutine init_h2plus_sep(d, n, n_remove, Z1, Z2, m, C, R, ximax, ximin, jz2, epsilon, eta_slp, save_step)
      !> @brief This subroutine initializes the B-spline coefficients for the H2+ molecule.
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of usable B-splines
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param ximax : real : the maximum position of the B-spline on xi-axis
      !> @param ximin : real : the minimum position of the B-spline on xi-axis
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param epsilon : real : the machine epsilon
      !> @param eta_slp : real : the parameter for the generation of the knot vector on eta
      !> @param save_step : boolean : whether to save matrices to a file
      type(mp_real), intent(in) :: Z1, Z2, m, C, R, ximax, ximin, epsilon, eta_slp
      integer, intent(in) :: d, n, n_remove, jz2
      logical, intent(in) :: save_step

      double precision :: tm0, tm1

      type(mp_real), dimension(:), allocatable :: knotxi, knoteta, knotxi_eps, knoteta_eps, knotxi_tmp
      type(mp_real), dimension(:, :), allocatable :: S11one, S22one, C11one, C11two, C22one, C22two, C11three, C22three,C12three, C21three, C_mat, S_mat, vect
      type(mp_real), dimension(:), allocatable :: w, fv1, fv2
      type(mp_real), dimension(:, :, :), allocatable :: bspline_xi, bspline_eta
      logical :: debug_bool = .false.

      integer :: ntot, i, j, ierr

      ntot = n + d + 2*n_remove
      zero = '0.d0'
      one = '1.d0'

      write (6, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d, " and Number of BSplines to remove: ", n_remove

      print *, "Generating B-spline knots..."
      tm0 = second()
      ! Generate the knot vectors for xi and eta
      allocate (knotxi_tmp(ntot+1), knotxi(ntot), knoteta(ntot))

      knotxi_tmp = knot_xi(d, n+1, n_remove, ximin, ximax)
      knotxi = knotxi_tmp(1:ntot) ! Remove the last knot to avoid singularities
      knoteta = knot_eta(d, n, n_remove, eta_slp)

      tm1 = second()
      print *, "Time taken to generate knots: ", tm1 - tm0, " seconds"

      print *, "Pre-Computing the powers of the knots..."
      tm0 = second()
      ! Pre-compute the powers of the knots for xi and eta
      allocate (pwr_knot_xi(ntot, 3*d), pwr_knot_eta(ntot, 3*d))

      !$OMP PARALLEL DO default(shared) private(i, j)
      do i = 1, ntot ! Loop over the number of knots
         do j = 1, 3*d ! Loop over the powers WARNING MEMORY LEAK ?
            if (j == 1) then
               pwr_knot_xi(i, j) = knotxi(i)
               pwr_knot_eta(i, j) = knoteta(i)
            else
               pwr_knot_xi(i, j) = pwr_knot_xi(i, j-1) * knotxi(i)
               pwr_knot_eta(i, j) = pwr_knot_eta(i, j-1) * knoteta(i)
            end if
         end do
      end do
      !$OMP END PARALLEL DO

      tm1 = second()
      print *, "Time taken to pre-compute powers of knots: ", tm1 - tm0, " seconds"

      print *, "Generating B-spline coefficients..."
      tm0 = second()
      ! Generate the B-spline coefficients for xi and eta
      allocate (bspline_xi(n, ntot, d), bspline_eta(n, ntot, d))

      do i = 1 + n_remove, n + n_remove ! Loop over the number of B-splines
         if (debug_bool) then
            print *, "Generating B-spline coefficients for B-spline ", i, "on xi-axis..."
         end if
         call init_bspine(d, i, knotxi_tmp, ntot, bspline_xi(i - n_remove, :, :), debug_bool)
         if (debug_bool) then
            print *, "Generating B-spline coefficients for B-spline ", i, "on eta-axis..."
         end if
         call init_bspine(d, i, knoteta, ntot, bspline_eta(i - n_remove, :, :), debug_bool)
      end do

      ! Modify the knots to avoid singularities.
      allocate (knotxi_eps(ntot), knoteta_eps(ntot))
      knotxi_eps = knotxi
      knoteta_eps = knoteta
      do i = 1, d
         knotxi_eps(i) = knotxi_eps(i) + epsilon
         knoteta_eps(i) = knoteta_eps(i) + epsilon
         knoteta_eps(ntot - i + 1) = knoteta_eps(ntot - i + 1) - epsilon
      end do

      if (debug_bool) then
         print *, "Adjusted B-spline knots for xi-axis: "
         call write_lists(knotxi_eps, 6, 25, 10)
         print *, "Adjusted B-spline knots for eta-axis: "
         call write_lists(knoteta_eps, 6, 25, 10)
      end if

      tm1 = second()
      print *, "Time taken to generate B-spline coefficients: ", tm1 - tm0, " seconds"

      print *, "Calculating the C11one integral..."
      tm0 = second()
      ! Calculate the C11one integral
      allocate (C11one(n**2, n**2))

      call int_C11one(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C11one)

      tm1 = second()
      print *, "Time taken to calculate C11one integral: ", tm1 - tm0

      open(10, file="C11one_sep.csv", status='replace')
      do i = 1, n**2
         call write_csv(C11one(i, :), 10, 30, 10)
      end do
      close(10)


   end subroutine init_h2plus_sep
      

end module h2plus_sep