module h2plus
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none

   type(mp_real), save :: one, zero

   private

   public :: init_h2plus

contains
   subroutine int_s11one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knot_xi, knot_eta, Z1, Z2, m, C, R, jz2, epsilon, S11)
      !> @brief This subroutine calculates the S11 integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param knot_xi : real(:) : the knot vector for the xi direction
      !> @param knot_eta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return S11 : real : the value of the S11 integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R, epsilon
      type(mp_real), intent(out) :: S11
      type(mp_real), dimension(:), intent(in) :: knot_xi, knot_eta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta
      integer, intent(in) :: jz2

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta, pol_init, pol_final, diff
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)),prod_eta(size(b_i_xi, 1), size(b_i_xi, 2)))

      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
      end do

      print *, size(prod_eta, 1), size(prod_eta, 2)
      print *, size(b_i_eta, 1), size(b_i_eta, 2)
      print *, size(b_j_eta, 1), size(b_j_eta, 2)
      do i = 1, size(b_i_eta, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      allocate (pol_init(size(b_i_xi, 1), size(b_i_eta, 1)), &
         pol_final(size(b_i_xi, 1), size(b_i_eta, 1)), &
         diff(size(b_i_xi, 1), size(b_i_eta, 1)))

      do i1 = 1, size(b_i_xi, 1) ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) ! Loop over the polynomials of eta
            pol_init(i1, i2) = zero
            pol_final(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 1) - j1
                  beta = size(prod_eta, 1) - j2
                  pol_init(i1, i2) = pol_init(i1, i2) + 2*mppi()*(R**3)*((prod_xi(i1, j1)*knot_xi(i1))**(alpha+1))* ((prod_eta(i2, j2)*knot_eta(i2))**(beta+1)) *((knot_eta(i2)**2)/((alpha+1)*(beta+3))-(knot_xi(i1)**2)/((alpha+3)*(beta+1)))
                  pol_final(i1, i2) = pol_final(i1, i2) + 2*mppi()*(R**3)*((prod_xi(i1, j1)*knot_xi(i1+1))**(alpha+1))* ((prod_eta(i2, j2)*knot_eta(i2+1))**(beta+1)) *((knot_eta(i2+1)**2)/((alpha+1)*(beta+3))-(knot_xi(i1+1)**2)/((alpha+3)*(beta+1)))
               end do
            end do
            diff(i1, i2) = pol_final(i1, i2) - pol_init(i1, i2)
         end do
      end do

      S11 = zero
      do i1 = 1, size(b_i_xi, 1) ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) ! Loop over the polynomials of eta
            S11 = S11 + prod_xi(i1, i2)*prod_eta(i1, i2)*diff(i1, i2)
         end do
      end do

   end subroutine int_s11one

   subroutine init_h2plus(d, n, n_remove, Z1, Z2, m, C, R, ximax, ximin, jz2, epsilon, eta_slp)
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
      type(mp_real), intent(in) :: Z1, Z2, m, C, R, ximax, ximin, epsilon, eta_slp
      integer, intent(in) :: d, n, n_remove, jz2

      type(mp_real), dimension(:), allocatable :: knotxi, knoteta
      type(mp_real), dimension(:, :, :), allocatable :: bspline_xi, bspline_eta
      type(mp_real) :: S11one
      integer :: ntot, i

      ntot = n + d + 2 + 2*n_remove

      ! Generate the knot vectors for xi and eta
      allocate (knotxi(ntot), knoteta(ntot))

      knotxi = knot_xi(d, n, n_remove, ximin, ximax)
      knoteta = knot_eta(d, n, n_remove, eta_slp)

      ! Generate the B-spline coefficients for xi and eta
      allocate (bspline_xi(n, ntot, d), &
         bspline_eta(n, ntot, d))

      do i = 1 + n_remove, n - n_remove
         call init_bspine(d, i, knotxi, bspline_xi(i, :, :), .false.)
         call init_bspine(d, i, knoteta, bspline_eta(i, :, :), .false.)
      end do

      ! Calculate the S11 integral
      call int_s11one(bspline_xi(5, :, :), bspline_eta(8, :, :), &
         bspline_xi(5, :, :), bspline_eta(8, :, :), &
         knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon, S11one)

      ! Print the result
      call mpwrite(6, 30, 10, S11one)

   end subroutine init_h2plus

end module h2plus
