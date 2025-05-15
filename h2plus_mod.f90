module h2plus
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none

   type(mp_real), save :: one, zero

   private

   public :: init_h2plus

contains
   subroutine int_s11one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, n_remove, knot_xi, knot_eta, Z1, Z2, m, C, R, jz2, epsilon, result)
      !> @brief This subroutine calculates the S11 integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knot_xi : real(:) : the knot vector for the xi direction
      !> @param knot_eta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R, epsilon
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knot_xi, knot_eta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta
      integer, intent(in) :: n_remove, jz2

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta, val_max_max, val_min_min, val_max_min, val_min_max, diff
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)+size(b_j_xi, 2)-1), prod_eta(size(b_i_xi, 1), size(b_i_eta, 2)+size(b_j_eta, 2)-1))

      ! Calculate the product of the B-spline coefficients
      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      allocate (val_max_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_max_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                diff(size(b_i_xi, 1), size(b_i_eta, 1)))

      ! Calculate the value of the S11 integral at each knot and take the difference
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 1) - j1
                  beta = size(prod_eta, 1) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*(R**3)*(knot_xi(i1)**(alpha+1)) * (knot_eta(i2)**(beta+1)) *((knot_eta(i2)**2)/((alpha+1)*(beta+3))-(knot_xi(i1)**2)/((alpha+3)*(beta+1)))
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*(R**3)*((knot_xi(i1+1))**(alpha+1))* ((knot_eta(i2+1))**(beta+1)) *((knot_eta(i2+1)**2)/((alpha+1)*(beta+3))-(knot_xi(i1+1)**2)/((alpha+3)*(beta+1)))
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*(R**3)*((knot_xi(i1))**(alpha+1))* ((knot_eta(i2+1))**(beta+1)) *((knot_eta(i2+1)**2)/((alpha+1)*(beta+3))-(knot_xi(i1)**2)/((alpha+3)*(beta+1)))
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*(R**3)*((knot_xi(i1+1))**(alpha+1))* ((knot_eta(i2))**(beta+1)) *((knot_eta(i2)**2)/((alpha+1)*(beta+3))-(knot_xi(i1+1)**2)/((alpha+3)*(beta+1)))
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

      result = zero
      do i1 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of xi
         do i2 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of eta
            result = result + diff(i1, i2)
         end do
      end do

   end subroutine int_s11one

   subroutine int_s22one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, n_remove, knot_xi, knot_eta, Z1, Z2, m, C, R, jz2, epsilon, result)
      !> @brief This subroutine calculates the S22 integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knot_xi : real(:) : the knot vector for the xi direction
      !> @param knot_eta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R, epsilon
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knot_xi, knot_eta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta
      integer, intent(in) :: n_remove, jz2

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta, val_max_max, val_min_min, val_max_min, val_min_max, diff
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)+size(b_j_xi, 2)-1), prod_eta(size(b_i_xi, 1), size(b_i_eta, 2)+size(b_j_eta, 2)-1))

      ! Calculate the product of the B-spline coefficients
      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      allocate (val_max_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_max_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                diff(size(b_i_xi, 1), size(b_i_eta, 1)))

      ! Calculate the value of the result integral at each knot and take the difference
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 1) - j1
                  beta = size(prod_eta, 1) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R**3*knot_eta(i2)**(1 + beta)*knot_xi(i1)**(1 + alpha)*(((1 + beta)*knot_eta(i2)**2*(5 + beta - (3 + beta)*knot_eta(i2)**2))/(1 + alpha) + ((3 + beta)*(-5 + knot_eta(i2)**4 + beta*(-1 + knot_eta(i2)**4))*knot_xi(i1)**2)/(3 + alpha) - ((5 + beta)*(-3 + knot_eta(i2)**2 + beta*(-1 + knot_eta(i2)**2))*knot_xi(i1)**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R**3*knot_eta(i2+1)**(1 + beta)*knot_xi(i1+1)**(1 + alpha)*(((1 + beta)*knot_eta(i2+1)**2*(5 + beta - (3 + beta)*knot_eta(i2+1)**2))/(1 + alpha) + ((3 + beta)*(-5 + knot_eta(i2+1)**4 + beta*(-1 + knot_eta(i2+1)**4))*knot_xi(i1+1)**2)/(3 + alpha) - ((5 + beta)*(-3 + knot_eta(i2+1)**2 + beta*(-1 + knot_eta(i2+1)**2))*knot_xi(i1+1)**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R**3*knot_eta(i2+1)**(1 + beta)*knot_xi(i1)**(1 + alpha)*(((1 + beta)*knot_eta(i2+1)**2*(5 + beta - (3 + beta)*knot_eta(i2+1)**2))/(1 + alpha) + ((3 + beta)*(-5 + knot_eta(i2+1)**4 + beta*(-1 + knot_eta(i2+1)**4))*knot_xi(i1)**2)/(3 + alpha) - ((5 + beta)*(-3 + knot_eta(i2+1)**2 + beta*(-1 + knot_eta(i2+1)**2))*knot_xi(i1)**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R**3*knot_eta(i2)**(1 + beta)*knot_xi(i1+1)**(1 + alpha)*(((1 + beta)*knot_eta(i2)**2*(5 + beta - (3 + beta)*knot_eta(i2)**2))/(1 + alpha) + ((3 + beta)*(-5 + knot_eta(i2)**4 + beta*(-1 + knot_eta(i2)**4))*knot_xi(i1+1)**2)/(3 + alpha) - ((5 + beta)*(-3 + knot_eta(i2)**2 + beta*(-1 + knot_eta(i2)**2))*knot_xi(i1+1)**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

      result = zero
      do i1 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of xi
         do i2 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of eta
            result = result + diff(i1, i2)
         end do
      end do

   end subroutine int_s22one

   subroutine int_C11one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, n_remove, knot_xi, knot_eta, Z1, Z2, m, C, R, jz2, epsilon, result)
      !> @brief This subroutine calculates the C11one integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
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
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knot_xi, knot_eta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta
      integer, intent(in) :: n_remove, jz2

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta, val_max_max, val_min_min, val_max_min, val_min_max, diff
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)+size(b_j_xi, 2)-1), prod_eta(size(b_i_xi, 1), size(b_i_eta, 2)+size(b_j_eta, 2)-1))

      ! Calculate the product of the B-spline coefficients
      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      allocate (val_max_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_max_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                diff(size(b_i_xi, 1), size(b_i_eta, 1)))

      ! Calculate the value of the S11 integral at each knot and take the difference
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 1) - j1
                  beta = size(prod_eta, 1) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2*mppi()*R*knot_eta(i2)**(1 + beta)*knot_xi(i1)**(1 + alpha)*(((1 + beta)*knot_eta(i2)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knot_eta(i2)))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*knot_xi(i1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knot_xi(i1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2*mppi()*R*knot_eta(i2+1)**(1 + beta)*knot_xi(i1+1)**(1 + alpha)*(((1 + beta)*knot_eta(i2+1)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knot_eta(i2+1)))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*knot_xi(i1+1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knot_xi(i1+1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2*mppi()*R*knot_eta(i2+1)**(1 + beta)*knot_xi(i1)**(1 + alpha)*(((1 + beta)*knot_eta(i2+1)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knot_eta(i2+1)))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*knot_xi(i1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knot_xi(i1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2*mppi()*R*knot_eta(i2)**(1 + beta)*knot_xi(i1+1)**(1 + alpha)*(((1 + beta)*knot_eta(i2)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knot_eta(i2)))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*knot_xi(i1+1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knot_xi(i1+1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

      result = zero
      do i1 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of xi
         do i2 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of eta
            result = result + diff(i1, i2)
         end do
      end do

   end subroutine int_C11one

   subroutine int_C11two(b_i_xi, b_i_eta, b_j_xi, b_j_eta, n_remove, knot_xi, knot_eta, Z1, Z2, m, C, R, jz2, epsilon, result)
      !> @brief This subroutine calculates the C11two integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knot_xi : real(:) : the knot vector for the xi direction
      !> @param knot_eta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R, epsilon
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knot_xi, knot_eta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta
      integer, intent(in) :: n_remove, jz2

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta, val_max_max, val_min_min, val_max_min, val_min_max, diff
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)+size(b_j_xi, 2)-1), prod_eta(size(b_i_xi, 1), size(b_i_eta, 2)+size(b_j_eta, 2)-1))

      ! Calculate the product of the B-spline coefficients
      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      allocate (val_max_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_max_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                diff(size(b_i_xi, 1), size(b_i_eta, 1)))

      ! Calculate the value of the S11 integral at each knot and take the difference
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 1) - j1
                  beta = size(prod_eta, 1) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R*knot_eta(i2)**(1 + beta)*knot_xi(i1)**(1 + alpha)*(((1 + beta)*knot_eta(i2)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knot_eta(i2)))/(1 + alpha) - ((Z1 + Z2)*(2 + beta)*(3 + beta)*knot_xi(i1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knot_xi(i1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R*knot_eta(i2+1)**(1 + beta)*knot_xi(i1+1)**(1 + alpha)*(((1 + beta)*knot_eta(i2+1)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knot_eta(i2+1)))/(1 + alpha) - ((Z1 + Z2)*(2 + beta)*(3 + beta)*knot_xi(i1+1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knot_xi(i1+1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R*knot_eta(i2+1)**(1 + beta)*knot_xi(i1)**(1 + alpha)*(((1 + beta)*knot_eta(i2+1)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knot_eta(i2+1)))/(1 + alpha) - ((Z1 + Z2)*(2 + beta)*(3 + beta)*knot_xi(i1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knot_xi(i1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R*knot_eta(i2)**(1 + beta)*knot_xi(i1+1)**(1 + alpha)*(((1 + beta)*knot_eta(i2)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knot_eta(i2)))/(1 + alpha) - ((Z1 + Z2)*(2 + beta)*(3 + beta)*knot_xi(i1+1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knot_xi(i1+1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

      result = zero
      do i1 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of xi
         do i2 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of eta
            result = result + diff(i1, i2)
         end do
      end do

   end subroutine int_C11two

   subroutine int_C22one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, n_remove, knot_xi, knot_eta, Z1, Z2, m, C, R, jz2, epsilon, result)
      !> @brief This subroutine calculates the C22one integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knot_xi : real(:) : the knot vector for the xi direction
      !> @param knot_eta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R, epsilon
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knot_xi, knot_eta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta
      integer, intent(in) :: n_remove, jz2

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta, val_max_max, val_min_min, val_max_min, val_min_max, diff
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)+size(b_j_xi, 2)-1), prod_eta(size(b_i_xi, 1), size(b_i_eta, 2)+size(b_j_eta, 2)-1))

      ! Calculate the product of the B-spline coefficients
      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      allocate (val_max_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_max_min(size(b_i_xi, 1), size(b_i_eta, 1)), &
                val_min_max(size(b_i_xi, 1), size(b_i_eta, 1)), &
                diff(size(b_i_xi, 1), size(b_i_eta, 1)))

      ! Calculate the value of the S11 integral at each knot and take the difference
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 1) - j1
                  beta = size(prod_eta, 1) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*R*knot_eta(i2)**(1 + beta)*knot_xi(i1)**(1 + alpha)*((knot_eta(i2)*((-Z1 + Z2)/(2 + beta) + (c**2*m*R*knot_eta(i2))/(3 + beta) + ((Z1 - Z2)*knot_eta(i2)**2)/(4 + beta) - (c**2*m*R*knot_eta(i2)**3)/(5 + beta)))/(1 + alpha) - ((Z1 + Z2)*(-3 + knot_eta(i2)**2 + beta*(-1 + knot_eta(i2)**2))*knot_xi(i1))/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*knot_eta(i2)*(1/(2 + beta) - knot_eta(i2)**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + knot_eta(i2)**4/(5 + beta)))*knot_xi(i1)**2)/(3 + alpha) + ((Z1 + Z2)*(-3 + knot_eta(i2)**2 + beta*(-1 + knot_eta(i2)**2))*knot_xi(i1)**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*knot_eta(i2)**2)*knot_xi(i1)**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*R*knot_eta(i2+1)**(1 + beta)*knot_xi(i1+1)**(1 + alpha)*((knot_eta(i2+1)*((-Z1 + Z2)/(2 + beta) + (c**2*m*R*knot_eta(i2+1))/(3 + beta) + ((Z1 - Z2)*knot_eta(i2+1)**2)/(4 + beta) - (c**2*m*R*knot_eta(i2+1)**3)/(5 + beta)))/(1 + alpha) - ((Z1 + Z2)*(-3 + knot_eta(i2+1)**2 + beta*(-1 + knot_eta(i2+1)**2))*knot_xi(i1+1))/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*knot_eta(i2+1)*(1/(2 + beta) - knot_eta(i2+1)**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + knot_eta(i2+1)**4/(5 + beta)))*knot_xi(i1+1)**2)/(3 + alpha) + ((Z1 + Z2)*(-3 + knot_eta(i2+1)**2 + beta*(-1 + knot_eta(i2+1)**2))*knot_xi(i1+1)**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*knot_eta(i2+1)**2)*knot_xi(i1+1)**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*R*knot_eta(i2+1)**(1 + beta)*knot_xi(i1)**(1 + alpha)*((knot_eta(i2+1)*((-Z1 + Z2)/(2 + beta) + (c**2*m*R*knot_eta(i2+1))/(3 + beta) + ((Z1 - Z2)*knot_eta(i2+1)**2)/(4 + beta) - (c**2*m*R*knot_eta(i2+1)**3)/(5 + beta)))/(1 + alpha) - ((Z1 + Z2)*(-3 + knot_eta(i2+1)**2 + beta*(-1 + knot_eta(i2+1)**2))*knot_xi(i1))/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*knot_eta(i2+1)*(1/(2 + beta) - knot_eta(i2+1)**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + knot_eta(i2+1)**4/(5 + beta)))*knot_xi(i1)**2)/(3 + alpha) + ((Z1 + Z2)*(-3 + knot_eta(i2+1)**2 + beta*(-1 + knot_eta(i2+1)**2))*knot_xi(i1)**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*knot_eta(i2+1)**2)*knot_xi(i1)**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*R*knot_eta(i2)**(1 + beta)*knot_xi(i1+1)**(1 + alpha)*((knot_eta(i2)*((-Z1 + Z2)/(2 + beta) + (c**2*m*R*knot_eta(i2))/(3 + beta) + ((Z1 - Z2)*knot_eta(i2)**2)/(4 + beta) - (c**2*m*R*knot_eta(i2)**3)/(5 + beta)))/(1 + alpha) - ((Z1 + Z2)*(-3 + knot_eta(i2)**2 + beta*(-1 + knot_eta(i2)**2))*knot_xi(i1+1))/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*knot_eta(i2)*(1/(2 + beta) - knot_eta(i2)**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + knot_eta(i2)**4/(5 + beta)))*knot_xi(i1+1)**2)/(3 + alpha) + ((Z1 + Z2)*(-3 + knot_eta(i2)**2 + beta*(-1 + knot_eta(i2)**2))*knot_xi(i1+1)**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*knot_eta(i2)**2)*knot_xi(i1+1)**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

      result = zero
      do i1 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of xi
         do i2 = 1 + n_remove, size(b_i_xi, 1) - n_remove ! Loop over the polynomials of eta
            result = result + diff(i1, i2)
         end do
      end do

   end subroutine int_C22one

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
      type(mp_real), dimension(:, :), allocatable :: S11one, S22one, C11one, C11two, C22one
      type(mp_real), dimension(:, :, :), allocatable :: bspline_xi, bspline_eta

      integer :: ntot, i, j
      integer, dimension(2) :: i2, j2

      ntot = n + d + 2*n_remove
      zero = '0.d0'

      print *, "Generating B-spline knots..."
      ! Generate the knot vectors for xi and eta
      allocate (knotxi(ntot), knoteta(ntot))

      knotxi = knot_xi(d, n, n_remove, ximin, ximax)
      knoteta = knot_eta(d, n, n_remove, eta_slp)

      print *, "Generating B-spline coefficients..."
      ! Generate the B-spline coefficients for xi and eta
      allocate (bspline_xi(n, ntot, d), bspline_eta(n, ntot, d))

      do i = 1, n ! Loop over the number of B-splines
         call init_bspine(d, i, knotxi, bspline_xi(i, :, :), .false.)
         call init_bspine(d, i, knoteta, bspline_eta(i, :, :), .false.)
      end do

      print *, "Calculating the S11 integral..."
      ! Calculate the S11 integral
      allocate (S11one(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, S11one, n, d, n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            S11one(i, j) = zero
            call int_s11one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon, S11one(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      print *, "Calculating the S22 integral..."
      ! Calculate the S22 integral
      allocate (S22one(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, S22one, n, d, n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            S22one(i, j) = zero
            call int_s22one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon, S22one(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      print *, "Calculating the C11one integral..."
      ! Calculate the C11one integral
      allocate (C11one(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C11one, n, d, n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C11one(i, j) = zero
            call int_C11one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon, C11one(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      print *, "Calculating the C11two integral..."
      ! Calculate the C11two integral
      allocate (C11two(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C11two, n, d, n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C11two(i, j) = zero
            call int_C11two(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon, C11two(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      print *, "Calculating the C22one integral..."
      ! Calculate the C22one integral
      allocate (C22one(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C22one, n, d, n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C22one(i, j) = zero
            call int_C22one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            n_remove, knotxi, knoteta, Z1, Z2, m, C, R, jz2, epsilon, C22one(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Save matrices to separate files
      open (unit=2, file='S11one.csv', status='replace')
      do i = 1, n**2
         call write_csv(S11one(i, :), 2, 30, 10)
      end do
      close (2)

      open (unit=3, file='S22one.csv', status='replace')
      do i = 1, n**2
         call write_csv(S22one(i, :), 3, 30, 10)
      end do
      close (3)

      open (unit=4, file='C11one.csv', status='replace')
      do i = 1, n**2
         call write_csv(C11one(i, :), 4, 30, 10)
      end do
      close (4)

      open (unit=5, file='C11two.csv', status='replace')
      do i = 1, n**2
         call write_csv(C11two(i, :), 5, 30, 10)
      end do
      close (5)

      open (unit=6, file='C22one.csv', status='replace')
      do i = 1, n**2
         call write_csv(C22one(i, :), 6, 30, 10)
      end do
      close (6)

      ! Save logs
      open (unit=1, file='log_file', status='replace')
      write (1, '(a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d
      write (1, '(a)') "Speed of light: "
      call mpwrite(1, 30, 10, C)
      write (1, '(a)') "Mass of the electron: "
      call mpwrite(1, 30, 10, m)
      write (1, '(a)') "Slope for eta: "
      call mpwrite(1, 30, 10, eta_slp)
      write (1, '(a)') "Xi knot vector: "
      call write_lists(knotxi, 1, 30, 5)
      write (1, '(a)') "Eta knot vector: "
      call write_lists(knoteta, 1, 30, 5)
      write (1, '(a)') "B-Spline basis function Order: "
      do i = 1, n**2
         i2 = indexToPair(i, n)
         write (1, '(i4, a, i4, a, i4)', advance='no') i, " ", i2(1), " ", i2(2)
         write (1, '(a)') " "
      end do
      write (1, '(a)') "--------------------------------------------------------------"
      write (1, '(a)') "S11 integral: "
      do i = 1, n**2
         call write_lists(S11one(i, :), 1, 30, 10)
      end do
      write (1, '(a)') "--------------------------------------------------------------"
      write (1, '(a)') "S22 integral: "
      do i = 1, n**2
         call write_lists(S22one(i, :), 1, 30, 10)
      end do
      write (1, '(a)') "--------------------------------------------------------------"
      write (1, '(a)') "C11one integral: "
      do i = 1, n**2
         call write_lists(C11one(i, :), 1, 30, 10)
      end do
      write (1, '(a)') "--------------------------------------------------------------"
      write (1, '(a)') "C11two integral: "
      do i = 1, n**2
         call write_lists(C11two(i, :), 1, 30, 10)
      end do
      write (1, '(a)') "--------------------------------------------------------------"
      write (1, '(a)') "C22one integral: "
      do i = 1, n**2
         call write_lists(C22one(i, :), 1, 30, 10)
      end do
      write (1, '(a)') "--------------------------------------------------------------"
      close (1)

   end subroutine init_h2plus

end module h2plus
