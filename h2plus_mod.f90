module h2plus
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none

   type(mp_real), save :: one, zero

   private

   public :: init_h2plus

contains
   function fun_s11one(xi, eta, Z1, Z2, m, C, R, alpha, beta) result(s11one)
      !> @brief This subroutine calculates the S11one function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: jz2, alpha, beta
      type(mp_real) :: s11one

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      s11one = 2*mppi()*(R**3)*(xi**(alpha+1)) * (eta**(beta+1)) *(-(eta**2)/((alpha+1)*(beta+3))+(xi**2)/((alpha+3)*(beta+1)))
   end function fun_s11one

   subroutine int_s11one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the S11 integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

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

      allocate (val_max_max(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_min_min(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_max_min(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_min_max(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                diff(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1))

      result = zero
      ! Calculate the value of the S11 integral at each knot and take the difference
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            diff(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 2) - j1
                  beta = size(prod_eta, 2) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*(R**3)*(knotxi(i1)**(alpha+1)) * (knoteta(i2)**(beta+1)) *(-(knoteta(i2)**2)/((alpha+1)*(beta+3))+(knotxi(i1)**2)/((alpha+3)*(beta+1)))
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*(R**3)*((knotxi(i1+1))**(alpha+1))* ((knoteta(i2+1))**(beta+1)) *(-(knoteta(i2+1)**2)/((alpha+1)*(beta+3))+(knotxi(i1+1)**2)/((alpha+3)*(beta+1)))
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*(R**3)*((knotxi(i1))**(alpha+1))* ((knoteta(i2+1))**(beta+1)) *(-(knoteta(i2+1)**2)/((alpha+1)*(beta+3))+(knotxi(i1)**2)/((alpha+3)*(beta+1)))
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*2*mppi()*(R**3)*((knotxi(i1+1))**(alpha+1))* ((knoteta(i2))**(beta+1)) *(-(knoteta(i2)**2)/((alpha+1)*(beta+3))+(knotxi(i1+1)**2)/((alpha+3)*(beta+1)))
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
            result = result + val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

   end subroutine int_s11one

   function fun_s22one(xi, eta, Z1, Z2, m, C, R, alpha, beta) result(s22one)
      !> @brief This subroutine calculates the S22one function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: jz2, alpha, beta
      type(mp_real) :: s22one

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      s22one = (2*mppi()*R**3*eta**(1 + beta)*xi**(1 + alpha)*(((1 + beta)*eta**2*(5 + beta - (3 + beta)*eta**2))/(1 + alpha) + ((3 + beta)*(-5 + eta**4 + beta*(-1 + eta**4))*xi**2)/(3 + alpha) - ((5 + beta)*(-3 + eta**2 + beta*(-1 + eta**2))*xi**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
   end function fun_s22one

   subroutine int_s22one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the S22 integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

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

      allocate (val_max_max(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_min_min(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_max_min(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_min_max(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                diff(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1))

      ! Calculate the value of the result integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_eta, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            diff(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 2) - j1
                  beta = size(prod_eta, 2) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R**3*knoteta(i2)**(1 + beta)*knotxi(i1)**(1 + alpha)*(((1 + beta)*knoteta(i2)**2*(5 + beta - (3 + beta)*knoteta(i2)**2))/(1 + alpha) + ((3 + beta)*(-5 + knoteta(i2)**4 + beta*(-1 + knoteta(i2)**4))*knotxi(i1)**2)/(3 + alpha) - ((5 + beta)*(-3 + knoteta(i2)**2 + beta*(-1 + knoteta(i2)**2))*knotxi(i1)**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R**3*knoteta(i2+1)**(1 + beta)*knotxi(i1+1)**(1 + alpha)*(((1 + beta)*knoteta(i2+1)**2*(5 + beta - (3 + beta)*knoteta(i2+1)**2))/(1 + alpha) + ((3 + beta)*(-5 + knoteta(i2+1)**4 + beta*(-1 + knoteta(i2+1)**4))*knotxi(i1+1)**2)/(3 + alpha) - ((5 + beta)*(-3 + knoteta(i2+1)**2 + beta*(-1 + knoteta(i2+1)**2))*knotxi(i1+1)**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R**3*knoteta(i2+1)**(1 + beta)*knotxi(i1)**(1 + alpha)*(((1 + beta)*knoteta(i2+1)**2*(5 + beta - (3 + beta)*knoteta(i2+1)**2))/(1 + alpha) + ((3 + beta)*(-5 + knoteta(i2+1)**4 + beta*(-1 + knoteta(i2+1)**4))*knotxi(i1)**2)/(3 + alpha) - ((5 + beta)*(-3 + knoteta(i2+1)**2 + beta*(-1 + knoteta(i2+1)**2))*knotxi(i1)**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(2*mppi()*R**3*knoteta(i2)**(1 + beta)*knotxi(i1+1)**(1 + alpha)*(((1 + beta)*knoteta(i2)**2*(5 + beta - (3 + beta)*knoteta(i2)**2))/(1 + alpha) + ((3 + beta)*(-5 + knoteta(i2)**4 + beta*(-1 + knoteta(i2)**4))*knotxi(i1+1)**2)/(3 + alpha) - ((5 + beta)*(-3 + knoteta(i2)**2 + beta*(-1 + knoteta(i2)**2))*knotxi(i1+1)**4)/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
            result = result + val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

   end subroutine int_s22one

   function fun_c11one(xi, eta, Z1, Z2, m, C, R, alpha, beta) result(c11one)
      !> @brief This subroutine calculates the c11one function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: jz2, alpha, beta
      type(mp_real) :: c11one

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c11one = (-2*mppi()*R*eta**(1 + beta)*xi**(1 + alpha)*(((1 + beta)*eta*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*eta))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*xi)/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*xi**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
   end function fun_c11one

   subroutine int_C11one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C11one integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return S11 : real : the value of the S11 integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

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

      allocate (val_max_max(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_min_min(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_max_min(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_min_max(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                diff(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1))

      ! Calculate the value of the S11 integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            diff(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 2) - j1
                  beta = size(prod_eta, 2) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2*mppi()*R*knoteta(i2)**(1 + beta)*knotxi(i1)**(1 + alpha)*(((1 + beta)*knoteta(i2)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knoteta(i2)))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*knotxi(i1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knotxi(i1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2*mppi()*R*knoteta(i2+1)**(1 + beta)*knotxi(i1+1)**(1 + alpha)*(((1 + beta)*knoteta(i2+1)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knoteta(i2+1)))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*knotxi(i1+1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knotxi(i1+1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2*mppi()*R*knoteta(i2+1)**(1 + beta)*knotxi(i1)**(1 + alpha)*(((1 + beta)*knoteta(i2+1)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knoteta(i2+1)))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*knotxi(i1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knotxi(i1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2*mppi()*R*knoteta(i2)**(1 + beta)*knotxi(i1+1)**(1 + alpha)*(((1 + beta)*knoteta(i2)*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*knoteta(i2)))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*knotxi(i1+1))/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*knotxi(i1+1)**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
            result = result + val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

   end subroutine int_C11one

   function fun_c11two(xi, eta, Z1, Z2, m, C, R, alpha, beta) result(c11two)
      !> @brief This subroutine calculates the c11two function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: jz2, alpha, beta
      type(mp_real) :: c11two

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c11two = (2*mppi()*R*eta**(1 + beta)*xi**(1 + alpha)*(((1 + beta)*eta*((Z1 - Z2)*(3 + beta) + c**2*m*R*(2 + beta)*eta))/(1 + alpha) - ((Z1 + Z2)*(2 + beta)*(3 + beta)*xi)/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*xi**2)/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
   end function fun_c11two

   subroutine int_C11two(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C11two integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

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

      allocate (val_max_max(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_min_min(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_max_min(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                val_min_max(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1), &
                diff(size(b_i_xi, 1)-1, size(b_i_eta, 1)-1))

      ! Calculate the value of the integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            diff(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 2) - j1
                  beta = size(prod_eta, 2) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11two(knotxi(i1), knoteta(i2), Z1, Z2, m, C, R, alpha, beta)
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11two(knotxi(i1+1), knoteta(i2+1), Z1, Z2, m, C, R, alpha, beta)
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11two(knotxi(i1), knoteta(i2+1), Z1, Z2, m, C, R, alpha, beta)
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11two(knotxi(i1+1), knoteta(i2), Z1, Z2, m, C, R, alpha, beta)
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
            result = result + val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

   end subroutine int_C11two

   function fun_c22one(xi, eta, Z1, Z2, m, C, R, alpha, beta) result(c22one)
      !> @brief This subroutine calculates the C22one function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: jz2, alpha, beta
      type(mp_real) :: c22one

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c22one = 2*mppi()*R**2*eta**(1 + beta)*xi**(1 + alpha)*((eta*((-Z1 + Z2)/(2 + beta) + (c**2*m*R*eta)/(3 + beta) + ((Z1 - Z2)*eta**2)/(4 + beta) - (c**2*m*R*eta**3)/(5 + beta)))/(1 + alpha) - ((Z1 + Z2)*(-3 + eta**2 + beta*(-1 + eta**2))*xi)/((2 + alpha)*(1 + beta)*(3 + beta)) + (((Z1 - Z2)*eta*(1/(2 + beta) - eta**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + eta**4/(5 + beta)))*xi**2)/(3 + alpha) + ((Z1 + Z2)*(-3 + eta**2 + beta*(-1 + eta**2))*xi**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*eta**2)*xi**4)/((5 + alpha)*(1 + beta)*(3 + beta)))

   end function fun_c22one

   subroutine int_C22one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C22one integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta
      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)+size(b_j_xi, 2)-1), prod_eta(size(b_i_xi, 1), size(b_i_eta, 2)+size(b_j_eta, 2)-1))

      ! Calculate the product of the B-spline coefficients
      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      ! Calculate the value of the S11 integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_eta, 1) - 1 ! Loop over the polynomials of eta
            val_min_min = zero
            val_max_max = zero
            val_min_max = zero
            val_max_min = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 2) - j1
                  beta = size(prod_eta, 2) - j2
                  val_min_min = val_min_min + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22one(knotxi(i1), knoteta(i2), Z1, Z2, m, C, R, alpha, beta)
                  val_max_max = val_max_max + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22one(knotxi(i1+1), knoteta(i2+1), Z1, Z2, m, C, R, alpha, beta)
                  val_min_max = val_min_max + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22one(knotxi(i1), knoteta(i2+1), Z1, Z2, m, C, R, alpha, beta)
                  val_max_min = val_max_min + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22one(knotxi(i1+1), knoteta(i2), Z1, Z2, m, C, R, alpha, beta)
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

      ! ! Test
      ! print *, "Begin test "
      ! do i1 = 1, size(b_i_xi, 1) - 1
      !    call write_lists(diff(i1, :), 6, 30, 10)
      ! end do
      ! print *, "End test "
      ! ! Test end

   end subroutine int_C22one

   function fun_c22two(xi, eta, Z1, Z2, m, C, R, alpha, beta) result(c22two)
      !> @brief This subroutine calculates the C22two function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: jz2, alpha, beta
      type(mp_real) :: c22two

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c22two = (-2)*mppi()*R*eta**(1 + beta)*xi**(1 + alpha)*((eta*((Z1 - Z2)/(2 + beta) + (c**2*m*R*eta)/(3 + beta) + ((-Z1 + Z2)*eta**2)/(4 + beta) - (c**2*m*R*eta**3)/(5 + beta)))/(1 + alpha) + ((Z1 + Z2)*(-3 + eta**2 + beta*(-1 + eta**2))*xi)/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*eta*(-(1/(2 + beta)) + eta**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + eta**4/(5 + beta)))*xi**2)/(3 + alpha) - ((Z1 + Z2)*(-3 + eta**2 + beta*(-1 + eta**2))*xi**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*eta**2)*xi**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
   end function fun_c22two

   subroutine int_C22two(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C22two integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta
      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)+size(b_j_xi, 2)-1), prod_eta(size(b_i_xi, 1), size(b_i_eta, 2)+size(b_j_eta, 2)-1))

      ! Calculate the product of the B-spline coefficients
      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      ! Calculate the value of the S11 integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min = zero
            val_max_max = zero
            val_min_max = zero
            val_max_min = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 2) - j1
                  beta = size(prod_eta, 2) - j2
                  val_min_min = val_min_min + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2)*mppi()*R*knoteta(i2)**(1 + beta)*knotxi(i1)**(1 + alpha)*((knoteta(i2)*((Z1 - Z2)/(2 + beta) + (c**2*m*R*knoteta(i2))/(3 + beta) + ((-Z1 + Z2)*knoteta(i2)**2)/(4 + beta) - (c**2*m*R*knoteta(i2)**3)/(5 + beta)))/(1 + alpha) + ((Z1 + Z2)*(-3 + knoteta(i2)**2 + beta*(-1 + knoteta(i2)**2))*knotxi(i1))/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*knoteta(i2)*(-(1/(2 + beta)) + knoteta(i2)**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + knoteta(i2)**4/(5 + beta)))*knotxi(i1)**2)/(3 + alpha) - ((Z1 + Z2)*(-3 + knoteta(i2)**2 + beta*(-1 + knoteta(i2)**2))*knotxi(i1)**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*knoteta(i2)**2)*knotxi(i1)**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
                  val_max_max = val_max_max + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2)*mppi()*R*knoteta(i2+1)**(1 + beta)*knotxi(i1+1)**(1 + alpha)*((knoteta(i2+1)*((Z1 - Z2)/(2 + beta) + (c**2*m*R*knoteta(i2+1))/(3 + beta) + ((-Z1 + Z2)*knoteta(i2+1)**2)/(4 + beta) - (c**2*m*R*knoteta(i2+1)**3)/(5 + beta)))/(1 + alpha) + ((Z1 + Z2)*(-3 + knoteta(i2+1)**2 + beta*(-1 + knoteta(i2+1)**2))*knotxi(i1+1))/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*knoteta(i2+1)*(-(1/(2 + beta)) + knoteta(i2+1)**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + knoteta(i2+1)**4/(5 + beta)))*knotxi(i1+1)**2)/(3 + alpha) - ((Z1 + Z2)*(-3 + knoteta(i2+1)**2 + beta*(-1 + knoteta(i2+1)**2))*knotxi(i1+1)**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*knoteta(i2+1)**2)*knotxi(i1+1)**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
                  val_min_max = val_min_max + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2)*mppi()*R*knoteta(i2+1)**(1 + beta)*knotxi(i1)**(1 + alpha)*((knoteta(i2+1)*((Z1 - Z2)/(2 + beta) + (c**2*m*R*knoteta(i2+1))/(3 + beta) + ((-Z1 + Z2)*knoteta(i2+1)**2)/(4 + beta) - (c**2*m*R*knoteta(i2+1)**3)/(5 + beta)))/(1 + alpha) + ((Z1 + Z2)*(-3 + knoteta(i2+1)**2 + beta*(-1 + knoteta(i2+1)**2))*knotxi(i1))/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*knoteta(i2+1)*(-(1/(2 + beta)) + knoteta(i2+1)**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + knoteta(i2+1)**4/(5 + beta)))*knotxi(i1)**2)/(3 + alpha) - ((Z1 + Z2)*(-3 + knoteta(i2+1)**2 + beta*(-1 + knoteta(i2+1)**2))*knotxi(i1)**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*knoteta(i2+1)**2)*knotxi(i1)**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
                  val_max_min = val_max_min + prod_xi(i1, j1)*prod_eta(i2, j2)*(-2)*mppi()*R*knoteta(i2)**(1 + beta)*knotxi(i1+1)**(1 + alpha)*((knoteta(i2)*((Z1 - Z2)/(2 + beta) + (c**2*m*R*knoteta(i2))/(3 + beta) + ((-Z1 + Z2)*knoteta(i2)**2)/(4 + beta) - (c**2*m*R*knoteta(i2)**3)/(5 + beta)))/(1 + alpha) + ((Z1 + Z2)*(-3 + knoteta(i2)**2 + beta*(-1 + knoteta(i2)**2))*knotxi(i1+1))/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*knoteta(i2)*(-(1/(2 + beta)) + knoteta(i2)**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + knoteta(i2)**4/(5 + beta)))*knotxi(i1+1)**2)/(3 + alpha) - ((Z1 + Z2)*(-3 + knoteta(i2)**2 + beta*(-1 + knoteta(i2)**2))*knotxi(i1+1)**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*knoteta(i2)**2)*knotxi(i1+1)**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

   end subroutine int_C22two

   function fun_c11three(xi, eta, Z1, Z2, m, C, R, alpha, beta, chi, delta) result(c11three)
      !> @brief This subroutine calculates the C11three function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: alpha, beta, chi, delta
      type(mp_real) :: c11three

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c11three = (2*mppi()*eta**(beta + chi)*xi**(alpha + delta)*(-((delta*eta**2*(beta + chi))/(alpha + delta)) + (xi**2*(2*chi + (beta + chi)*(eta**2*(delta - chi) + chi)))/(2 + alpha + delta)))/((beta + chi)*(2 + beta + chi))
   end function fun_c11three

   subroutine int_C11three(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C11three integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max
      integer :: i1, i2, j1, j2, j3, j4, alpha, beta, chi, delta

      zero = '0.0d0'

      ! Calculate the value of the S11 integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min = zero
            val_max_max = zero
            val_min_max = zero
            val_max_min = zero
            do j1 = 1, size(b_i_xi, 2) ! Loop over the order of xi i
               do j2 = 1, size(b_i_eta, 2) ! Loop over the order of eta i
                  do j3 = 1, size(b_j_xi, 2) ! Loop over the order of xi j
                     do j4 = 1, size(b_j_eta, 2) ! Loop over the order of eta j
                        alpha = size(b_i_xi, 2) - j1
                        beta = size(b_i_eta, 2) - j2
                        chi = size(b_j_xi, 2) - j3
                        delta = size(b_j_eta, 2) - j4
                        val_min_min = val_min_min + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*fun_c11three(knotxi(i1), knoteta(i2), Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_max = val_max_max + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*fun_c11three(knotxi(i1+1), knoteta(i2+1), Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_min_max = val_min_max + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*fun_c11three(knotxi(i1), knoteta(i2+1), Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_min = val_max_min + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*fun_c11three(knotxi(i1+1), knoteta(i2), Z1, Z2, m, C, R, alpha, beta, chi, delta)
                     end do
                  end do
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

   end subroutine int_C11three

   function fun_c22three(xi, eta, Z1, Z2, m, C, R, alpha, beta, chi, delta) result(c22three)
      !> @brief This subroutine calculates the C22three function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: alpha, beta, chi, delta
      type(mp_real) :: c22three

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      if (beta ==0 .and. chi == 0) then
         c22three = zero
      else
         c22three = (2*mppi()*R**2*eta**(beta + chi)*xi**(alpha + delta)*((eta**2*xi**2*(2 + delta + delta*eta**2 - 2*(1 + delta)*eta**4 + (-1 + eta**2)**2*chi - (2*(2 + beta - delta))/(2 + beta + chi) + (2*(6 + beta + 2*delta)*eta**4)/(6 + beta + chi)))/(2 + alpha + delta) + (delta*eta**4*(-1 + (eta**2*(4 + beta + chi))/(6 + beta + chi)))/(alpha + delta)-(xi**4*(-((chi*(4 + beta + chi))/(beta + chi)) - (eta**2*(2*delta - chi)*(4 + beta + chi))/(2 + beta + chi) + (eta**6*(delta - chi)*(4 + beta + chi))/(6 + beta + chi) + eta**4*(delta + chi)))/(4 + alpha + delta)-(xi**6*(8*chi + (beta + chi)*(delta*eta**2*(4 + beta - (2 + beta)*eta**2) - (-1 + eta**2)*(6 + beta + (-2 - beta + delta)*eta**2)*chi + (-1 + eta**2)**2*chi**2)))/((6 + alpha + delta)*(beta + chi)*(2 + beta + chi))))/(4 + beta + chi)
      end if
   end function fun_c22three

   subroutine int_C22three(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C22three integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max
      integer :: i1, i2, j1, j2, j3, j4, alpha, beta, chi, delta

      zero = '0.0d0'

      ! Calculate the value of the S11 integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min = zero
            val_max_max = zero
            val_min_max = zero
            val_max_min = zero
            do j1 = 1, size(b_i_xi, 2) ! Loop over the order of xi i
               do j2 = 1, size(b_i_eta, 2) ! Loop over the order of eta i
                  do j3 = 1, size(b_j_xi, 2) ! Loop over the order of xi j
                     do j4 = 1, size(b_j_eta, 2) ! Loop over the order of eta j
                        alpha = size(b_i_xi, 2) - j1
                        beta = size(b_i_eta, 2) - j2
                        chi = size(b_j_xi, 2) - j3
                        delta = size(b_j_eta, 2) - j4
                        val_min_min = val_min_min + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*fun_c22three(knotxi(i1), knoteta(i2), Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_max = val_max_max + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*fun_c22three(knotxi(i1+1), knoteta(i2+1), Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_min_max = val_min_max + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*fun_c22three(knotxi(i1), knoteta(i2+1), Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_min = val_max_min + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*fun_c22three(knotxi(i1+1), knoteta(i2), Z1, Z2, m, C, R, alpha, beta, chi, delta)
                     end do
                  end do
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

   end subroutine int_C22three

   function fun_c12three(xi, eta, Z1, Z2, m, C, R, alpha, beta, chi, delta) result(c12three)
      !> @brief This subroutine calculates the C12three function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: jz2, alpha, beta, chi, delta
      type(mp_real) :: c12three

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c12three = 2*mppi()*R**2*eta**(1 + beta + chi)*xi**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + xi**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + eta**2/(3 + beta + chi))
   end function fun_c12three

   subroutine int_C12three(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C12three integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max
      integer :: i, i1, i2, j1, j2, j3, j4, alpha, beta, chi, delta

      zero = '0.0d0'

      ! Calculate the value of the S11 integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min = zero
            val_max_max = zero
            val_min_max = zero
            val_max_min = zero
            do j1 = 1, size(b_i_xi, 2) ! Loop over the order of xi i
               do j2 = 1, size(b_i_eta, 2) ! Loop over the order of eta i
                  do j3 = 1, size(b_j_xi, 2) ! Loop over the order of xi j
                     do j4 = 1, size(b_j_eta, 2) ! Loop over the order of eta j
                        alpha = size(b_i_xi, 2) - j1
                        beta = size(b_i_eta, 2) - j2
                        chi = size(b_j_xi, 2) - j3
                        delta = size(b_j_eta, 2) - j4
                        val_min_min = val_min_min + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*2*mppi()*R**2*knoteta(i2)**(1 + beta + chi)*knotxi(i1)**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + knotxi(i1)**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + knoteta(i2)**2/(3 + beta + chi))
                        val_max_max = val_max_max + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*2*mppi()*R**2*knoteta(i2+1)**(1 + beta + chi)*knotxi(i1+1)**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + knotxi(i1+1)**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + knoteta(i2+1)**2/(3 + beta + chi))
                        val_min_max = val_min_max + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*2*mppi()*R**2*knoteta(i2+1)**(1 + beta + chi)*knotxi(i1)**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + knotxi(i1)**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + knoteta(i2+1)**2/(3 + beta + chi))
                        val_max_min = val_max_min + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*2*mppi()*R**2*knoteta(i2)**(1 + beta + chi)*knotxi(i1+1)**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + knotxi(i1+1)**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + knoteta(i2)**2/(3 + beta + chi))
                     end do
                  end do
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

   end subroutine int_C12three

   function fun_c21three(xi, eta, Z1, Z2, m, C, R, alpha, beta, chi, delta) result(c21three)
      !> @brief This subroutine calculates the C21three function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: xi, eta, m, C, R, Z1, Z2
      integer :: jz2, alpha, beta, chi, delta
      type(mp_real) :: c21three

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c21three = 2*mppi()*R**2*eta**(1 + beta + chi)*xi**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + xi**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + eta**2/(3 + beta + chi))
   end function fun_c21three

   subroutine int_C21three(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C21three integral for the H2+ molecule.
      !> @param b_i_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_i_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param b_j_xi : real(:, :) : the B-spline coefficients for the xi direction
      !> @param b_j_eta : real(:, :) : the B-spline coefficients for the eta direction
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param knotxi : real(:) : the knot vector for the xi direction
      !> @param knoteta : real(:) : the knot vector for the eta direction
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei
      !> @param jz2 : real : the quantum number (2*jz)
      !> @return result : real : the value of the result integral
      type(mp_real), intent(in) :: Z1, Z2, m, C, R
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
      type(mp_real), dimension(:, :), intent(in) :: b_i_xi, b_i_eta, b_j_xi, b_j_eta

      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max
      integer :: i1, i2, j1, j2, j3, j4, alpha, beta, chi, delta

      zero = '0.0d0'

      ! Calculate the value of the S11 integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min = zero
            val_max_max = zero
            val_min_max = zero
            val_max_min = zero
            do j1 = 1, size(b_i_xi, 2) ! Loop over the order of xi i
               do j2 = 1, size(b_i_eta, 2) ! Loop over the order of eta i
                  do j3 = 1, size(b_j_xi, 2) ! Loop over the order of xi j
                     do j4 = 1, size(b_j_eta, 2) ! Loop over the order of eta j
                        alpha = size(b_i_xi, 2) - j1
                        beta = size(b_i_eta, 2) - j2
                        chi = size(b_j_xi, 2) - j3
                        delta = size(b_j_eta, 2) - j4
                        val_min_min = val_min_min + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*2*mppi()*R**2*knoteta(i2)**(1 + beta + chi)*knotxi(i1)**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + knotxi(i1)**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + knoteta(i2)**2/(3 + beta + chi))
                        val_max_max = val_max_max + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*2*mppi()*R**2*knoteta(i2+1)**(1 + beta + chi)*knotxi(i1+1)**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + knotxi(i1+1)**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + knoteta(i2+1)**2/(3 + beta + chi))
                        val_min_max = val_min_max + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*2*mppi()*R**2*knoteta(i2+1)**(1 + beta + chi)*knotxi(i1)**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + knotxi(i1)**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + knoteta(i2+1)**2/(3 + beta + chi))
                        val_max_min = val_max_min + b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)*2*mppi()*R**2*knoteta(i2)**(1 + beta + chi)*knotxi(i1+1)**(1 + alpha + delta)*(-(1/(1 + alpha + delta)) + knotxi(i1+1)**2/(3 + alpha + delta))*(-delta + chi)*(-(1/(1 + beta + chi)) + knoteta(i2)**2/(3 + beta + chi))
                     end do
                  end do
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do
      
   end subroutine int_C21three

   subroutine init_h2plus(d, n, n_remove, Z1, Z2, m, C, R, ximax, ximin, jz2, epsilon, eta_slp, save_step)
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

      type(mp_real), dimension(:), allocatable :: knotxi, knoteta, knotxi_eps, knoteta_eps
      type(mp_real), dimension(:, :), allocatable :: S11one, S22one, C11one, C11two, C22one, C22two, C11three, C22three,C12three, C21three, C_mat, S_mat, vect
      type(mp_real), dimension(:), allocatable :: w, fv1, fv2
      type(mp_real), dimension(:, :, :), allocatable :: bspline_xi, bspline_eta
      logical :: debug_bool = .true.

      integer :: ntot, i, j, ierr
      integer, dimension(2) :: i2, j2

      ntot = n + d + 2*n_remove
      zero = '0.d0'

      write (6, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d, " and Number of BSplines to remove: ", n_remove

      print *, "Generating B-spline knots..."
      ! Generate the knot vectors for xi and eta
      allocate (knotxi(ntot), knoteta(ntot))

      knotxi = knot_xi(d, n, n_remove, ximin, ximax)
      knoteta = knot_eta(d, n, n_remove, eta_slp)

      print *, "Generating B-spline coefficients..."
      ! Generate the B-spline coefficients for xi and eta
      allocate (bspline_xi(n, ntot, d), bspline_eta(n, ntot, d))

      do i = 1 + n_remove, n + n_remove ! Loop over the number of B-splines
         if (debug_bool) then
            print *, "Generating B-spline coefficients for B-spline ", i, "on xi-axis..."
         end if
         call init_bspine(d, i, knotxi, bspline_xi(i-n_remove, :, :), debug_bool)
         if (debug_bool) then
            print *, "Generating B-spline coefficients for B-spline ", i, "on eta-axis..."
         end if
         call init_bspine(d, i, knoteta, bspline_eta(i-n_remove, :, :), debug_bool)
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

      print *, "Calculating the S11 integral..."
      ! Calculate the S11 integral
      allocate (S11one(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, S11one, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            S11one(i, j) = zero
            call int_s11one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, S11one(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving S11one matrix to file..."
         open (unit=2, file='S11one.csv', status='replace')
         do i = 1, n**2
            call write_csv(S11one(i, :), 2, 30, 10)
         end do
         close (2)
      end if

      print *, "Calculating the S22 integral..."
      ! Calculate the S22 integral
      allocate (S22one(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, S22one, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            S22one(i, j) = zero
            call int_s22one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, S22one(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving S22one matrix to file..."
         open (unit=3, file='S22one.csv', status='replace')
         do i = 1, n**2
            call write_csv(S22one(i, :), 3, 30, 10)
         end do
         close (3)
      end if

      print *, "Calculating the C11one integral..."
      ! Calculate the C11one integral
      allocate (C11one(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C11one, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C11one(i, j) = zero
            call int_C11one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, C11one(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C11one matrix to file..."
         open (unit=4, file='C11one.csv', status='replace')
         do i = 1, n**2
            call write_csv(C11one(i, :), 4, 30, 10)
         end do
         close (4)
      end if

      print *, "Calculating the C11two integral..."
      ! Calculate the C11two integral
      allocate (C11two(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C11two, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C11two(i, j) = zero
            call int_C11two(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, C11two(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C11two matrix to file..."
         open (unit=5, file='C11two.csv', status='replace')
         do i = 1, n**2
            call write_csv(C11two(i, :), 5, 30, 10)
         end do
         close (5)
      end if

      print *, "Calculating the C22one integral..."
      ! Calculate the C22one integral
      allocate (C22one(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C22one, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C22one(i, j) = zero
            call int_C22one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, C22one(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C22one matrix to file..."
         open (unit=16, file='C22one.csv', status='replace')
         do i = 1, n**2
            call write_csv(C22one(i, :), 16, 30, 10)
         end do
         close (16)
      end if

      print *, "Calculating the C22two integral..."
      ! Calculate the C22two integral
      allocate (C22two(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C22two, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C22two(i, j) = zero
            call int_C22two(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, C22two(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C22two matrix to file..."
         open (unit=7, file='C22two.csv', status='replace')
         do i = 1, n**2
            call write_csv(C22two(i, :), 7, 30, 10)
         end do
         close (7)
      end if

      print *, "Calculating the C11three integral..."
      ! Calculate the C11three integral
      allocate (C11three(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C11three, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C11three(i, j) = zero
            call int_C11three(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                              bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                              knotxi, knoteta, Z1, Z2, m, C, R, C11three(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C11three matrix to file..."
         open (unit=8, file='C11three.csv', status='replace')
         do i = 1, n**2
            call write_csv(C11three(i, :), 8, 30, 10)
         end do
         close (8)
      end if

      print *, "Calculating the C22three integral..."
      ! Calculate the C22three integral
      allocate (C22three(n**2, n**2))

      ! !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C22three, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C22three(i, j) = zero
            call int_C22three(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                              bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                              knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C22three(i, j))
         end do
      end do
      ! !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C22three matrix to file..."
         open (unit=9, file='C22three.csv', status='replace')
         do i = 1, n**2
            call write_csv(C22three(i, :), 9, 30, 10)
         end do
         close (9)
      end if

      print *, "Calculating the C12three integral..."
      ! Calculate the C12three integral
      allocate (C12three(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C12three, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C12three(i, j) = zero
            call int_C12three(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                              bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                              knotxi, knoteta, Z1, Z2, m, C, R, C12three(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C12three matrix to file..."
         open (unit=10, file='C12three.csv', status='replace')
         do i = 1, n**2
            call write_csv(C12three(i, :), 10, 30, 10)
         end do
         close (10)
      end if

      print *, "Calculating the C21three integral..."
      ! Calculate the C21three integral
      allocate (C21three(n**2, n**2))

      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C21three, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = 1, n**2
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            C21three(i, j) = zero
            call int_C21three(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                              bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                              knotxi, knoteta, Z1, Z2, m, C, R, C21three(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C21three matrix to file..."
         open (unit=11, file='C21three.csv', status='replace')
         do i = 1, n**2
            call write_csv(C21three(i, :), 11, 30, 10)
         end do
         close (11)
      end if

      print *, "Generating the C Matrix..."
      ! Generate the C matrix
      allocate (C_mat(4*n**2, 4*n**2))

      C_mat = zero
      do i = 1, n**2
         do j = 1, n**2
            ! Diagonal blocks
            C_mat(i, j) = C11one(i, j)
            C_mat(i + n**2, j + n**2) = C22one(i, j)
            C_mat(i + 2*n**2, j + 2*n**2) = C11two(i, j)
            C_mat(i + 3*n**2, j + 3*n**2) = C22two(i, j)
            ! Off-diagonal blocks
            C_mat(i, j + 2*n**2) = C11three(i, j)
            C_mat(i, j + 3*n**2) = C12three(i, j)
            C_mat(i + n**2, j + 2*n**2) = C21three(i, j)
            C_mat(i + n**2, j + 3*n**2) = C22three(i, j)
            C_mat(i + 2*n**2, j) = C11three(j, i)
            C_mat(i + 2*n**2, j + n**2) = C21three(j, i)
            C_mat(i + 3*n**2, j) = C12three(j, i)
            C_mat(i + 3*n**2, j + n**2) = C22three(j, i)
         end do
      end do

      print *, "Generating the S Matrix..."
      ! Generate the S matrix
      allocate (S_mat(4*n**2, 4*n**2))

      S_mat = zero
      do i = 1, n**2
         do j = 1, n**2
            ! Diagonal blocks
            S_mat(i, j) = S11one(i, j)
            S_mat(i + n**2, j + n**2) = S22one(i, j)
            S_mat(i + 2*n**2, j + 2*n**2) = S11one(i, j)
            S_mat(i + 3*n**2, j + 3*n**2) = S22one(i, j)
         end do
      end do

      if (debug_bool) then
         print *, "Check if C is hermitian..."
         do i = 1, 4*n**2
            do j = 1, 4*n**2
               if (C_mat(i, j) - C_mat(j, i) > epsilon) then
                  print *, "C is not hermitian at (", i, ",", j, ")"
               end if
            end do
         end do

         print *, "Check if S is hermitian..."
         do i = 1, 4*n**2
            do j = 1, 4*n**2
               if (S_mat(i, j) - S_mat(j, i) > epsilon) then
                  print *, "S is not hermitian at (", i, ",", j, ")"
               end if
            end do
         end do
      end if 

      print *, "Calculating the eigenvalues..."
      ! Calculate the eigenvalues and eigenvectors
      allocate (w(4*n**2), fv1(4*n**2), fv2(4*n**2))

      call rsg(4*n**2, 4*n**2, C_mat, S_mat, w, 0, vect, fv1, fv2, ierr)

      print *, "Error code: ", ierr

      print *, "Saving logs to file..."
      ! Save logs
      open (unit=1, file='log_file', status='replace')
      write (1, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d, " and Number of BSplines to remove: ", n_remove
      write (1, '(a)') "Speed of light: "
      call mpwrite(1, 30, 10, C)
      write (1, '(a)') "Mass of the electron: "
      call mpwrite(1, 30, 10, m)
      write (1, '(a)') "Slope for eta: "
      call mpwrite(1, 30, 10, eta_slp)
      write (1, '(a)') "Non-adjusted Xi knot vector: "
      call write_lists(knotxi, 1, 30, 10)
      write (1, '(a)') "Non-adjusted Eta knot vector: "
      call write_lists(knoteta, 1, 30, 10)
      write (1, '(a)') "B-Spline basis function Order: "
      do i = 1, n**2
         i2 = indexToPair(i, n)
         write (1, '(i4, a, i4, a, i4)', advance='no') i, " ", i2(1), " ", i2(2)
         write (1, '(a)') " "
      end do
      write (1, '(a)') "Adjusted Xi knot vector: "
      call write_lists(knotxi_eps, 1, 30, 10)
      write (1, '(a)') "Adjusted Eta knot vector: "
      call write_lists(knoteta_eps, 1, 30, 10)
      if (debug_bool) then
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
         write (1, '(a)') "C22two integral: "
         do i = 1, n**2
            call write_lists(C22two(i, :), 1, 30, 10)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C11three integral: "
         do i = 1, n**2
            call write_lists(C11three(i, :), 1, 30, 10)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C22three integral: "
         do i = 1, n**2
            call write_lists(C22three(i, :), 1, 30, 10)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C12three integral: "
         do i = 1, n**2
            call write_lists(C12three(i, :), 1, 30, 10)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C21three integral: "
         do i = 1, n**2
            call write_lists(C21three(i, :), 1, 30, 10)
         end do
      end if
      write (1, '(a)') "--------------------------------------------------------------"
      write (1, '(a)') "Eigenvalues: "
      do i = 1, 4*n**2
         write (1, '(i4, a, i4)', advance='no') i, " "
         call mpwrite(1, 30, 10, w(i)-m*c*c)
      end do
      close (1)

   end subroutine init_h2plus

end module h2plus
