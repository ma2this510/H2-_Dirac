module h2plus
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none

   type(mp_real), save :: one, zero

   type(mp_real), dimension(:, :), allocatable, save :: pwr_knot_xi, pwr_knot_eta

   private

   public :: init_h2plus

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

   function fun_s11one(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta) result(s11one)
      !> @brief This subroutine calculates the S11one function for the H2+ molecule.
      !> @param xi_i : integer : th i-th xi knot
      !> @param eta_i : integer : the i-th eta knot
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) ::  m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta
      type(mp_real) :: s11one

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      s11one = 2*mppi()*(R**3)*(get_pw_xi(xi_i, alpha + 1))*(get_pw_eta(eta_i,beta + 1))*(-(get_pw_eta(eta_i,2))/((alpha + 1)*(beta + 3)) + (get_pw_xi(xi_i,2))/((alpha + 3)*(beta + 1)))
   end function fun_s11one

   subroutine int_s11one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the S11 integral for the H2+ molecule.
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

      allocate (val_max_max(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_min_min(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_max_min(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_min_max(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                diff(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1))

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
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_s11one(i1, i2, Z1, Z2, m, C, R, alpha, beta)
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_s11one(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta)
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_s11one(i1, i2+1, Z1, Z2, m, C, R, alpha, beta)
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_s11one(i1+1, i2, Z1, Z2, m, C, R, alpha, beta)
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
            result = result + val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

   end subroutine int_s11one

   function fun_s22one(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta) result(s22one)
      !> @brief This subroutine calculates the S22one function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta
      type(mp_real) :: s22one

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      s22one = (2*mppi()*R**3*get_pw_eta(eta_i, 1 + beta)*get_pw_xi(xi_i, 1 + alpha)*(((1 + beta)*get_pw_eta(eta_i, 2)*(5 + beta - (3 + beta)*get_pw_eta(eta_i, 2)))/(1 + alpha) + ((3 + beta)*(-5 + get_pw_eta(eta_i, 4) + beta*(-1 + get_pw_eta(eta_i, 4)))*get_pw_xi(xi_i, 2))/(3 + alpha) - ((5 + beta)*(-3 + get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))*get_pw_xi(xi_i, 4))/(5 + alpha)))/((1 + beta)*(3 + beta)*(5 + beta))
   end function fun_s22one

   subroutine int_s22one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the S22 integral for the H2+ molecule.
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

      allocate (val_max_max(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_min_min(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_max_min(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_min_max(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                diff(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1))

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
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_s22one(i1, i2, Z1, Z2, m, C, R, alpha, beta)
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_s22one(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta)
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_s22one(i1, i2+1, Z1, Z2, m, C, R, alpha, beta)
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_s22one(i1+1, i2, Z1, Z2, m, C, R, alpha, beta)
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
            result = result + val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

   end subroutine int_s22one

   function fun_c11one(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta) result(c11one)
      !> @brief This subroutine calculates the c11one function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta
      type(mp_real) :: c11one

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c11one = (-2*mppi()*R**2*get_pw_eta(eta_i, 1 + beta)*get_pw_xi(xi_i, 1 + alpha)*(((1 + beta)*eta_i*(-(Z1*(3 + beta)) + Z2*(3 + beta) + c**2*m*R*(2 + beta)*eta_i))/(1 + alpha) + ((Z1 + Z2)*(2 + beta)*(3 + beta)*xi_i)/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*get_pw_xi(xi_i, 2))/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
   end function fun_c11one

   subroutine int_C11one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C11one integral for the H2+ molecule.
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

      allocate (val_max_max(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_min_min(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_max_min(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_min_max(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                diff(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1))

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
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11one(i1, i2, Z1, Z2, m, C, R, alpha, beta)
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11one(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta)
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11one(i1, i2+1, Z1, Z2, m, C, R, alpha, beta)
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11one(i1+1, i2, Z1, Z2, m, C, R, alpha, beta)
               end do
            end do
            diff(i1, i2) = val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
            result = result + val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

   end subroutine int_C11one

   function fun_c11two(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta) result(c11two)
      !> @brief This subroutine calculates the c11two function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta
      type(mp_real) :: c11two

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c11two = (2*mppi()*R**2*get_pw_eta(eta_i, 1 + beta)*get_pw_xi(xi_i, 1 + alpha)*(((1 + beta)*eta_i*((Z1 - Z2)*(3 + beta) + c**2*m*R*(2 + beta)*eta_i))/(1 + alpha) - ((Z1 + Z2)*(2 + beta)*(3 + beta)*xi_i)/(2 + alpha) - (c**2*m*R*(2 + beta)*(3 + beta)*get_pw_xi(xi_i, 2))/(3 + alpha)))/((1 + beta)*(2 + beta)*(3 + beta))
   end function fun_c11two

   subroutine int_C11two(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C11two integral for the H2+ molecule.
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

      type(mp_real), dimension(:, :), allocatable :: prod_xi, prod_eta, val_max_max, val_min_min, val_max_min, val_min_max
      integer :: i, i1, i2, j1, j2, alpha, beta

      zero = '0.0d0'

      allocate (prod_xi(size(b_i_xi, 1), size(b_i_xi, 2)+size(b_j_xi, 2)-1), prod_eta(size(b_i_xi, 1), size(b_i_eta, 2)+size(b_j_eta, 2)-1))

      ! Calculate the product of the B-spline coefficients
      do i = 1, size(b_i_xi, 1) ! Loop over the polynomials
         ! Calculate the product of the B-spline coefficients
         prod_xi(i, :) = fusion_coef(b_i_xi(i, :), b_j_xi(i, :))
         prod_eta(i, :) = fusion_coef(b_i_eta(i, :), b_j_eta(i, :))
      end do

      allocate (val_max_max(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_min_min(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_max_min(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1), &
                val_min_max(size(b_i_xi, 1) - 1, size(b_i_eta, 1) - 1))

      ! Calculate the value of the integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min(i1, i2) = zero
            val_max_max(i1, i2) = zero
            val_min_max(i1, i2) = zero
            val_max_min(i1, i2) = zero
            do j1 = 1, size(prod_xi, 2) ! Loop over the order of xi
               do j2 = 1, size(prod_eta, 2) ! Loop over the order of eta
                  alpha = size(prod_xi, 2) - j1
                  beta = size(prod_eta, 2) - j2
                  val_min_min(i1, i2) = val_min_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11two(i1, i2, Z1, Z2, m, C, R, alpha, beta)
                  val_max_max(i1, i2) = val_max_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11two(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta)
                  val_min_max(i1, i2) = val_min_max(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11two(i1, i2+1, Z1, Z2, m, C, R, alpha, beta)
                  val_max_min(i1, i2) = val_max_min(i1, i2) + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c11two(i1+1, i2, Z1, Z2, m, C, R, alpha, beta)
               end do
            end do
            result = result + val_max_max(i1, i2) + val_min_min(i1, i2) - val_max_min(i1, i2) - val_min_max(i1, i2)
         end do
      end do

   end subroutine int_C11two

   function fun_c22one(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta) result(c22one)
      !> @brief This subroutine calculates the C22one function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta
      type(mp_real) :: c22one

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      ! c22one = 2*mppi()*R**2*eta**(1 + beta)*xi**(1 + alpha)*((eta*((-Z1 + Z2)/(2 + beta) + (c**2*m*R*eta)/(3 + beta) + ((Z1 - Z2)*eta**2)/(4 + beta) - (c**2*m*R*eta**3)/(5 + beta)))/(1 + alpha) - ((Z1 + Z2)*(-3 + eta**2 + beta*(-1 + eta**2))*xi)/((2 + alpha)*(1 + beta)*(3 + beta)) + (((Z1 - Z2)*eta*(1/(2 + beta) - eta**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + eta**4/(5 + beta)))*xi**2)/(3 + alpha) + ((Z1 + Z2)*(-3 + eta**2 + beta*(-1 + eta**2))*xi**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*eta**2)*xi**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
      c22one = 2*mppi()*R**2*get_pw_eta(eta_i, 1 + beta)*get_pw_xi(xi_i, 1 + alpha)*((eta_i*(-(c**2*m*R*(8 + 6*beta + beta**2)*eta_i*(-5 + 3*get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))) + Z1*(15 + 8*beta + beta**2)*(2*(-2 + get_pw_eta(eta_i, 2)) + beta*(-1 + get_pw_eta(eta_i, 2))) - Z2*(15 + 8*beta + beta**2)*(2*(-2 + get_pw_eta(eta_i, 2)) + beta*(-1 + get_pw_eta(eta_i, 2)))))/((1 + alpha)*(2 + beta)*(3 + beta)*(4 + beta)*(5 + beta)) - ((Z1 + Z2)*(-3 + get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))*xi_i)/((2 + alpha)*(1 + beta)*(3 + beta)) + ((-((Z1 - Z2)*(5 + 6*beta + beta**2)*eta_i*(2*(-2 + get_pw_eta(eta_i, 2)) + beta*(-1 + get_pw_eta(eta_i, 2)))) + c**2*m*R*(8 + 6*beta + beta**2)*(-5 + get_pw_eta(eta_i, 4) + beta*(-1 + get_pw_eta(eta_i, 4))))*get_pw_xi(xi_i, 2))/((3 + alpha)*(1 + beta)*(2 + beta)*(4 + beta)*(5 + beta)) + ((Z1 + Z2)*(-3 + get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))*get_pw_xi(xi_i, 3))/((4 + alpha)*(1 + beta)*(3 + beta)) - (c**2*m*R*(-3 + get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))*get_pw_xi(xi_i, 4))/((5 + alpha)*(1 + beta)*(3 + beta)))

   end function fun_c22one

   subroutine int_C22one(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C22one integral for the H2+ molecule.
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
      val_min_min = val_min_min + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22one(i1, i2, Z1, Z2, m, C, R, alpha, beta)
  val_max_max = val_max_max + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22one(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta)
  val_min_max = val_min_max + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22one(i1, i2+1, Z1, Z2, m, C, R, alpha, beta)
  val_max_min = val_max_min + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22one(i1+1, i2, Z1, Z2, m, C, R, alpha, beta)
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

   end subroutine int_C22one

   function fun_c22two(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta) result(c22two)
      !> @brief This subroutine calculates the C22two function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta
      type(mp_real) :: c22two

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      ! c22two = (-2)*mppi()*R*eta**(1 + beta)*xi**(1 + alpha)*((eta*((Z1 - Z2)/(2 + beta) + (c**2*m*R*eta)/(3 + beta) + ((-Z1 + Z2)*eta**2)/(4 + beta) - (c**2*m*R*eta**3)/(5 + beta)))/(1 + alpha) + ((Z1 + Z2)*(-3 + eta**2 + beta*(-1 + eta**2))*xi)/((2 + alpha)*(1 + beta)*(3 + beta)) - (((Z1 - Z2)*eta*(-(1/(2 + beta)) + eta**2/(4 + beta)) + c**2*m*R*(-(1/(1 + beta)) + eta**4/(5 + beta)))*xi**2)/(3 + alpha) - ((Z1 + Z2)*(-3 + eta**2 + beta*(-1 + eta**2))*xi**3)/((4 + alpha)*(1 + beta)*(3 + beta)) + (c**2*m*R*(3 + beta - (1 + beta)*eta**2)*xi**4)/((5 + alpha)*(1 + beta)*(3 + beta)))
      c22two = -2*mppi()*R**2*get_pw_eta(eta_i, 1 + beta)*get_pw_xi(xi_i, 1 + alpha)*((eta_i*(-(c**2*m*R*(8 + 6*beta + beta**2)*eta_i*(-5 + 3*get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))) - Z1*(15 + 8*beta + beta**2)*(2*(-2 + get_pw_eta(eta_i, 2)) + beta*(-1 + get_pw_eta(eta_i, 2))) + Z2*(15 + 8*beta + beta**2)*(2*(-2 + get_pw_eta(eta_i, 2)) + beta*(-1 + get_pw_eta(eta_i, 2)))))/((1 + alpha)*(2 + beta)*(3 + beta)*(4 + beta)*(5 + beta)) + ((Z1 + Z2)*(-3 + get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))*xi_i)/((2 + alpha)*(1 + beta)*(3 + beta)) + (((Z1 - Z2)*(5 + 6*beta + beta**2)*eta_i*(2*(-2 + get_pw_eta(eta_i, 2)) + beta*(-1 + get_pw_eta(eta_i, 2))) + c**2*m*R*(8 + 6*beta + beta**2)*(-5 + get_pw_eta(eta_i, 4) + beta*(-1 + get_pw_eta(eta_i, 4))))*get_pw_xi(xi_i, 2))/((3 + alpha)*(1 + beta)*(2 + beta)*(4 + beta)*(5 + beta)) - ((Z1 + Z2)*(-3 + get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))*get_pw_xi(xi_i, 3))/((4 + alpha)*(1 + beta)*(3 + beta)) - (c**2*m*R*(-3 + get_pw_eta(eta_i, 2) + beta*(-1 + get_pw_eta(eta_i, 2)))*get_pw_xi(xi_i, 4))/((5 + alpha)*(1 + beta)*(3 + beta)))
   end function fun_c22two

   subroutine int_C22two(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C22two integral for the H2+ molecule.
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
      val_min_min = val_min_min + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22two(i1, i2, Z1, Z2, m, C, R, alpha, beta)
  val_max_max = val_max_max + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22two(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta)
  val_min_max = val_min_max + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22two(i1, i2+1, Z1, Z2, m, C, R, alpha, beta)
  val_max_min = val_max_min + prod_xi(i1, j1)*prod_eta(i2, j2)*fun_c22two(i1+1, i2, Z1, Z2, m, C, R, alpha, beta)
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

   end subroutine int_C22two

   function fun_c11three(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta, chi, delta) result(c11three)
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
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta, chi, delta
      type(mp_real) :: c11three

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      if (beta == 0 .and. delta == 0 .and. alpha == 0 .and. chi == 0) then
         c11three = zero ! doesn't affect result
      else if (alpha == 0 .and. chi == 0) then ! No idea why this case break everything
         c11three = c*mppi()*R**2*delta*get_pw_eta(eta_i, beta + delta)*(one/(beta + delta) - get_pw_eta(eta_i, 2)/(2 + beta + delta))*get_pw_xi(xi_i, 2)
      else if (beta == 0 .and. delta == 0) then
         c11three = c*mppi()*R**2*get_pw_eta(eta_i, 2)*chi*(-one*(get_pw_xi(xi_i, alpha + chi)/(alpha + chi)) + get_pw_xi(xi_i, 2 + alpha + chi)/(2 + alpha + chi))
      else
         c11three = c*(-2*mppi()*R**2*get_pw_eta(eta_i, beta + delta)*get_pw_xi(xi_i, alpha + chi)*(((beta + delta)*get_pw_eta(eta_i, 2)*chi)/(alpha + chi) - (get_pw_xi(xi_i, 2)*(-one*(delta**2*(get_pw_eta(eta_i, 2) - one)) + one*beta*get_pw_eta(eta_i, 2)*chi + delta*(2 + beta - beta*get_pw_eta(eta_i, 2) + get_pw_eta(eta_i, 2)*chi)))/(2 + alpha + chi)))/(one*(beta + delta)*(2 + beta + delta))
      end if
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

      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max, prefactor
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
               alpha = size(b_i_xi, 2) - j1
               do j2 = 1, size(b_i_eta, 2) ! Loop over the order of eta i
                  beta = size(b_i_eta, 2) - j2
                  do j3 = 1, size(b_j_xi, 2) ! Loop over the order of xi j
                     chi = size(b_j_xi, 2) - j3
                     do j4 = 1, size(b_j_eta, 2) ! Loop over the order of eta j
                        delta = size(b_j_eta, 2) - j4
                        prefactor = b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)
                        val_min_min = val_min_min + prefactor*fun_c11three(i1, i2, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_max = val_max_max + prefactor*fun_c11three(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_min_max = val_min_max + prefactor*fun_c11three(i1, i2+1, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_min = val_max_min + prefactor*fun_c11three(i1+1, i2, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                     end do
                  end do
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

   end subroutine int_C11three

   function fun_c22three(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta, chi, delta) result(c22three)
      !> @brief This subroutine calculates the C22three function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta, chi, delta
      type(mp_real) :: c22three

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      if (delta == 0 .and. chi == 0) then
         c22three = zero ! doesn't affect result
      else if (chi == 0) then
         c22three = -c*2*mppi()*R**2*delta*get_pw_eta(eta_i, beta + delta)*(one/(beta + delta) - (2*get_pw_eta(eta_i, 2))/(2 + beta + delta) + get_pw_eta(eta_i, 4)/(4 + beta + delta))*(-one*(get_pw_xi(xi_i, 2 + alpha)/(2 + alpha)) + get_pw_xi(xi_i, 4 + alpha)/(4 + alpha))
      else if (delta == 0) then
         c22three = c*2*mppi()*R**2*(-one*(get_pw_eta(eta_i, 2 + beta)/(2 + beta)) + get_pw_eta(eta_i, 4 + beta)/(4 + beta))*get_pw_xi(xi_i, alpha + chi)*chi*(one/(alpha + chi) - (2*get_pw_xi(xi_i, 2))/(2 + alpha + chi) + get_pw_xi(xi_i, 4)/(4 + alpha + chi))
      else
         c22three = c*2*mppi()*R**2*((-get_pw_eta(eta_i, 2 + beta + delta)/(2 + beta + delta) + get_pw_eta(eta_i, 4 + beta + delta)/(4 + beta + delta))*get_pw_xi(xi_i, alpha + chi)*chi*(one/(alpha + chi) - 2*get_pw_xi(xi_i, 2)/(2 + alpha + chi) + get_pw_xi(xi_i, 4)/(4 + alpha + chi)) - delta*get_pw_eta(eta_i, beta + delta)*(one/(beta + delta) - 2*get_pw_eta(eta_i, 2)/(2 + beta + delta) + get_pw_eta(eta_i, 4)/(4 + beta + delta))*(get_pw_xi(xi_i, 4 + alpha + chi)/(4 + alpha + chi) - get_pw_xi(xi_i, 2 + alpha + chi)/(2 + alpha + chi)))
      end if
   end function fun_c22three

   subroutine int_C22three(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C22three integral for the H2+ molecule.
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

      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max, prefactor
      integer :: i1, i2, j1, j2, j3, j4, alpha, beta, chi, delta

      zero = '0.0d0'

      ! Calculate the value of the C22 integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min = zero
            val_max_max = zero
            val_min_max = zero
            val_max_min = zero
            do j1 = 1, size(b_i_xi, 2) ! Loop over the order of xi i
               alpha = size(b_i_xi, 2) - j1
               do j2 = 1, size(b_i_eta, 2) ! Loop over the order of eta i
                  beta = size(b_i_eta, 2) - j2
                  do j3 = 1, size(b_j_xi, 2) ! Loop over the order of xi j
                     chi = size(b_j_xi, 2) - j3
                     do j4 = 1, size(b_j_eta, 2) ! Loop over the order of eta j
                        delta = size(b_j_eta, 2) - j4
                        prefactor = b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)
                        val_min_min = val_min_min + prefactor*fun_c22three(i1, i2, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_max = val_max_max + prefactor*fun_c22three(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_min_max = val_min_max + prefactor*fun_c22three(i1, i2+1, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_min = val_max_min + prefactor*fun_c22three(i1+1, i2, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                     end do
                  end do
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do
   end subroutine int_C22three

   function fun_c12three(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta, chi, delta) result(c12three)
      !> @brief This subroutine calculates the C12three function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta, chi, delta
      type(mp_real) :: c12three

      type(mp_real) :: zero, one
      zero = '0.0d0'
      one = '1.0d0'

      c12three = (2*c*mppi()*R**2*get_pw_eta(eta_i, 1 + beta + delta)*get_pw_xi(xi_i, 1 + alpha + chi)*((3 + alpha + chi)*(delta*(3 + beta + delta - (1 + beta + delta)*get_pw_eta(eta_i, 2)) + (1 + beta + delta)*get_pw_eta(eta_i, 2)*(-2 + chi) - (3 + beta + delta)*chi) + get_pw_xi(xi_i, 2)*(1 + alpha + chi)*(2*(3 + beta + delta) - delta*(3 + beta + delta) + delta*(1 + beta + delta)*get_pw_eta(eta_i, 2) + (3 + beta + delta)*chi - (1 + beta + delta)*get_pw_eta(eta_i, 2)*chi)))/((1 + beta + delta)*(3 + beta + delta)*(1 + alpha + chi)*(3 + alpha + chi))
   end function fun_c12three

   subroutine int_C12three(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C12three integral for the H2+ molecule.
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

      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max, prefactor
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
               alpha = size(b_i_xi, 2) - j1
               do j2 = 1, size(b_i_eta, 2) ! Loop over the order of eta i
                  beta = size(b_i_eta, 2) - j2
                  do j3 = 1, size(b_j_xi, 2) ! Loop over the order of xi j
                     chi = size(b_j_xi, 2) - j3
                     do j4 = 1, size(b_j_eta, 2) ! Loop over the order of eta j
                        delta = size(b_j_eta, 2) - j4
                        prefactor = b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)
                        val_min_min = val_min_min + prefactor*fun_c12three(i1, i2, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_max = val_max_max + prefactor*fun_c12three(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_min_max = val_min_max + prefactor*fun_c12three(i1, i2+1, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_min = val_max_min + prefactor*fun_c12three(i1+1, i2, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                     end do
                  end do
               end do
            end do
            result = result + val_max_max + val_min_min - val_max_min - val_min_max
         end do
      end do

   end subroutine int_C12three

   function fun_c21three(xi_i, eta_i, Z1, Z2, m, C, R, alpha, beta, chi, delta) result(c21three)
      !> @brief This subroutine calculates the C21three function for the H2+ molecule.
      !> @param xi : real : the xi coordinate
      !> @param eta : real : the eta coordinate
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param m : real : the mass of the electron
      !> @param C : real : the speed of light
      !> @param R : real : the distance between the two nuclei

      !> @param alpha : real : the alpha parameter
      !> @param beta : real : the beta parameter
      type(mp_real) :: m, C, R, Z1, Z2
      integer :: xi_i, eta_i, alpha, beta, chi, delta
      type(mp_real) :: c21three

      zero = '0.0d0'
      one = '1.0d0'

      c21three = c*2*mppi()*R**2*get_pw_eta(eta_i, 1 + beta + delta)*(-one*(one/(one + beta + delta)) + get_pw_eta(eta_i, 2)/(3 + beta + delta))*get_pw_xi(xi_i, 1 + alpha + chi)*(delta - chi)*(-one*(one/(one + alpha + chi)) + get_pw_xi(xi_i, 2)/(3 + alpha + chi))
   end function fun_c21three

   subroutine int_C21three(b_i_xi, b_i_eta, b_j_xi, b_j_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
      !> @brief This subroutine calculates the C21three integral for the H2+ molecule.
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

      type(mp_real) :: val_max_max, val_min_min, val_max_min, val_min_max, prefactor
      integer :: i1, i2, j1, j2, j3, j4, alpha, beta, chi, delta

      zero = '0.0d0'

      ! Calculate the value of the C21(3) integral at each knot and take the difference
      result = zero
      do i1 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of xi
         do i2 = 1, size(b_i_xi, 1) - 1 ! Loop over the polynomials of eta
            val_min_min = zero
            val_max_max = zero
            val_min_max = zero
            val_max_min = zero
            do j1 = 1, size(b_i_xi, 2) ! Loop over the order of xi i
               alpha = size(b_i_xi, 2) - j1
               do j2 = 1, size(b_i_eta, 2) ! Loop over the order of eta i
                  beta = size(b_i_eta, 2) - j2
                  do j3 = 1, size(b_j_xi, 2) ! Loop over the order of xi j
                     chi = size(b_j_xi, 2) - j3
                     do j4 = 1, size(b_j_eta, 2) ! Loop over the order of eta j
                        delta = size(b_j_eta, 2) - j4
                        prefactor = b_i_xi(i1, j1)*b_i_eta(i2, j2)*b_j_xi(i1, j3)*b_j_eta(i2, j4)
                        val_min_min = val_min_min + prefactor*fun_c21three(i1, i2, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_max = val_max_max + prefactor*fun_c21three(i1+1, i2+1, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_min_max = val_min_max + prefactor*fun_c21three(i1, i2+1, Z1, Z2, m, C, R, alpha, beta, chi, delta)
                        val_max_min = val_max_min + prefactor*fun_c21three(i1+1, i2, Z1, Z2, m, C, R, alpha, beta, chi, delta)
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

      double precision :: tm0, tm1

      type(mp_real), dimension(:), allocatable :: knotxi, knoteta, knotxi_eps, knoteta_eps, knotxi_tmp
      type(mp_real), dimension(:, :), allocatable :: S11one, S22one, C11one, C11two, C22one, C22two, C11three, C22three,C12three, C21three, C_mat, S_mat, vect
      type(mp_real), dimension(:), allocatable :: w, fv1, fv2
      type(mp_real), dimension(:, :, :), allocatable :: bspline_xi, bspline_eta
      logical :: debug_bool = .false.

      integer :: ntot, i, j, ierr
      integer, dimension(2) :: i2, j2

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

      print *, "Calculating the S11 integral..."
      tm0 = second()
      ! Calculate the S11 integral
      allocate (S11one(n**2, n**2))

      S11one(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, S11one, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = i, min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_s11one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, S11one(i, j))
            if (i /= j) then
               S11one(j, i) = S11one(i, j)
            end if
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

      tm1 = second()
      print *, "Time taken to calculate S11one integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the S22 integral..."
      tm0 = second()
      ! Calculate the S22 integral
      allocate (S22one(n**2, n**2))

      S22one(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, S22one, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = i, min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_s22one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, S22one(i, j))
            if (i /= j) then
               S22one(j, i) = S22one(i, j)
            end if
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

      tm1 = second()
      print *, "Time taken to calculate S22one integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the C11one integral..."
      tm0 = second()
      ! Calculate the C11one integral
      allocate (C11one(n**2, n**2))

      C11one(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C11one, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = i, min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_C11one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, C11one(i, j))
            if (i /= j) then
               C11one(j, i) = C11one(i, j)
            end if
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

      tm1 = second()
      print *, "Time taken to calculate C11one integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the C11two integral..."
      tm0 = second()
      ! Calculate the C11two integral
      allocate (C11two(n**2, n**2))

      C11two(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C11two, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = i, min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_C11two(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi, knoteta, Z1, Z2, m, C, R, C11two(i, j))
            if (i /= j) then
               C11two(j, i) = C11two(i, j)
            end if
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

      tm1 = second()
      print *, "Time taken to calculate C11two integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the C22one integral..."
      tm0 = second()
      ! Calculate the C22one integral
      allocate (C22one(n**2, n**2))

      C22one(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C22one, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = i, min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_C22one(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C22one(i, j))
            if (i /= j) then
               C22one(j, i) = C22one(i, j)
            end if
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

      tm1 = second()
      print *, "Time taken to calculate C22one integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the C22two integral..."
      tm0 = second()
      ! Calculate the C22two integral
      allocate (C22two(n**2, n**2))

      C22two(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C22two, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = i, min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_C22two(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                            bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                            knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C22two(i, j))
            if (i /= j) then
               C22two(j, i) = C22two(i, j)
            end if
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

      tm1 = second()
      print *, "Time taken to calculate C22two integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the C11three integral..."
      tm0 = second()
      ! Calculate the C11three integral
      allocate (C11three(n**2, n**2))

      C11three(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C11three, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = max(1, i-d*n), min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_C11three(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                              bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                              knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C11three(i, j))
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

      tm1 = second()
      print *, "Time taken to calculate C11three integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the C22three integral..."
      tm0 = second()
      ! Calculate the C22three integral
      allocate (C22three(n**2, n**2))

      C22three(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C22three, n, d, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = max(1, i-d*n), min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_C22three(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                              bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                              knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C22three(i, j))
         end do
      end do
      !$OMP END PARALLEL DO

      if (save_step) then
         print *, "Saving C22three matrix to file..."
         open (unit=9, file='C22three.csv', status='replace')
         do i = 1, n**2
            call write_csv(C22three(i, :), 9, 30, 10)
         end do
         close (9)
      end if

      tm1 = second()
      print *, "Time taken to calculate C22three integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the C12three integral..."
      tm0 = second()
      ! Calculate the C12three integral
      allocate (C12three(n**2, n**2))

      C12three(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C12three, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = max(1, i-d*n), min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_C12three(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                              bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                              knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C12three(i, j))
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

      tm1 = second()
      print *, "Time taken to calculate C12three integral: ", tm1 - tm0, " seconds"

      print *, "Calculating the C21three integral..."
      tm0 = second()
      ! Calculate the C21three integral
      allocate (C21three(n**2, n**2))

      C21three(:, :) = zero

      !$OMP PARALLEL DO PRIVATE(i, j, i2, j2) SHARED(bspline_xi, bspline_eta, C21three, n, d, knotxi, knoteta, Z1, Z2, m, C, R)
      do i = 1, n**2 ! Loop over the number of B-splines
         do j = max(1, i-d*n), min(i+d*n, n**2) ! Band of width kn
            i2 = indexToPair(i, n)
            j2 = indexToPair(j, n)

            call int_C21three(bspline_xi(i2(1), :, :), bspline_eta(i2(2), :, :), &
                              bspline_xi(j2(1), :, :), bspline_eta(j2(2), :, :), &
                              knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C21three(i, j))
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

      tm1 = second()
      print *, "Time taken to calculate C21three integral: ", tm1 - tm0, " seconds"

      print *, "Generating the C Matrix..."
      tm0 = second()
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

      if (save_step) then
         print *, "Saving C_mat matrix to file..."
         open (unit=13, file='C_mat.csv', status='replace')
         do i = 1, 4*n**2
            call write_csv(C_mat(i, :), 13, 50, 30)
         end do
         close (13)
      end if

      tm1 = second()
      print *, "Time taken to generate C matrix: ", tm1 - tm0, " seconds"

      print *, "Generating the S Matrix..."
      tm0 = second()
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

      if (save_step) then
         print *, "Saving S_mat matrix to file..."
         open (unit=14, file='S_mat.csv', status='replace')
         do i = 1, 4*n**2
            call write_csv(S_mat(i, :), 14, 50, 30)
         end do
         close (14)
      end if

      tm1 = second()
      print *, "Time taken to generate S matrix: ", tm1 - tm0, " seconds"
      print *, "C matrix and S matrix generated successfully."

      if (debug_bool) then
         print *, "Check if C is hermitian..."
         do i = 1, 4*n**2
            do j = 1, 4*n**2
               if (C_mat(i, j) - C_mat(j, i) > epsilonn(one)) then
                  print *, "C is not hermitian at (", i, ",", j, ") :"
                  call mpwrite(6, 35, 15, C_mat(i, j) - C_mat(j, i))
               end if
            end do
         end do

         print *, "Check if S is hermitian..."
         do i = 1, 4*n**2
            do j = 1, 4*n**2
               if (S_mat(i, j) - S_mat(j, i) > epsilonn(one)) then
                  print *, "S is not hermitian at (", i, ",", j, ") :"
                  call mpwrite(6, 35, 15, S_mat(i, j) - S_mat(j, i))
               end if
            end do
         end do
      end if

      print *, "Calculating the eigenvalues..."
      tm0 = second()
      ! Calculate the eigenvalues and eigenvectors
      allocate (w(4*n**2), fv1(4*n**2), fv2(4*n**2), vect(4*n**2, 4*n**2))

      call rsg(4*n**2, 4*n**2, C_mat, S_mat, w, 0, vect, fv1, fv2, ierr)

      print *, "Error code: ", ierr

      tm1 = second()
      print *, "Time taken to calculate eigenvalues: ", tm1 - tm0, " seconds"

      print *, "Saving logs to file..."
      ! Save logs
      open (unit=1, file='log_file', status='replace')
      write (1, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d, " and Number of BSplines to remove: ", n_remove
      write (1, '(a)') "Speed of light: "
      call mpwrite(1, 35, 15, C)
      write (1, '(a)') "Mass of the electron: "
      call mpwrite(1, 35, 15, m)
      write (1, '(a)') "Slope for eta: "
      call mpwrite(1, 35, 15, eta_slp)
      write (1, '(a)') "Non-adjusted Xi knot vector: "
      call write_lists(knotxi, 1, 35, 15)
      write (1, '(a)') "Non-adjusted Eta knot vector: "
      call write_lists(knoteta, 1, 35, 15)
      write (1, '(a)') "Internuclear distance: "
      call mpwrite(1, 35, 15, R)
      write (1, '(a)') "B-Spline basis function Order: "
      do i = 1, n**2
         i2 = indexToPair(i, n)
         write (1, '(i4, a, i4, a, i4)', advance='no') i, " ", i2(1), " ", i2(2)
         write (1, '(a)') " "
      end do
      write (1, '(a)') "Adjusted Xi knot vector: "
      call write_lists(knotxi_eps, 1, 35, 15)
      write (1, '(a)') "Adjusted Eta knot vector: "
      call write_lists(knoteta_eps, 1, 35, 15)
      if (debug_bool) then
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "S11 integral: "
         do i = 1, n**2
            call write_lists(S11one(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "S22 integral: "
         do i = 1, n**2
            call write_lists(S22one(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C11one integral: "
         do i = 1, n**2
            call write_lists(C11one(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C11two integral: "
         do i = 1, n**2
            call write_lists(C11two(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C22one integral: "
         do i = 1, n**2
            call write_lists(C22one(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C22two integral: "
         do i = 1, n**2
            call write_lists(C22two(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C11three integral: "
         do i = 1, n**2
            call write_lists(C11three(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C22three integral: "
         do i = 1, n**2
            call write_lists(C22three(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C12three integral: "
         do i = 1, n**2
            call write_lists(C12three(i, :), 1, 35, 15)
         end do
         write (1, '(a)') "--------------------------------------------------------------"
         write (1, '(a)') "C21three integral: "
         do i = 1, n**2
            call write_lists(C21three(i, :), 1, 35, 15)
         end do
      end if
      write (1, '(a)') "--------------------------------------------------------------"
      write (1, '(a)') "Value of mc²:"
      call mpwrite(1, 35, 15, m*c*c)
      write (1, '(a)') "--------------------------------------------------------------"
      write (1, '(a)') "Eigenvalues: "
      do i = 1, 4*n**2
         write (1, '(i4, a, i4)', advance='no') i, " "
         call mpwrite(1, 35, 15, w(4*n**2 - i + 1))
         write (1, '(i4, a, i4)', advance='no') i, " Translated by -mc^2 : "
         call mpwrite(1, 35, 15, w(4*n**2 - i + 1) - m*c*c)
         write (1, '(a)') " "
      end do
      close (1)

      print *, "Logs saved to log_file."
      print *, "Saving eigenvalues to file..."

      open (unit=12, file='eigenvalues.txt', status='replace')
      write (12, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d, " and Number of BSplines to remove: ", n_remove
      write (12, '(a)') "Speed of light: "
      call mpwrite(12, 35, 15, C)
      write (12, '(a)') "Mass of the electron: "
      call mpwrite(12, 35, 15, m)
      write (12, '(a)') "Slope for eta: "
      call mpwrite(12, 35, 15, eta_slp)
      write (12, '(a)') "Internuclear distance: "
      call mpwrite(12, 35, 15, R)
      write (12, '(a)') "Non-adjusted Xi knot vector: "
      call write_lists(knotxi, 12, 35, 15)
      write (12, '(a)') "Non-adjusted Eta knot vector: "
      call write_lists(knoteta, 12, 35, 15)
      write (12, '(a)') "Epsilon: "
      call mpwrite(12, 35, 15, epsilon)
      do i = 1, 4*n**2
         if (w(4*n**2 - i + 1) > -m*c*c .and. w(4*n**2 - i + 1) < m*c*c) then
            write (12, '(i4, a, i4)', advance='no') i, " "
            call mpwrite(12, 35, 15, w(4*n**2 - i + 1))
            write (12, '(i4, a, i4)', advance='no') i, " Translated by -mc^2 : "
            call mpwrite(12, 35, 15, w(4*n**2 - i + 1) - m*c*c)
            write (12, '(a)') " "
         end if
      end do
      close (12)

      print *, "Eigenvalues saved to eigenvalues.txt."

   end subroutine init_h2plus

   function epsilonn(alpha)
      !> @brief Calculate the machine epsilon
      !> @param alpha The value to calculate the machine epsilon
      USE mpmodule
      implicit type(mp_real) (a - h, o - z)

      ten = '10.d0'
      epsilonn = ten**(-mpipl)

      return
   end function epsilonn

end module h2plus
