module h2plus
   use mpmodule
   use bspline_gen
   use tools_mp
   use TwoD_Piecewise
   implicit none

   type(mp_real), save :: one, zero

   private

   public :: init_h2plus

contains
   subroutine integral(b1, order1, b2, order2, knot, result)
      !> @brief Calculate the integral of the product of two Piecewise polynomial of any order
      !> @param b1 : real(:,:) : the coef of the first Piecewise polynomial
      !> @param order1 : integer : the max order of the first Piecewise polynomial
      !> @param b2 : real(:,:) : the coef of the second Piecewise polynomial
      !> @param order2 : integer : the max order of the second Piecewise polynomial
      !> @param knot : real(:) : the knot vector
      !> @param result : real : the result of the integral
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      integer, intent(in) :: order1, order2
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(out) :: result

      type(mp_real), dimension(size(b1, 1)) :: int_per_knot

      integer :: order_prod, i_tmp, j_tmp
      type(mp_real), dimension(:, :), allocatable :: produit, int_init, int_final

      zero = '0.d0'
      one = '1.d0'

      allocate (produit(size(b1, 1), 2*size(b1, 2) - 1))
      allocate (int_final(size(b1, 1), 2*size(b1, 2)))
      allocate (int_init(size(b1, 1), 2*size(b1, 2)))

      produit = zero
      order_prod = order1 + order2
      int_final = zero
      int_init = zero
      int_per_knot = zero

      do i_tmp = 1, size(b1, 1) ! Loop over the number of Piecewise polynomial
         produit(i_tmp, :) = fusion_coef(b1(i_tmp, :), b2(i_tmp, :))
      end do

      do i_tmp = 1, size(produit, 1) - 1 ! Loop over the number of Piecewise polynomial
         do j_tmp = 1, size(produit, 2) ! Loop over the different orders
            if (order_prod - j_tmp + 2 /= 0) then
               if (knot(i_tmp + 1) /= zero) then
                int_final(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*(knot(i_tmp + 1)**(order_prod - j_tmp + 2))/(order_prod - j_tmp + 2)
                  ! +2 because of the integral that add one order
               end if
               if (knot(i_tmp) /= zero) then
                  int_init(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*(knot(i_tmp)**(order_prod - j_tmp + 2))/(order_prod - j_tmp + 2)
               end if
            else
               if (knot(i_tmp + 1) /= zero) then
                  int_final(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*log(knot(i_tmp + 1))
               end if
               if (knot(i_tmp) /= zero) then
                  int_init(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*log(knot(i_tmp))
               end if
            end if
         end do
      end do

      do i_tmp = 1, size(b1, 1) - 1 ! Loop over the number of Piecewise polynomial
         int_per_knot(i_tmp) = zero
         do j_tmp = 1, size(produit, 2) ! Loop over the different orders
            int_per_knot(i_tmp) = int_per_knot(i_tmp) + int_final(i_tmp, j_tmp) - int_init(i_tmp, j_tmp)
         end do
      end do

      result = zero
      do i_tmp = 1, size(b1, 1) - 1
         result = result + int_per_knot(i_tmp)
      end do

   end subroutine integral

   subroutine integral_r_z(b1, b2, n_remove, result)
      !> @brief This subroutine calculates the integral of the B-spline coefficients for the H2+ molecule.
      !> @param b1 : TwoDpiecewise : the B-spline coefficients for the left side
      !> @param b2 : TwoDpiecewise : the B-spline coefficients for the right side
      !> @param result : real : the result of the integral
      type(TwoDpiecewise), intent(in) :: b1, b2
      integer, intent(in) :: n_remove
      type(mp_real), dimension(:, :, :), intent(out) :: result

      integer :: i, j, k, tmp1, tmp2
      type(mp_real), dimension(:, :), allocatable :: result1, result2, result3, resultr

      allocate (result1(size(b1%coef1, 1) - 2*n_remove, size(b2%coef1, 1) - 2*n_remove))
      allocate (result2(size(b1%coef2, 1) - 2*n_remove, size(b2%coef2, 1) - 2*n_remove))
      allocate (result3(size(b1%coef3, 1) - 2*n_remove, size(b2%coef3, 1) - 2*n_remove))
      allocate (resultr(size(b1%coefr, 1) - 2*n_remove, size(b2%coefr, 1) - 2*n_remove))

      result1 = zero
      result2 = zero
      result3 = zero
      resultr = zero

      print *, size(b1%coef2, 1), size(b2%coef2, 1), n_remove

      do i = 1 + n_remove, size(b1%coef1, 1) - n_remove
         do j = 1 + n_remove, size(b2%coef1, 1) - n_remove
            call integral(b1%coef1(i, :, :), b1%d, b2%coef1(j, :, :), b2%d, b1%knot1, result1(i - n_remove, j - n_remove))
         end do
      end do
      do i = 1 + n_remove, size(b1%coef2, 1) - n_remove
         do j = 1 + n_remove, size(b2%coef2, 1) - n_remove
            call integral(b1%coef2(i, :, :), b1%d, b2%coef2(j, :, :), b2%d, b1%knot2, result2(i - n_remove, j - n_remove))
         end do
      end do
      do i = 1 + n_remove, size(b1%coef3, 1) - n_remove
         do j = 1 + n_remove, size(b2%coef3, 1) - n_remove
            call integral(b1%coef3(i, :, :), b1%d, b2%coef3(j, :, :), b2%d, b1%knot3, result3(i - n_remove, j - n_remove))
         end do
      end do
      do i = 1 + n_remove, size(b1%coefr, 1) - n_remove
         do j = 1 + n_remove, size(b2%coefr, 1) - n_remove
            call integral(b1%coefr(i, :, :), b1%d, b2%coefr(j, :, :), b2%d, b1%knotr, resultr(i - n_remove, j - n_remove))
         end do
      end do

      do i = 1, size(result, 3)
         tmp1 = 1
         tmp2 = 1
         result(tmp1:size(result2, 1), tmp2:size(result2, 2), i) = result2
         tmp1 = size(result2, 1) + 1
         tmp2 = size(result2, 2) + 1
         result(tmp1:size(result1, 1)+tmp1-1, tmp2:size(result1, 2)+tmp2-1, i) = result1
         tmp1 = size(result1, 1) + size(result2, 1) + 1
         tmp2 = size(result1, 2) + size(result2, 2) + 1
         result(tmp1:size(result3, 1)+tmp1-1, tmp2:size(result3, 2)+tmp2-1, i) = result3
      end do
   end subroutine integral_r_z

   subroutine Sll(b1, b2, Z, Zmax, R, n_remove, result)
      !> @brief This subroutine calculates Sll matrix element for the H2+ molecule.
      !> @param b1 : TwoDpiecewise : the B-spline coefficients for the left side
      !> @param b2 : TwoDpiecewise : the B-spline coefficients for the right side
      !> @param Z : real : the position of the nucleus
      !> @param Zmax : real : the maximum position of the B-spline on z-axis
      !> @param R : real : the position of the B-spline on r-axis
      !> @param result : real(:,:) : the result of the Sll matrix element
      type(mp_real), intent(in) :: Z, Zmax, R
      type(TwoDpiecewise), intent(in) :: b1, b2
      integer, intent(in) :: n_remove
      type(mp_real), dimension(:, :, :), intent(out) :: result

      call integral_r_z(b1, b2, n_remove, result)
   end subroutine Sll

   subroutine init_h2plus(d, n, n_remove, Z1, Z2, Z, Zmax, Zmin, Rmin, Rmax)
      !> @brief This subroutine initializes the B-spline coefficients for the H2+ molecule.
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of knots
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param Z1 : real : the number of protons for the first atom
      !> @param Z2 : real : the number of protons for the second atom
      !> @param Z : real : position of both nuclei
      !> @param Zmax : real : maximum position of the B-spline on z-axis
      !> @param Zmin : real : minimum position of the B-spline on z-axis
      !> @param Rmin : real : minimum position of the B-spline on r-axis
      !> @param Rmax : real : maximum position of the B-spline on r-axis
      type(mp_real), intent(in) :: Z1, Z2, Z, Zmax, Zmin, Rmin, Rmax
      integer, intent(in) :: d, n, n_remove

      type(TwoDpiecewise) :: b1, b2
      type(mp_real), dimension(:, :, :), allocatable :: result_Sll

      zero = '0.d0'
      one = '1.d0'

      print *, "Initializing Knot vectors and B-splines for H2+ molecule"
      call b1%init(d, n, Z, Zmax, Zmin, Rmin, Rmax)
      call b2%init(d, n, Z, Zmax, Zmin, Rmin, Rmax)

      print *, "B-splines generated"

      print *, "Calculating Sll matrix element"
      allocate (result_Sll(4*n - 6*n_remove, 4*n - 6*n_remove, n - 2*n_remove))
      call Sll(b1, b2, Z, Zmax, Rmin, n_remove, result_Sll)
      print *, "Sll matrix element calculated"

   end subroutine init_h2plus

end module h2plus
