module bspline_gen
   use mpmodule
   implicit none

   type(mp_real), save :: one, zero

   private

   public :: fusion_coef
   public :: print_table
   public :: init_bspine
   public :: knot_xi, knot_eta, knot_eta_lin

contains

   function fusion_coef(coef1, coef2)
      !> @brief This function calculates the product of two polynoms
      !>
      !> @param coef1 : real(:) : the coef of the first polynom by increasing order
      !> @param coef2 : real(:) : the coef of the second polynom by increasing order
      !> @return fusion_coef : real(:) : the coef of the product of the two polynoms by increasing order
      type(mp_real), intent(in) :: coef1(:), coef2(:)
      type(mp_real) :: fusion_coef(size(coef1) + size(coef2) - 1)

      integer :: s1, s2, i, j

      fusion_coef = zero

      s1 = size(coef1)
      s2 = size(coef2)

      do i = 1, s1
         do j = 1, s2
            fusion_coef(i + j - 1) = fusion_coef(i + j - 1) + coef1(i)*coef2(j)
         end do
      end do

   end function fusion_coef

   recursive subroutine rec_coef(d, i, knot, table, sol_int, tot, index)
      !> @brief Recursive function to calculate the coefficients of the B-spline
      !> @param d : integer : the degree of the B-spline
      !> @param i : integer : the index of the B-spline
      !> @param knot : real(:) : the knot vector
      !> @param table : real(:,:) : the coef of the B-spline
      !> @param sol_int : real(:) : Variable to keep track of the filling of tot.
      !> @param tot : real(:) : Variable TOTAL fill when reach end of recursion
      !> @param index : integer : the index of the current B-spline

      integer, intent(in) :: d, i
      type(mp_real), intent(in), dimension(:) :: knot
      type(mp_real), intent(in), dimension(:) :: table
      type(mp_real), intent(out), dimension(:, :) :: tot
      integer, intent(out) :: sol_int
      integer, intent(out), dimension(:) :: index

      type(mp_real) :: coef1(2), coef2(2)
      type(mp_real) :: denum1, denum2

      type(mp_real), dimension(:), allocatable :: res1, res2

      denum1 = knot(i + d - 1) - knot(i)
      denum2 = knot(i + d) - knot(i + 1)

      if (denum1 == zero) then
         coef1(1) = zero
         coef1(2) = zero
      else
         coef1(1) = 1/denum1
         coef1(2) = -knot(i)/denum1
      end if

      if (denum2 == zero) then
         coef2(1) = zero
         coef2(2) = zero
      else
         coef2(1) = -1/denum2
         coef2(2) = knot(i + d)/denum2
      end if

      allocate (res1(size(table) + 1), res2(size(table) + 1))

      res1 = fusion_coef(table, coef1)
      res2 = fusion_coef(table, coef2)

      if (d > 2) then
         call rec_coef(d - 1, i, knot, res1, sol_int, tot, index)
         call rec_coef(d - 1, i + 1, knot, res2, sol_int, tot, index)
      else
         sol_int = sol_int + 1
         index(sol_int) = i
         tot(sol_int, :) = res1

         sol_int = sol_int + 1
         index(sol_int) = i + 1
         tot(sol_int, :) = res2
      end if

   end subroutine rec_coef

   subroutine print_table(d, knot, table, file)
      !> @brief Print the polynome between each node of a given spline
      !> @param d : integer : the degree of the B-spline
      !> @param knot : real(:) : the knot vector
      !> @param table : real(:,:) : the coef of the B-spline
      !> @param file : integer : the file to write the result
      implicit none
      integer, intent(in) :: d
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(in), dimension(size(knot), d) :: table
      integer, intent(in), optional :: file

      character(25) :: str_tmp(1), str_tmp_1(1)
      integer :: i_tmp, j_tmp

      if (present(file)) then
         do i_tmp = 1, size(knot) - 1
            write (file, *) "------------------------------------------------------------------"
            call mpeform(knot(i_tmp), 25, 10, str_tmp)
            call mpeform(knot(i_tmp + 1), 25, 10, str_tmp_1)
            write (file, *) "for node", i_tmp, str_tmp, "<= x < ", str_tmp_1
            do j_tmp = 1, d
               call mpeform(table(i_tmp, j_tmp), 25, 10, str_tmp)
               write (file, *) "x^", d - j_tmp, " : ", str_tmp
            end do
         end do
      else
         do i_tmp = 1, size(knot) - 1
            print *, "------------------------------------------------------------------"
            call mpeform(knot(i_tmp), 25, 10, str_tmp)
            call mpeform(knot(i_tmp + 1), 25, 10, str_tmp_1)
            print *, "for node", i_tmp, str_tmp, "<= x < ", str_tmp_1
            do j_tmp = 1, d
               call mpeform(table(i_tmp, j_tmp), 25, 10, str_tmp)
               print *, "x^", d - j_tmp, " : ", str_tmp
            end do
         end do
      end if
   end subroutine print_table

   function calcul_double(table, index, s) result(total)
      !> @brief Sum the different branch linked to the same internal node and return the final coefs per node
      !> @param table : real(:,:) : the un-proccessed coef of the B-spline
      !> @param index : integer(:) : the index of the current contribution
      !> @param s : integer : the number of nodes
      !> @return total : real(:,:) : the proccesed coef of the B-spline
      implicit none
      type(mp_real), intent(in) :: table(:, :)
      integer, intent(in) :: index(:)
      integer, intent(in) :: s
      type(mp_real), allocatable :: total(:, :)

      integer :: i, j

      allocate (total(s, size(table, 2)))

      total = zero

      do i = 1, size(table, 1)
         do j = 1, size(table, 2)
            total(index(i), j) = total(index(i), j) + table(i, j)
         end do
      end do

   end function calcul_double

   subroutine init_bspine(d, i, knot, ntot, result, display)
      !> @brief Main function to calculate the B-spline coefficients.
      !> @warning The degree and the index are the Mathematica values + 1
      !> @param d : integer : the degree of the B-spline
      !> @param i : integer : the index of the B-spline
      !> @param knot : mp_real(:) : the knot vector
      !> @param ntot : integer : the number of usable knots
      !> @param result : mp_real(:,:) : the final coef of the B-spline
      !> @param display : logical : display the result
      implicit none
      integer, intent(in) :: d, i, ntot
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(inout), dimension(:, :) :: result
      logical, intent(in) :: display

      type(mp_real), dimension(1) :: table
      integer :: sol_int
      type(mp_real), dimension(2**(d - 1), d) :: tot
      integer, dimension(2**(d - 1)) :: index

      zero = '0.d0'
      one = '1.d0'

      if (d == 1) then
         print *, "The degree of the B-spline must be greater than 1"
         stop
      else if (i + d > size(knot)) then
         print *, "For the given degree and index, the knot vector needs to be at least of size ", i + d
         stop
      end if

      table = zero
      table(1) = one
      sol_int = 0
      tot = zero
      index = 0

      call rec_coef(d, i, knot, table, sol_int, tot, index)

      result = calcul_double(tot, index, ntot)

      if (display) then
         call print_table(d, knot, result)
      end if

   end subroutine init_bspine

   function knot_xi(d, n, n_remove, ximin, ximax, xi_slp) result(result)
      !> @brief This function generates the knot vector for xi.
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of usable B-splines
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param ximin : real : the minimum position of the B-spline on xi-axis
      !> @param ximax : real : the maximum position of the B-spline on xi-axis
      !> @param xi_slp : real : the parameter for the generation of the knot vector on xi
      !> @return result : real(:) : the knot vector of the B-spline
      type(mp_real), intent(in) :: ximin, ximax, xi_slp
      integer, intent(in) :: d, n, n_remove
      type(mp_real), dimension(:), allocatable :: result

      type(mp_real) :: one, zero

      integer :: ntot, i, itot

      zero = '0.d0'
      one = '1.d0'

      ntot = n + d + 2*n_remove

      allocate (result(ntot))

      itot = 0

      do i = 1, d
         itot = itot + 1
         result(itot) = ximin
      end do

      do i = 1, n - d + 2*n_remove
         itot = itot + 1
         result(itot) = ximin + (ximax - ximin) * (exp(xi_slp * i / (n - d + 2*n_remove + 1)) - 1) / (exp(xi_slp) - 1)
      end do

      do i = 1, d
         itot = itot + 1
         result(itot) = ximax
      end do

   end function knot_xi

   function knot_eta(d, n, n_remove, eta_slp) result(result)
      !> @brief This function generates the knot vector for eta.
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of usable B-splines
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param eta_slp : real : the parameter for the generation of the knot vector on eta
      !> @return result : real(:) : the knot vector of the B-spline
      integer, intent(in) :: d, n, n_remove
      type(mp_real), intent(in) :: eta_slp
      type(mp_real), dimension(:), allocatable :: result

      type(mp_real) :: one, zero
      integer :: ntot, i, itot

      zero = '0.d0'
      one = '1.d0'
      ntot = n + d + 2*n_remove

      allocate (result(ntot))
      itot = 0

      do i = 1, d
         itot = itot + 1
         result(itot) = -1*one
      end do

      if (mod(n - d + 2*n_remove, 4) == 0) then
         do i = 1, (n - d + 2*n_remove)/4 ! Exponential regime - negative side
            itot = itot + 1
            result(itot) = -one * (eta_slp)**(4 * one * (i) / (n - d + 2*n_remove))
         end do

         do i = 1, (n - d + 2*n_remove)/4 ! Linear regime - negative side
            itot = itot + 1
            result(itot) = -eta_slp + 4 * (i) * eta_slp / (n - d + 2*n_remove)
         end do

         do i = (n - d + 2*n_remove)/4, 1, -1 ! Linear regime - positive side
            itot = itot + 1
            result(itot) = eta_slp - 4 * (i) * eta_slp / (n - d + 2*n_remove)
         end do

         do i = (n - d + 2*n_remove)/4, 1, -1 ! Exponential regime - positive side
            itot = itot + 1
            result(itot) = one * (eta_slp)**(4 * one * (i) / (n - d + 2*n_remove))
         end do

      else

         do i = 1, (n - d + 2*n_remove)/4 + 1 ! Exponential regime - negative side
            itot = itot + 1
            result(itot) = -one * (eta_slp)**((i * one) / ((n - d + 2*n_remove)/4 + one))
         end do

         do i = 1, (n - d + 2*n_remove)/4 ! Linear regime - negative side
            itot = itot + 1
            result(itot) = -eta_slp + (i) * eta_slp / ((n - d + 2*n_remove)/4)
         end do

         do i = (n - d + 2*n_remove)/4, 1, -1 ! Linear regime - positive side
            itot = itot + 1
            result(itot) = eta_slp - (i) * eta_slp / ((n - d + 2*n_remove)/4)
         end do

         do i = (n - d + 2*n_remove)/4 + 1, 1, -1 ! Exponential regime - positive side
            itot = itot + 1
            result(itot) = one * (eta_slp)**((i * one) / ((n - d + 2*n_remove)/4 + one))
         end do

      end if

      do i = 1, d
         itot = itot + 1
         result(itot) = one
      end do

   end function knot_eta

   function knot_eta_lin(d, n, n_remove, eta_slp) result(result)
      !> @brief This function generates the knot vector for eta.
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of usable B-splines
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param eta_slp : real : the parameter for the generation of the knot vector on eta
      !> @return result : real(:) : the knot vector of the B-spline
      integer, intent(in) :: d, n, n_remove
      type(mp_real), intent(in) :: eta_slp
      type(mp_real), dimension(:), allocatable :: result

      type(mp_real) :: one, zero
      integer :: ntot, i, itot

      zero = '0.d0'
      one = '1.d0'
      ntot = n + d + 2*n_remove

      allocate (result(ntot))
      itot = 0

      do i = 1, d
         itot = itot + 1
         result(itot) = -1*one
      end do

      do i = 1, n - d + 2*n_remove
         itot = itot + 1
         result(itot) = -1*one + 2*one*i / ( (n - d + 2*n_remove + 1)*one )
      end do

      do i = 1, d
         itot = itot + 1
         result(itot) = one
      end do

   end function knot_eta_lin
end module bspline_gen
