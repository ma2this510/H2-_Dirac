module bspline_gen
    use mpmodule
    implicit none

    type(mp_real), save :: one, zero

    private 

    public :: fusion_coef
    public :: print_table
    public :: init_bspine
    public :: knot_reg1, knot_reg2, knot_reg3, knot_r


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

      integer :: i_tmp, j_tmp

      if (present(file)) then
         do i_tmp = 1, size(knot) - 1
            write (file, *) "------------------------------------------------------------------"
            write (file, *) "for node", i_tmp, knot(i_tmp), "<= x < ", knot(i_tmp + 1)
            do j_tmp = 1, d
               write (file, *) "x^", d - j_tmp, " : ", table(i_tmp, j_tmp)
            end do
         end do
      else
         do i_tmp = 1, size(knot) - 1
            print *, "------------------------------------------------------------------"
            print *, "for node", i_tmp, knot(i_tmp), "<= x < ", knot(i_tmp + 1)
            do j_tmp = 1, d
               print *, "x^", d - j_tmp, " : ", table(i_tmp, j_tmp)
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

   subroutine init_bspine(d, i, knot, result, display)
      !> @brief Main function to calculate the B-spline coefficients.
      !> @warning The degree and the index are the Mathematica values + 1
      !> @param d : integer : the degree of the B-spline
      !> @param i : integer : the index of the B-spline
      !> @param knot : mp_real(:) : the knot vector
      !> @param result : mp_real(:,:) : the final coef of the B-spline
      !> @param display : logical : display the result
      implicit none
      integer, intent(in) :: d, i
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(inout), dimension(size(knot), d) :: result
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

      result = calcul_double(tot, index, size(knot))

      if (display) then
         call print_table(d, knot, result)
      end if

   end subroutine init_bspine

   function knot_reg1(d, n, Z, Zmin) result(result)
      !> @brief This function generates the knot vector for the B-spline. Region 1 : Z-axis from 0 (Zmin) to Z
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of knots
      !> @param Z : real : the position of the nucleus
      !> @param Zmin : real : the minimum position of the B-spline on z-axis
      !> @return result : real(:) : the knot vector of the B-spline
      type(mp_real), intent(in) :: Z, Zmin
      integer, intent(in) :: d, n
      type(mp_real), dimension(:), allocatable :: result, tmp1, tmp2

      integer :: i
      zero = '0.d0'
      one = '1.d0'

      allocate (tmp1(n), tmp2(n))

      do i = 1, d
         tmp1(i) = -Z
      end do
      do i = d + 1, n
         tmp1(i) = -Z + Zmin*(Z/Zmin)**(((i - (d + 1))*one)/((n - d - 1)*one))
      end do

      do i = 1, d
         tmp2(i) = Z
      end do
      do i = d + 1, n
         tmp2(i) = Z - Zmin*(Z/Zmin)**(((i - (d + 1))*one)/((n - d - 1)*one))
      end do

      result = [tmp1, tmp2(n:1:-1)]

   end function knot_reg1

   function knot_reg2(d, n, Z, Zmin, Zmax) result(result)
      !> @brief This function generates the knot vector for the B-spline. Region 2 : Z-axis from Z to Zmax
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of knots
      !> @param Z : real : the position of the nucleus
      !> @param Zmin : real : the minimum position of the B-spline on z-axis
      !> @param Zmax : real : the maximum position of the B-spline on z-axis
      !> @return result : real(:) : the knot vector of the B-spline
      type(mp_real), intent(in) :: Z, Zmin, Zmax
      integer, intent(in) :: d, n
      type(mp_real), dimension(:), allocatable :: result

      integer :: i

      allocate (result(n))

      do i = 1, d
         result(i) = Z
      end do
      do i = d + 1, n - d + 1
         result(i) = Z + Zmin*((Zmax - Z)/Zmin)**(((i - (d + 1))*one)/((n - 2*d)*one))
      end do
      do i = n - d + 2, n
         result(i) = Zmax
      end do
   end function knot_reg2

   function knot_reg3(d, n, Z, Zmin, Zmax) result(result)
      !> @brief This function generates the knot vector for the B-spline. Region 3 : R-axis from 0 (Rmin) to Rmax
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of knots
      !> @param Z : real : the position of the nucleus
      !> @param Zmin : real : the minimum position of the B-spline on z-axis
      !> @param Zmax : real : the maximum position of the B-spline on z-axis
      !> @return result : real(:) : the knot vector of the B-spline
      type(mp_real), intent(in) :: Z, Zmin, Zmax
      integer, intent(in) :: d, n
      type(mp_real), dimension(:), allocatable :: result, tmp

      integer :: i

      allocate (result(n), tmp(n))

      tmp = knot_reg2(d, n, Z, Zmin, Zmax)
      do i = 1, n
         result(n - i + 1) = -tmp(i)
      end do
   end function knot_reg3

   function knot_r(d, n, Rmin, Rmax) result(result)
      !> @brief This function generates the knot vector for the B-spline. Region 4 : R-axis from Rmin to Rmax on r-axis
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of knots
      !> @param Rmin : real : the minimum position of the B-spline on r-axis
      !> @param Rmax : real : the maximum position of the B-spline on r-axis
      !> @return result : real(:) : the knot vector of the B-spline
      integer, intent(in) :: d, n
      type(mp_real), intent(in) :: Rmin, Rmax
      type(mp_real), dimension(:), allocatable :: result
      integer :: i

      zero = '0.d0'

      allocate (result(n))

      do i = 1, d
         result(i) = zero
      end do
      do i = d + 1, n - d + 1
         result(i) = Rmin*((Rmax - Rmin)/Rmin)**(((i - (d + 1))*one)/((n - 2*d)*one))
      end do
      do i = n - d + 2, n
         result(i) = Rmax
      end do
   end function knot_r

end module bspline_gen