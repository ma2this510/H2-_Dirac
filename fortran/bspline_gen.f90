module bspline_gen
   use mpmodule
   use tools_mp
   implicit none

   type(mp_real), save :: one, zero

   private

   public :: fusion_coef
   public :: print_table
   public :: init_bspine
   public :: knot_xi, knot_eta, knot_eta_lin
   public :: gen_new_knots

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
      integer :: ntot, i, itot, nhalf, nexp, nlin

      zero = '0.d0'
      one = '1.d0'
      ntot = n + d + 2*n_remove

      nhalf = (n - d + 2*n_remove) / 2
      nexp = 3 * nhalf / 4
      nlin = nhalf - nexp

      allocate (result(ntot))
      itot = 0

      do i = 1, d
         itot = itot + 1
         result(itot) = -1*one
      end do

      do i = 1, nexp ! Exponential regime - negative side
         itot = itot + 1
         result(itot) = -one * (eta_slp)**(one * (i) / nexp)
      end do

      do i = 1, nlin ! Linear regime - negative side
         itot = itot + 1
         result(itot) = -eta_slp + (i) * eta_slp / nlin
      end do

      do i = nlin, 1, -1 ! Linear regime - positive side
         itot = itot + 1
         result(itot) = eta_slp - (i) * eta_slp / nlin
      end do

      do i = nexp, 1, -1 ! Exponential regime - positive side
         itot = itot + 1
         result(itot) = one * (eta_slp)**(one * (i) / nexp)
      end do

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

   subroutine gen_new_knots(d, n, n_remove, rmin, rmax, rslp, Ncircle, Ntheta, ximin, ximax, xislp, eta_cut, nexp, knot_xi_vec_tmp, knot_xi_vec, knot_eta_vec)
      !> @brief This subroutine generates the knot vectors for xi and eta using two different regimes (concentric crcle and exponential)
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of usable B-splines
      !> @param n_remove : integer : the number of knots to remove from each end
      !> @param rmin : mp_real : the minimum radius for the circular regime
      !> @param rmax : mp_real : the maximum radius for the circular regime
      !> @param rslp : mp_real : the parameter for the generation of the radius in the circular regime
      !> @param Ncircle : integer : the number of usable B-splines in the circular regime PER NUCLEUS
      !> @param Ntheta : integer : the number of knots on each circle
      !> @param ximin : mp_real : the minimum position of the B-spline on xi-axis
      !> @param ximax : mp_real : the maximum position of the B-spline on xi-axis
      !> @param xislp : mp_real : the parameter for the generation of the knot vector on xi
      !> @param eta_cut : mp_real : the parameter for the generation of the knot vector on eta
      !> @param nexp : integer : number of exponential points in eta direction
      !> @param knot_xi_vec : mp_real(:) : the knot vector of the B-spline on xi (output)
      !> @param knot_eta_vec : mp_real(:) : the knot vector of the B-spline on eta (output)
      integer, intent(in) :: d, n, n_remove, Ncircle, Ntheta, nexp
      type(mp_real), intent(in) :: rmin, rmax, rslp, ximin, ximax, xislp, eta_cut
      type(mp_real), dimension(:), intent(inout) :: knot_xi_vec_tmp, knot_xi_vec, knot_eta_vec

      type(mp_real) :: zero, one, theta
      integer :: Nnormal, Nr, i, j
      type(mp_real), dimension(:), allocatable :: rlist, thetalist, x_pts, z_pts_minus, z_pts_plus, xi_circle, eta_circle_minus, eta_circle_plus, xi_normal, eta_normal, xi_tot, eta_tot

      zero = '0.d0'
      one = '1.d0'

      print *, "Generating new knot vectors..."
      Nr = Ncircle / Ntheta ! Number of radius points
      Nnormal = n - d + 2 * n_remove - 2*Ncircle + 2 ! Number of points in normal regime
      print *, "Number of normal points: ", Nnormal
      print *, "Number of circles: ", Nr

      ! Generate radius list
      allocate (rlist(Nr))
      do i = 1, Nr
         rlist(i) = rmin + (rmax - rmin) * (exp(rslp * (i - 1) / (Nr - 1)) - one) / (exp(rslp) - one)
      end do

      ! Generate theta list
      print *, "Generating theta points..."
      allocate (thetalist(Ntheta))
      do j = 1, Ntheta
         thetalist(j) = mppi() * j / (Ntheta + 1)
      end do

      ! Generate (x,z) points from (r,theta)
      print *, "Generating circle points..."
      allocate (x_pts(Ncircle), z_pts_plus(Ncircle), z_pts_minus(Ncircle))
      do i = 1, Nr
         do j = 1, Ntheta
            if (i==1) then
               theta = mppi() * (j - 1) / (Ntheta - 1)
            else 
               theta = thetalist(j)
            end if
            x_pts( (i - 1)*Ntheta + j ) = rlist(i) * sin(theta)
            z_pts_plus( (i - 1)*Ntheta + j ) = one + rlist(i) * cos(theta)
            z_pts_minus( (i - 1)*Ntheta + j ) = -one + rlist(i) * cos(theta)
         end do
      end do

      ! Generate circle knot vector in prolate
      print *, "Generating circle knot vectors..."
      allocate (xi_circle(Ncircle), eta_circle_minus(Ncircle), eta_circle_plus(Ncircle))
      do i = 1, Ncircle
         xi_circle(i) = (sqrt( x_pts(i)**2 + (z_pts_minus(i) + one)**2 ) + sqrt( x_pts(i)**2 + (z_pts_minus(i) - one)**2 )) / (2*one)
         eta_circle_minus(i) = (sqrt( x_pts(i)**2 + (z_pts_minus(i) + one)**2 ) - sqrt( x_pts(i)**2 + (z_pts_minus(i) - one)**2 )) / (2*one)
         xi_circle(i) = (sqrt( x_pts(i)**2 + (z_pts_plus(i) + one)**2 ) + sqrt( x_pts(i)**2 + (z_pts_plus(i) - one)**2 )) / (2*one)
         eta_circle_plus(i) = (sqrt( x_pts(i)**2 + (z_pts_plus(i) + one)**2 ) - sqrt( x_pts(i)**2 + (z_pts_plus(i) - one)**2 )) / (2*one)
      end do

      ! Generate normal knot vector in prolate
      print *, "Generating normal knot vectors..."
      allocate (xi_normal(Nnormal + Ncircle + 1), eta_normal(Nnormal))
      do i = 1, Nnormal + Ncircle + 1
         xi_normal(i) = ximin + (ximax - ximin) * (exp(xislp * (i - one) / (Nnormal + Ncircle)) - one) / (exp(xislp) - one)
      end do

      do i = 1, Nnormal
         if (i <= nexp / 2) then
            eta_normal(i) = -one * (eta_cut)**(one * (i) / (nexp / 2)) ! i from 1 to nexp/2
         else if (i > Nnormal - nexp / 2) then
            eta_normal(i) = one * (eta_cut)**(one * (Nnormal - i + one) / (nexp / 2)) ! i from Nnormal - nexp/2 + 1 to Nnormal
         else
            eta_normal(i) = -eta_cut + (2*eta_cut) * (i - nexp/2) / (Nnormal - 2*(nexp/2) + one) ! i from nexp/2 + 1 to Nnormal - nexp/2
         end if
      end do
      
      ! Combine circle and normal knot vectors
      print *, "Combining knot vectors..."
      allocate (xi_tot(2*Ncircle + Nnormal + 1))
      allocate (eta_tot(2*Ncircle + Nnormal))

      do i = 1, Ncircle
         xi_tot(i) = xi_circle(i)
         
         eta_tot(i) = eta_circle_minus(i)
         eta_tot(i + Ncircle) = eta_circle_plus(i)
      end do

      do i = 1, Nnormal + Ncircle + 1
         xi_tot(i + Ncircle) = xi_normal(i)
      end do

      do i = 1, Nnormal
         eta_tot(i + 2*Ncircle) = eta_normal(i)
      end do

      ! Sort the combined knot vectors
      print *, "Sorting knot vectors..."
      call sort_mp_real(xi_tot)
      call sort_mp_real(eta_tot)


      ! Finally fill the output knot vectors
      print *, "Filling output knot vectors..."
      knot_xi_vec_tmp = zero
      knot_xi_vec = zero
      knot_eta_vec = zero

      ! First d knots
      do i = 1, d - 1
         knot_xi_vec_tmp(i) = xi_tot(1)
         knot_xi_vec(i) = xi_tot(1)
         knot_eta_vec(i) = eta_tot(1)
      end do

      ! Middle knots
      do i = 1, n - d + 2*n_remove + 2
         knot_xi_vec_tmp(i + d - 1) = xi_tot(i)
         knot_xi_vec(i + d - 1) = xi_tot(i)
         knot_eta_vec(i + d - 1) = eta_tot(i)
      end do

      ! Last d knots
      do i = 1, d - 1
         knot_xi_vec_tmp(i + n + n_remove + 1) = xi_tot(size(xi_tot))
         knot_xi_vec(i + n + n_remove + 1) = xi_tot(size(xi_tot))
         knot_eta_vec(i + n + n_remove + 1) = eta_tot(size(eta_tot))
      end do

      knot_xi_vec_tmp(size(knot_xi_vec_tmp)) = xi_tot(size(xi_tot))

      ! Save knot vectors to files
      !open(unit=10, file = 'tmp_ex/knot_xi.txt', status='replace', action='write')
      !do i = 1, size(knot_xi_vec)
      !   call mpwrite(10, 35, 15, knot_xi_vec(i))
      !end do
      !close(10)

      !open(unit=11, file = 'tmp_ex/knot_eta.txt', status='replace', action='write')
      !do i = 1, size(knot_eta_vec)
      !   call mpwrite(11, 35, 15, knot_eta_vec(i))
      !end do
      !close(11)

      !print *, "Knot vectors saved to knot_xi.txt and knot_eta.txt"

   end subroutine gen_new_knots
end module bspline_gen
