module h2plus
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none

   type(mp_real), save :: one, zero

   private

   public :: init_h2plus

contains

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

   subroutine Sll(bl1, bl2, br1, br2, knot1, knot2, Z, Zmax, R, result)
        !> @brief This subroutine calculates Sll matrix element for the H2+ molecule.
        !> @param bl1 : real(:,:) : the B-spline coefficients for the left side and for the region 1
        !> @param bl2 : real(:,:) : the B-spline coefficients for the left side and for the region 2 (et 3)
        !> @param br1 : real(:,:) : the B-spline coefficients for the right side and for the region 1
        !> @param br2 : real(:,:) : the B-spline coefficients for the right side and for the region 2 (et 3)
        !> @param knot1 : real(:) : the knot vector for the region 1
        !> @param knot2 : real(:) : the knot vector for the region 2
        !> @param Z : real : the position of the nucleus
        !> @param Zmax : real : the maximum position of the B-spline on z-axis
        !> @param R : real : the position of the B-spline on r-axis
        !> @param result : real(:,:) : the result of the Sll matrix element
        type(mp_real), intent(in) :: Z, Zmax, R
        type(mp_real), dimension(:,:), intent(in) :: bl1, bl2, br1, br2
        type(mp_real), dimension(:), intent(in) :: knot1, knot2
        type(mp_real), intent(out) :: result

        type(mp_real), dimension(size(bl1, 1), size(bl2, 2)) :: term1, term2
        type(mp_real) :: result1, result2

        call integral(bl2, size(bl2, 1) -1, br2, size(br2, 1), knot2, result2)
        call integral(bl1, size(bl1, 1) -1, br1, size(br1, 1), knot1, result1)

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

      integer :: i
      type(mp_real), dimension(:), allocatable :: knot1, knot2, knot3, knotr
      type(mp_real), dimension(:, :, :), allocatable :: bspline_1, bspline_2, bspline_3, bspline_r

      zero = '0.d0'
      one = '1.d0'

      ! Generate the knot vector
      print *, "Generate the knot vector"
      allocate (knot1(2*n + d), knot2(n + d), knot3(n + d), knotr(n + d))

      knot1 = knot_reg1(d, n + d, Z, Zmin)
      knot2 = knot_reg2(d, n + d, Z, Zmin, Zmax)
      knot3 = knot_reg3(d, n + d, Z, Zmin, Zmax)
      knotr = knot_r(d, n + d, Rmin, Rmax)

      ! Generate the B-splines
      print *, "Generate the B-splines"

      allocate (bspline_1(2*n, size(knot1), d), &
                bspline_2(n, size(knot2), d), &
                bspline_3(n, size(knot3), d), &
                bspline_r(n, size(knotr), d))

      do i = 1, 2*n
         call init_bspine(d, i, knot1, bspline_1(i, :, :), .false.)
      end do

      do i = 1, n
         call init_bspine(d, i, knot2, bspline_2(i, :, :), .false.)
      end do

      do i = 1, n
         call init_bspine(d, i, knot3, bspline_3(i, :, :), .false.)
      end do

      do i = 1, n
         call init_bspine(d, i, knotr, bspline_r(i, :, :), .false.)
      end do

      print *, "B-splines generated"
   end subroutine init_h2plus

end module h2plus
