module TwoD_Piecewise
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none
   private

   type, public :: TwoDpiecewise
      integer, public :: d
      type(mp_real), dimension(:), allocatable, public :: knot1, knot2, knot3, knotr
      type(mp_real), dimension(:, :, :), allocatable, public :: coef1, coef2, coef3, coefr
   contains
      procedure, public :: init => init_2Dpiecewise
   end type TwoDpiecewise
contains

   subroutine init_2Dpiecewise(this, d, n, Z, Zmax, Zmin, Rmin, Rmax)
      class(TwoDpiecewise), intent(inout) :: this
      integer, intent(in) :: d, n
      type(mp_real), intent(in) :: Z, Zmax, Zmin, Rmin, Rmax

      integer :: i

      allocate (this%knot1(2*n + d), this%knot2(n + d), this%knot3(n + d), this%knotr(n + d))

      this%knot1 = knot_reg1(d, n + d, Z, Zmin)
      this%knot2 = knot_reg2(d, n + d, Z, Zmin, Zmax)
      this%knot3 = knot_reg3(d, n + d, Z, Zmin, Zmax)
      this%knotr = knot_r(d, n + d, Rmin, Rmax)

      allocate (this%coef1(2*n, size(this%knot1), d), &
                this%coef2(n, size(this%knot2), d), &
                this%coef3(n, size(this%knot3), d), &
                this%coefr(n, size(this%knotr), d))

      do i = 1, 2*n
         call init_bspine(d, i, this%knot1, this%coef1(i, :, :), .false.)
      end do

      do i = 1, n
         call init_bspine(d, i, this%knot2, this%coef2(i, :, :), .false.)
      end do

      do i = 1, n
         call init_bspine(d, i, this%knot3, this%coef3(i, :, :), .false.)
      end do

      do i = 1, n
         call init_bspine(d, i, this%knotr, this%coefr(i, :, :), .false.)
      end do

   end subroutine init_2Dpiecewise

end module TwoD_Piecewise
