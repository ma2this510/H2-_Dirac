module tools_mp
   use mpmodule
   implicit none

   private

   public :: write_lists
   public :: write_csv
   public :: progress_bar
   public :: multiply_elem
   public :: indexToPair

   interface write_lists
      module procedure write_lists_real
      module procedure write_lists_cmplx
   end interface

contains

   subroutine write_lists_real(r1s, i1, i2, i3)
      !> @brief Write a list of mp_real
      !> @param r1s : mp_real(:) : the list of mp_real
      !> @param i1 : integer : the unit number
      !> @param i2 : integer : the field width
      !> @param i3 : integer : the number of decimal
      implicit none
      integer, intent(in) :: i1, i2, i3
      type(mp_real), intent(in), dimension(:) :: r1s

      character(i2) :: str_tmp(1)
      integer :: i_tmp

      do i_tmp = 1, size(r1s)
         str_tmp = ""
         call mpeform(r1s(i_tmp), i2, i3, str_tmp)
         write (i1, '(a)', advance='no') str_tmp
      end do
      write (i1, *)  ! End the line after the list is printed

   end subroutine write_lists_real

   subroutine write_lists_cmplx(r1s, i1, i2, i3)
      !> @brief Write a list of mp_complex
      !> @param r1s : mp_complex(:) : the list of mp_complex
      !> @param i1 : integer : the unit number
      !> @param i2 : integer : the field width
      !> @param i3 : integer : the number of decimal
      implicit none
      integer, intent(in) :: i1, i2, i3
      type(mp_complex), intent(in), dimension(:) :: r1s

      character(i2) :: str_tmp_1(1), str_tmp_2(1)
      integer :: i_tmp

      do i_tmp = 1, size(r1s)
         str_tmp_1 = ""
         str_tmp_2 = ""
         call mpeform(mpreal(r1s(i_tmp)), i2, i3, str_tmp_1)
         call mpeform(aimag(r1s(i_tmp)), i2, i3, str_tmp_2)
         write (i1, '(a1,a,a1,a,a1)', advance='no') '(', str_tmp_1, ', ', str_tmp_2, ')'
      end do
      write (i1, *)  ! End the line after the list is printed

   end subroutine write_lists_cmplx

   subroutine write_csv(r1s, i1, i2, i3)
      !> @brief Write a list of mp_real with a comma separator
      !> @param r1s : mp_real(:) : the list of mp_real
      !> @param i1 : integer : the unit number
      !> @param i2 : integer : the field width
      !> @param i3 : integer : the number of decimal
      implicit none
      integer, intent(in) :: i1, i2, i3
      type(mp_real), intent(in), dimension(:) :: r1s

      character(i2) :: str_tmp(1)
      integer :: i_tmp

      do i_tmp = 1, size(r1s) - 1
         str_tmp = ""
         call mpeform(r1s(i_tmp), i2, i3, str_tmp)
         write (i1, '(a)', advance='no') str_tmp
         write (i1, '(a)', advance='no') ','
      end do
      call mpeform(r1s(size(r1s)), i2, i3, str_tmp)
      write (i1, '(a)', advance='no') str_tmp
      write (i1, *)  ! End the line after the list is printed

   end subroutine write_csv

   subroutine progress_bar(iteration, maximum)
!
! Prints progress bar.
! Maciej Żok, 2010 MIT License
! https://github.com/macie/fortran-libs
!
! Args:
!     iteration - iteration number
!     maximum - total iterations
!
      implicit none
      integer :: iteration, maximum
      integer :: counter
      integer :: step, done

      step = nint(iteration*100/(1.0*maximum))
      done = floor(step/10.0)  ! mark every 10%

      do counter = 1, 36                    ! clear whole line - 36 chars
         write (6, '(a)', advance='no') '\b'  ! (\b - backslash)
      end do

      write (6, '(a)', advance='no') ' -> In progress... ['
      if (done .LE. 0) then
         do counter = 1, 10
            write (6, '(a)', advance='no') '='
         end do
      else if ((done .GT. 0) .and. (done .LT. 10)) then
         do counter = 1, done
            write (6, '(a)', advance='no') '#'
         end do
         do counter = done + 1, 10
            write (6, '(a)', advance='no') '='
         end do
      else
         do counter = 1, 10
            write (6, '(a)', advance='no') '#'
         end do
      end if
      write (6, '(a)', advance='no') '] '
      write (6, '(I3.1)', advance='no') step
      write (6, '(a)', advance='no') '%'
   end

   function multiply_elem(scal, mat) result(result)
      !> @brief Multiply a scalar with a matrix
      !> @param scal : mp_real : the scalar
      !> @param mat : mp_real(:,:) : the matrix
      !> @return result : mp_real(:,:) : the result of the multiplication
      implicit none
      type(mp_real), intent(in) :: scal
      type(mp_real), intent(in), dimension(:, :) :: mat
      type(mp_real), dimension(size(mat, 1), size(mat, 2)) :: result

      integer :: i_tmp, j_tmp

      do i_tmp = 1, size(mat, 1)
         do j_tmp = 1, size(mat, 2)
            result(i_tmp, j_tmp) = scal*mat(i_tmp, j_tmp)
         end do
      end do

   end function multiply_elem

   function indexToPair(index, n) result(pair)
      !> @brief Convert a linear index to a pair of indices
      !> @param index : integer : the linear index
      !> @param n : integer : the size of the second dimension
      !> @return pair : integer(2) : the pair of indices
      implicit none
      integer, intent(in) :: index, n
      integer, dimension(2) :: pair

      pair(1) = (index - 1)/n + 1
      pair(2) = mod(index - 1, n) + 1

   end function indexToPair

end module tools_mp
