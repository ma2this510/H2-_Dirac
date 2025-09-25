module partial_diag
    use mpmodule
    implicit none

contains

    subroutine pdiag(n, C_mat, S_mat, eig, maxit)
        integer, intent(in) :: n, maxit
        type(mp_real), intent(inout) :: eig
        type(mp_real), dimension(n, n), intent(in) :: C_mat, S_mat
        type(mp_real), dimension(n*(n+1)/2) :: lin_C, lin_S
        type(mp_real), dimension(3*n+1) :: wa
        type(mp_real), dimension(n) :: v
        type(mp_real) :: eps

        integer :: i, j, idx, ijob

        ! Initialize parameters
        ijob = 0
        eps = '0.0d-28'
        v(1) = '1.0d0'
        do i = 2, n
            v(i) = '0.0d0'
        end do

        ! Convert C_mat and S_mat to linear storage format
        idx = 1
        do i = 1, n
            do j = 1, i
                lin_C(idx) = C_mat(i, j)
                lin_S(idx) = S_mat(i, j)
                idx = idx + 1
            end do
        end do

        ! Call the partial diagonalization routine
        call invsg(lin_C, lin_S, n, eig, v, eps, ijob, maxit, wa)

    end subroutine pdiag

end module partial_diag