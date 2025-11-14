module partial_diag
    use mpmodule
    implicit none

contains

    subroutine pdiag(n, lin_C, lin_S, eig, maxit, v)
        integer, intent(in) :: n, maxit
        type(mp_real), intent(inout) :: eig
        type(mp_real), dimension(n), intent(inout) :: v
        type(mp_real), dimension(n*(n+1)/2), intent(in) :: lin_C, lin_S
        type(mp_real), dimension(n*(n+1)/2) :: rep_C, rep_S
        type(mp_real), dimension(3*n+1) :: wa
        type(mp_real) :: eps

        integer :: i, ijob

        ! Initialize parameters
        ijob = 0
        eps = '0.0d-28'
        v(1) = '1.0d0'
        do i = 2, n
            v(i) = '0.0d0'
        end do
        
        ! Call the partial diagonalization routine
        do i = 1, 3
            rep_C(:) = lin_C(:)
            rep_S(:) = lin_S(:)
            call invsg(rep_C, rep_S, n, eig, v, eps, ijob, maxit, wa)
            print *, 'Iteration ', i, ', Numerical eigenvalue = '
            call mpwrite(6, 35, 15, eig)
        end do

    end subroutine pdiag

end module partial_diag