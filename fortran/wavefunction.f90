module wavefunction
    use mpmodule
    use bspline_gen
    use tools_mp
    implicit none

    type(mp_real), save :: one, zero

    private

    public :: get_wavefunc

contains
    function bspline_val_xi(xi, bspline)
        type(mp_real), intent(in) :: xi
        type(mp_real), dimension(:, :), intent(in) :: bspline

        integer :: k, l1, alpha

        type(mp_real) :: bspline_val_xi

        bspline_val_xi = zero
        do k = 1, size(bspline, 1) - 1 ! Loop over the polynomials
            do l1 = 1, size(bspline, 2) ! loop over differents orders
                alpha = size(bspline, 2) - 1

                bspline_val_xi = bspline_val_xi + bspline(k, l1) * (xi**(alpha))
            end do
        end do

        return
    end function bspline_val_xi

    function bspline_val_eta(eta, bspline)
        type(mp_real), intent(in) :: eta
        type(mp_real), dimension(:, :), intent(in) :: bspline

        integer :: k, l1, alpha

        type(mp_real) :: bspline_val_eta

        bspline_val_eta = zero
        do k = 1, size(bspline, 1) - 1 ! Loop over the polynomials
            do l1 = 1, size(bspline, 2) ! loop over differents orders
                alpha = size(bspline, 2) - 1

                bspline_val_eta = bspline_val_eta + bspline(k, l1) * (eta**(alpha))
            end do
        end do

        return
    end function bspline_val_eta

    function prefactor2(xi, eta)
        type(mp_real), intent(in) :: xi, eta
        type(mp_real) :: prefactor2

        prefactor2 = (xi**2 - 1) * (1 - eta**2)
        return
    end function prefactor2

    subroutine get_wavefunc(b_xi, b_eta, point_xi, point_eta, ximax, ximin, eigvec, folder_name)
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta
        type(mp_real), dimension(:), intent(in) :: eigvec
        type(mp_real), dimension(:), intent(in) :: point_xi, point_eta
        type(mp_real), intent(in) :: ximax, ximin
        character(len=20), intent(in) :: folder_name

        type(mp_real), dimension(:, :, :), allocatable :: wavefun
        integer :: n_xi, n_eta, i, j, i1, i2, i3, i4, idx1, idx2

        zero = '0.d0'
        one = '1.d0'

        n_xi = size(point_xi) 
        n_eta = size(point_eta)

        allocate(wavefun(4, n_xi, n_eta))

        do i = 1, n_xi ! Loop over first coords
            do j = 1, n_eta ! Loop over seconds coords

                wavefun(1, i, j) = zero
                wavefun(2, i, j) = zero
                wavefun(3, i, j) = zero
                wavefun(4, i, j) = zero

                do i1 = 1, size(b_xi, 1) ! Loop over first dim basis fun
                    do i2 = 1, size(b_eta, 1) ! Loop over second dim basis fun
                        do i3 = 1, size(b_xi, 1) ! Loop over first dim basis fun
                            do i4 = 1, size(b_eta, 1) ! Loop over second dim basis fun 
                                idx1 = (i1 - 1) * size(b_xi, 1) + i2
                                idx2 = (i3 - 1) * size(b_eta, 1) + i4

                                wavefun(1, i, j) = wavefun(1, i, j) + eigvec(idx1) * bspline_val_xi(point_xi(i), b_xi(i1, :, :)) * bspline_val_eta(point_eta(j), b_xi(i2, :, :)) * eigvec(idx2) * bspline_val_xi(point_xi(i), b_xi(i3, :, :)) * bspline_val_eta(point_eta(j), b_xi(i4, :, :))

                                wavefun(2, i, j) = wavefun(2, i, j) + eigvec(size(b_xi, 1) + idx1) * prefactor2(point_xi(i), point_eta(j))**2 * bspline_val_eta(point_xi(i), b_xi(i1, :, :)) * bspline_val_xi(point_eta(j), b_xi(i2, :, :)) * eigvec(size(b_xi, 1) + idx2) * bspline_val_xi(point_xi(i), b_xi(i3, :, :)) * bspline_val_eta(point_eta(j), b_xi(i4, :, :))

                                wavefun(3, i, j) = wavefun(3, i, j) + eigvec(2 * size(b_xi, 1) + idx1) * bspline_val_xi(point_xi(i), b_xi(i1, :, :)) * bspline_val_eta(point_eta(j), b_xi(i2, :, :)) * eigvec(2 * size(b_xi, 1) + idx2) * bspline_val_xi(point_xi(i), b_xi(i3, :, :)) * bspline_val_eta(point_eta(j), b_xi(i4, :, :))

                                wavefun(4, i, j) = wavefun(4, i, j) + eigvec(3 * size(b_xi, 1) + idx1) * prefactor2(point_xi(i), point_eta(j))**2 * bspline_val_xi(point_xi(i), b_xi(i1, :, :)) * bspline_val_eta(point_eta(j), b_xi(i2, :, :)) * eigvec(3 * size(b_xi, 1) + idx2) * bspline_val_xi(point_xi(i), b_xi(i3, :, :)) * bspline_val_eta(point_eta(j), b_xi(i4, :, :))
                            end do
                        end do
                    end do
                end do

            end do
        end do

        open(unit = 15, file = trim(folder_name) // '/wavefun.txt', status = 'replace')
        write(15, '(a)') "Knots in xi :"
        call write_csv(point_xi, 15, 35, 15)

        write(15, '(a)') "Knots in eta :"
        call write_csv(point_eta, 15, 35, 15)

        write(15, '(a)') "Component 1 of wavefunction :"

        do i = 1, n_xi ! Loop over first coords
            call write_csv(wavefun(1, i, :), 15, 35, 15)
        end do

        write(15, '(a)') "Component 2 of wavefunction :"

        do i = 1, n_xi ! Loop over first coords
            call write_csv(wavefun(2, i, :), 15, 35, 15)
        end do

        write(15, '(a)') "Component 3 of wavefunction :"

        do i = 1, n_xi ! Loop over first coords
            call write_csv(wavefun(3, i, :), 15, 35, 15)
        end do

        write(15, '(a)') "Component 4 of wavefunction :"

        do i = 1, n_xi ! Loop over first coords
            call write_csv(wavefun(4, i, :), 15, 35, 15)
        end do

        close(15)

    end subroutine get_wavefunc
end module wavefunction