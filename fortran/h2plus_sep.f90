module h2plus_sep
    use mpmodule
    use bspline_gen
    use tools_mp
    use partial_diag
    implicit none

    type(mp_real), save :: one, zero

    private

    public :: init_h2plus_sep

contains
    subroutine int_C11one(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the C11one integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        type(mp_real), dimension(:, :, :, :), allocatable :: prod_xi, prod_eta
        integer :: alpha, beta, i, j, k, l
        integer, dimension(2) :: i2, j2

        allocate (prod_xi(size(b_xi, 1), size(b_xi, 1), size(b_xi, 2), size(b_xi, 3) + size(b_xi, 3) - 1), prod_eta(size(b_eta, 1), size(b_eta, 1), size(b_eta, 2), size(b_eta, 3) + size(b_eta, 3) - 1))

        do i = 1, size(b_xi, 1) ! Loop over the first basis functions
            do j = 1, size(b_xi, 1) ! Loop over the second basis functions
                do k = 1, size(b_xi, 2) ! Loop over the polynomials
                    ! Calculate the product of the B-spline coefficients
                    prod_xi(i, j, k, :) = fusion_coef(b_xi(i, k, :), b_xi(j, k, :))
                    prod_eta(i, j, k, :) = fusion_coef(b_eta(i, k, :), b_eta(j, k, :))
                end do
            end do
        end do

        allocate(xi_1(size(prod_xi, 1), size(prod_xi, 2)), xi_2(size(prod_xi, 1), size(prod_xi, 2)), eta_1(size(prod_eta, 1), size(prod_eta, 2)), eta_2(size(prod_eta, 1), size(prod_eta, 2)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(prod_xi, 1) ! Loop over first basis functions
            do j = 1, size(prod_xi, 2) ! Loop over second basis functions
                do k = 1, size(prod_xi, 3) - 1! Loop over the polynomials
                    do l = 1, size(prod_xi, 4) ! loop over order
                        alpha = size(prod_xi, 4) - l
                        xi_1(i, j) = xi_1(i, j) + prod_xi(i, j, k, l) * knotxi(k + 1)**(2 + alpha) * ((c**2 * m * knotxi(k + 1)) / (3 + alpha) - (Z1 + Z2) / (2 * R + R * alpha)) ! Finale
                        xi_1(i, j) = xi_1(i, j) - prod_xi(i, j, k, l) * knotxi(k)**(2 + alpha) * ((c**2 * m * knotxi(k)) / (3 + alpha) - (Z1 + Z2) / (2 * R + R * alpha)) ! Initial

                        xi_2(i, j) = xi_2(i, j) + prod_xi(i, j, k, l) * knotxi(k + 1)**(1 + alpha) / (1 + alpha) ! Finale
                        xi_2(i, j) = xi_2(i, j) - prod_xi(i, j, k, l) * knotxi(k)**(1 + alpha) / (1 + alpha) ! Initial
                    end do
                end do

                do k = 1, size(prod_eta, 3) - 1 ! Loop over the polynomials
                    do l = 1, size(prod_eta, 4) ! loop over order
                        beta = size(prod_eta, 4) - l
                        eta_1(i, j) = eta_1(i, j) + prod_eta(i, j, k, l) * knoteta(k + 1)**(1 + beta) / (1 + beta) ! Finale
                        eta_1(i, j) = eta_1(i, j) - prod_eta(i, j, k, l) * knoteta(k)**(1 + beta) / (1 + beta) ! Initial

                        eta_2(i, j) = eta_2(i, j) - prod_eta(i, j, k, l) * (knoteta(k + 1)**(2 + beta) * ((-Z1 + Z2) * one / (2 + beta) + (c**2 * m * R * knoteta(k + 1)) / (3 + beta))) / R ! Finale
                        eta_2(i, j) = eta_2(i, j) + prod_eta(i, j, k, l) * (knoteta(k)**(2 + beta) * ((-Z1 + Z2) * one / (2 + beta) + (c**2 * m * R * knoteta(k)) / (3 + beta))) / R ! Initial
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = 2 * mppi() * (R**3) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) + xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_C11one

    subroutine int_C11two(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the C22one integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        type(mp_real), dimension(:, :, :, :), allocatable :: prod_xi, prod_eta
        integer :: alpha, beta, i, j, k, l
        integer, dimension(2) :: i2, j2

        allocate (prod_xi(size(b_xi, 1), size(b_xi, 1), size(b_xi, 2), size(b_xi, 3) + size(b_xi, 3) - 1), prod_eta(size(b_eta, 1), size(b_eta, 1), size(b_eta, 2), size(b_eta, 3) + size(b_eta, 3) - 1))

        do i = 1, size(b_xi, 1) ! Loop over the first basis functions
            do j = 1, size(b_xi, 1) ! Loop over the second basis functions
                do k = 1, size(b_xi, 2) ! Loop over the polynomials
                    ! Calculate the product of the B-spline coefficients
                    prod_xi(i, j, k, :) = fusion_coef(b_xi(i, k, :), b_xi(j, k, :))
                    prod_eta(i, j, k, :) = fusion_coef(b_eta(i, k, :), b_eta(j, k, :))
                end do
            end do
        end do

        allocate(xi_1(size(prod_xi, 1), size(prod_xi, 2)), xi_2(size(prod_xi, 1), size(prod_xi, 2)), eta_1(size(prod_eta, 1), size(prod_eta, 2)), eta_2(size(prod_eta, 1), size(prod_eta, 2)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(prod_xi, 1) ! Loop over first basis functions
            do j = 1, size(prod_xi, 2) ! Loop over second basis functions
                do k = 1, size(prod_xi, 3) - 1! Loop over the polynomials
                    do l = 1, size(prod_xi, 4) ! loop over order
                        alpha = size(prod_xi, 4) - l
                        xi_1(i, j) = xi_1(i, j) + prod_xi(i, j, k, l) * (knotxi(k + 1)**(2 + alpha) * (((Z1 + Z2) * one) / (2 + alpha) + (c**2 * m * R * knotxi(k + 1)) / (3 + alpha))) / R ! Finale
                        xi_1(i, j) = xi_1(i, j) - prod_xi(i, j, k, l) * (knotxi(k)**(2 + alpha) * (((Z1 + Z2) * one) / (2 + alpha) + (c**2 * m * R * knotxi(k)) / (3 + alpha))) / R! Initial

                        xi_2(i, j) = xi_2(i, j) + prod_xi(i, j, k, l) * knotxi(k + 1)**(1 + alpha) / (1 + alpha) ! Finale
                        xi_2(i, j) = xi_2(i, j) - prod_xi(i, j, k, l) * knotxi(k)**(1 + alpha) / (1 + alpha) ! Initial
                    end do
                end do

                do k = 1, size(prod_eta, 3) - 1 ! Loop over the polynomials
                    do l = 1, size(prod_eta, 4) ! loop over order
                        beta = size(prod_eta, 4) - l
                        eta_1(i, j) = eta_1(i, j) + prod_eta(i, j, k, l) * knoteta(k + 1)**(1 + beta) / (1 + beta) ! Finale
                        eta_1(i, j) = eta_1(i, j) - prod_eta(i, j, k, l) * knoteta(k)**(1 + beta) / (1 + beta) ! Initial

                        eta_2(i, j) = eta_2(i, j) + prod_eta(i, j, k, l) * (knoteta(k + 1)**(2 + beta) * (((Z1 - Z2) * one) / (2 + beta) + (c**2 * m * R * knoteta(k + 1)) / (3 + beta))) / R ! Finale
                        eta_2(i, j) = eta_2(i, j) - prod_eta(i, j, k, l) * (knoteta(k)**(2 + beta) * (((Z1 - Z2) * one) / (2 + beta) + (c**2 * m * R * knoteta(k)) / (3 + beta))) / R ! Initial
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = -2 * mppi() * (R**3) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) - xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_C11two

    subroutine int_C22one(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the C22one integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        type(mp_real), dimension(:, :, :, :), allocatable :: prod_xi, prod_eta
        integer :: alpha, beta, i, j, k, l
        integer, dimension(2) :: i2, j2

        allocate (prod_xi(size(b_xi, 1), size(b_xi, 1), size(b_xi, 2), size(b_xi, 3) + size(b_xi, 3) - 1), prod_eta(size(b_eta, 1), size(b_eta, 1), size(b_eta, 2), size(b_eta, 3) + size(b_eta, 3) - 1))

        do i = 1, size(b_xi, 1) ! Loop over the first basis functions
            do j = 1, size(b_xi, 1) ! Loop over the second basis functions
                do k = 1, size(b_xi, 2) ! Loop over the polynomials
                    ! Calculate the product of the B-spline coefficients
                    prod_xi(i, j, k, :) = fusion_coef(b_xi(i, k, :), b_xi(j, k, :))
                    prod_eta(i, j, k, :) = fusion_coef(b_eta(i, k, :), b_eta(j, k, :))
                end do
            end do
        end do

        allocate(xi_1(size(prod_xi, 1), size(prod_xi, 2)), xi_2(size(prod_xi, 1), size(prod_xi, 2)), eta_1(size(prod_eta, 1), size(prod_eta, 2)), eta_2(size(prod_eta, 1), size(prod_eta, 2)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(prod_xi, 1) ! Loop over first basis functions
            do j = 1, size(prod_xi, 2) ! Loop over second basis functions
                do k = 1, size(prod_xi, 3) - 1! Loop over the polynomials
                    do l = 1, size(prod_xi, 4) ! loop over order
                        alpha = size(prod_xi, 4) - l
                        xi_1(i, j) = xi_1(i, j) + prod_xi(i, j, k, l) * (knotxi(k + 1)**(2 + alpha) * (-one * (((Z1 + Z2) * ((-one / (2 + alpha)) + knotxi(k + 1)**2 / (4 + alpha))) / R) + c**2 * m * knotxi(k + 1) * ((-one / (3 + alpha)) + knotxi(k + 1)**2 / (5 + alpha)))) ! Finale
                        xi_1(i, j) = xi_1(i, j) - prod_xi(i, j, k, l) * knotxi(k)**(2 + alpha) * ((((Z1 + Z2) * ((one / (2 + alpha)) - knotxi(k)**2 / (4 + alpha))) / R) - c**2 * m * knotxi(k) * ((one / (3 + alpha)) - knotxi(k)**2 / (5 + alpha))) ! Initial

                        xi_2(i, j) = xi_2(i, j) + prod_xi(i, j, k, l) * (knotxi(k + 1)**(3 + alpha) / (3 + alpha) - (knotxi(k + 1)**(1 + alpha) / (1 + alpha))) ! Finale
                        xi_2(i, j) = xi_2(i, j) - prod_xi(i, j, k, l) * (knotxi(k)**(3 + alpha) / (3 + alpha) - (knotxi(k)**(1 + alpha) / (1 + alpha))) ! Initial
                    end do
                end do

                do k = 1, size(prod_eta, 3) - 1 ! Loop over the polynomials
                    do l = 1, size(prod_eta, 4) ! loop over order
                        beta = size(prod_eta, 4) - l
                        eta_1(i, j) = eta_1(i, j) + prod_eta(i, j, k, l) * (knoteta(k + 1)**(1 + beta) / (1 + beta) - knoteta(k + 1)**(3 + beta) / (3 + beta)) ! Finale
                        eta_1(i, j) = eta_1(i, j) - prod_eta(i, j, k, l) * (knoteta(k)**(1 + beta) / (1 + beta) - knoteta(k)**(3 + beta) / (3 + beta)) ! Initial

                        eta_2(i, j) = eta_2(i, j) + prod_eta(i, j, k, l) * (knoteta(k + 1)**(2 + beta) * ((-Z1 + Z2) * ((-one / (2 + beta)) + knoteta(k + 1)**2 / (4 + beta)) + c**2 * m * R * knoteta(k + 1) * ((-one / (3 + beta)) + knoteta(k + 1)**2 / (5 + beta)))) / R ! Finale
                        eta_2(i, j) = eta_2(i, j) - prod_eta(i, j, k, l) * (knoteta(k)**(2 + beta) * ((-Z1 + Z2) * ((-one / (2 + beta)) + knoteta(k)**2 / (4 + beta)) + c**2 * m * R * knoteta(k) * ((-one / (3 + beta)) + knoteta(k)**2 / (5 + beta)))) / R ! Initial
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = 2 * mppi() * (R**3) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) + xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_C22one

    subroutine int_C22two(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the C22one integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        type(mp_real), dimension(:, :, :, :), allocatable :: prod_xi, prod_eta
        integer :: alpha, beta, i, j, k, l
        integer, dimension(2) :: i2, j2

        allocate (prod_xi(size(b_xi, 1), size(b_xi, 1), size(b_xi, 2), size(b_xi, 3) + size(b_xi, 3) - 1), prod_eta(size(b_eta, 1), size(b_eta, 1), size(b_eta, 2), size(b_eta, 3) + size(b_eta, 3) - 1))

        do i = 1, size(b_xi, 1) ! Loop over the first basis functions
            do j = 1, size(b_xi, 1) ! Loop over the second basis functions
                do k = 1, size(b_xi, 2) ! Loop over the polynomials
                    ! Calculate the product of the B-spline coefficients
                    prod_xi(i, j, k, :) = fusion_coef(b_xi(i, k, :), b_xi(j, k, :))
                    prod_eta(i, j, k, :) = fusion_coef(b_eta(i, k, :), b_eta(j, k, :))
                end do
            end do
        end do

        allocate(xi_1(size(prod_xi, 1), size(prod_xi, 2)), xi_2(size(prod_xi, 1), size(prod_xi, 2)), eta_1(size(prod_eta, 1), size(prod_eta, 2)), eta_2(size(prod_eta, 1), size(prod_eta, 2)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(prod_xi, 1) ! Loop over first basis functions
            do j = 1, size(prod_xi, 2) ! Loop over second basis functions
                do k = 1, size(prod_xi, 3) - 1! Loop over the polynomials
                    do l = 1, size(prod_xi, 4) ! loop over order
                        alpha = size(prod_xi, 4) - l
                        xi_1(i, j) = xi_1(i, j) + prod_xi(i, j, k, l) * (knotxi(k + 1)**(2 + alpha) * ((Z1 + Z2) * ((-one / (2 + alpha)) + knotxi(k + 1)**2 / (4 + alpha)) + c**2 * m * R * knotxi(k + 1) * ((-one / (3 + alpha)) + knotxi(k + 1)**2 / (5 + alpha)))) / R ! Finale
                        xi_1(i, j) = xi_1(i, j) - prod_xi(i, j, k, l) * (knotxi(k)**(2 + alpha) * ((Z1 + Z2) * ((-one / (2 + alpha)) + knotxi(k)**2 / (4 + alpha)) + c**2 * m * R * knotxi(k) * ((-one / (3 + alpha)) + knotxi(k)**2 / (5 + alpha)))) / R ! Initial

                        xi_2(i, j) = xi_2(i, j) + prod_xi(i, j, k, l) * (knotxi(k + 1)**(3 + alpha) / (3 + alpha) - (knotxi(k + 1)**(1 + alpha) / (1 + alpha))) ! Finale
                        xi_2(i, j) = xi_2(i, j) - prod_xi(i, j, k, l) * (knotxi(k)**(3 + alpha) / (3 + alpha) - (knotxi(k)**(1 + alpha) / (1 + alpha))) ! Initial
                    end do
                end do

                do k = 1, size(prod_eta, 3) - 1 ! Loop over the polynomials
                    do l = 1, size(prod_eta, 4) ! loop over order
                        beta = size(prod_eta, 4) - l
                        eta_1(i, j) = eta_1(i, j) + prod_eta(i, j, k, l) * (knoteta(k + 1)**(1 + beta) / (1 + beta) - knoteta(k + 1)**(3 + beta) / (3 + beta)) ! Finale
                        eta_1(i, j) = eta_1(i, j) - prod_eta(i, j, k, l) * (knoteta(k)**(1 + beta) / (1 + beta) - knoteta(k)**(3 + beta) / (3 + beta)) ! Initial

                        eta_2(i, j) = eta_2(i, j) + prod_eta(i, j, k, l) * ((knoteta(k + 1)**(2 + beta) * ((Z1 - Z2) * ((-one / (2 + beta)) + knoteta(k + 1)**2 / (4 + beta)) + c**2 * m * R * knoteta(k + 1) * ((-one / (3 + beta)) + knoteta(k + 1)**2 / (5 + beta)))) / R) ! Finale
                        eta_2(i, j) = eta_2(i, j) - prod_eta(i, j, k, l) * ((knoteta(k)**(2 + beta) * ((Z1 - Z2) * ((-one / (2 + beta)) + knoteta(k)**2 / (4 + beta)) + c**2 * m * R * knoteta(k) * ((-one / (3 + beta)) + knoteta(k)**2 / (5 + beta)))) / R) ! Initial
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = -2 * mppi() * (R**3) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) + xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_C22two

    subroutine int_C11three(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the C11three integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        integer :: alpha, beta, chi, delta, i, j, k, l1, l2
        integer, dimension(2) :: i2, j2

        allocate(xi_1(size(b_xi, 1), size(b_xi, 1)), xi_2(size(b_xi, 1), size(b_xi, 1)), eta_1(size(b_eta, 1), size(b_eta, 1)), eta_2(size(b_eta, 1), size(b_eta, 1)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(b_xi, 1) ! Loop over first basis functions
            do j = 1, size(b_xi, 1) ! Loop over second basis functions
                do k = 1, size(b_xi, 2) - 1! Loop over the polynomials
                    do l1 = 1, size(b_xi, 3) ! loop over order of first basis function
                        do l2 = 1, size(b_xi, 3) ! loop over order of second basis function
                            alpha = size(b_xi, 3) - l1
                            chi = size(b_xi, 3) - l2
                            if (chi /= 0) then
                                xi_1(i, j) = xi_1(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * knotxi(k + 1)**(alpha + chi) * chi * ((-one / (alpha + chi)) + knotxi(k + 1)**2 / (2 + alpha + chi))  ! Finale
                                xi_1(i, j) = xi_1(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * knotxi(k)**(alpha + chi) * chi * ((-one / (alpha + chi)) + knotxi(k)**2 / (2 + alpha + chi))  ! Initial
                            end if

                            xi_2(i, j) = xi_2(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * knotxi(k + 1)**(2 + alpha + chi) / (2 + alpha + chi) ! Finale
                            xi_2(i, j) = xi_2(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * knotxi(k)**(2 + alpha + chi) / (2 + alpha + chi) ! Initial
                        end do
                    end do
                end do

                do k = 1, size(b_eta, 2) - 1 ! Loop over the polynomials
                    do l1 = 1, size(b_eta, 3) ! loop over order
                        do l2 = 1, size(b_eta, 3) ! loop over order
                            beta = size(b_eta, 3) - l1
                            delta = size(b_eta, 3) - l2
                            eta_1(i, j) = eta_1(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * knoteta(k + 1)**(2 + beta + delta) / (2 + beta + delta) ! Finale
                            eta_1(i, j) = eta_1(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * knoteta(k)**(2 + beta + delta) / (2 + beta + delta) ! Initial

                            if (delta /= 0) then
                                eta_2(i, j) = eta_2(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * delta * knoteta(k + 1)**(beta + delta) * (one / (beta + delta) - knoteta(k + 1)**2 / (2 + beta + delta)) ! Finale
                                eta_2(i, j) = eta_2(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * delta * knoteta(k)**(beta + delta) * (one / (beta + delta) - knoteta(k)**2 / (2 + beta + delta)) ! Initial
                            end if
                        end do
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = 2 * c * mppi() * (R**2) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) + xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_C11three

    subroutine int_C22three(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the C22three integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        integer :: alpha, beta, chi, delta, i, j, k, l1, l2
        integer, dimension(2) :: i2, j2

        allocate(xi_1(size(b_xi, 1), size(b_xi, 1)), xi_2(size(b_xi, 1), size(b_xi, 1)), eta_1(size(b_eta, 1), size(b_eta, 1)), eta_2(size(b_eta, 1), size(b_eta, 1)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(b_xi, 1) ! Loop over first basis functions
            do j = 1, size(b_xi, 1) ! Loop over second basis functions
                do k = 1, size(b_xi, 2) - 1! Loop over the polynomials
                    do l1 = 1, size(b_xi, 3) ! loop over order of first basis function
                        do l2 = 1, size(b_xi, 3) ! loop over order of second basis function
                            alpha = size(b_xi, 3) - l1
                            chi = size(b_xi, 3) - l2
                            if (chi /= 0) then
                                xi_1(i, j) = xi_1(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * knotxi(k + 1)**(alpha + chi) * chi * (one / (alpha + chi) - (2 * knotxi(k + 1)**2) / (2 + alpha + chi) + knotxi(k + 1)**4 / (4 + alpha + chi))  ! Finale
                                xi_1(i, j) = xi_1(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * knotxi(k)**(alpha + chi) * chi * (one / (alpha + chi) - (2 * knotxi(k)**2) / (2 + alpha + chi) + knotxi(k)**4 / (4 + alpha + chi))  ! Initial
                            end if

                            xi_2(i, j) = xi_2(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * (knotxi(k + 1)**(4 + alpha + chi) / (4 + alpha + chi) - (knotxi(k + 1)**(2 + alpha + chi) / (2 + alpha + chi))) ! Finale
                            xi_2(i, j) = xi_2(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * (knotxi(k)**(4 + alpha + chi) / (4 + alpha + chi) - (knotxi(k)**(2 + alpha + chi) / (2 + alpha + chi))) ! Initial
                        end do
                    end do
                end do

                do k = 1, size(b_eta, 2) - 1 ! Loop over the polynomials
                    do l1 = 1, size(b_eta, 3) ! loop over order
                        do l2 = 1, size(b_eta, 3) ! loop over order
                            beta = size(b_eta, 3) - l1
                            delta = size(b_eta, 3) - l2
                            eta_1(i, j) = eta_1(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * (knoteta(k + 1)**(2 + beta + delta) / (2 + beta + delta) - knoteta(k + 1)**(4 + beta + delta) / (4 + beta + delta)) ! Finale
                            eta_1(i, j) = eta_1(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * (knoteta(k)**(2 + beta + delta) / (2 + beta + delta) - knoteta(k)**(4 + beta + delta) / (4 + beta + delta)) ! Initial

                            if (delta /= 0) then
                                eta_2(i, j) = eta_2(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * delta * knoteta(k + 1)**(beta + delta) * (one / (beta + delta) - (2 * knoteta(k + 1)**2) / (2 + beta + delta) + knoteta(k + 1)**4 / (4 + beta + delta)) ! Finale
                                eta_2(i, j) = eta_2(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * delta * knoteta(k)**(beta + delta) * (one / (beta + delta) - (2 * knoteta(k)**2) / (2 + beta + delta) + knoteta(k)**4 / (4 + beta + delta)) ! Initial
                            end if
                        end do
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = -2 * c * mppi() * (R**2) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) + xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_C22three

    subroutine int_C12three(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the C12three integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, xi_3, xi_4, eta_1, eta_2, eta_3, eta_4
        integer :: alpha, beta, chi, delta, i, j, k, l1, l2
        integer, dimension(2) :: i2, j2

        allocate(xi_1(size(b_xi, 1), size(b_xi, 1)), xi_2(size(b_xi, 1), size(b_xi, 1)), xi_3(size(b_xi, 1), size(b_xi, 1)), xi_4(size(b_xi, 1), size(b_xi, 1)), eta_1(size(b_eta, 1), size(b_eta, 1)), eta_2(size(b_eta, 1), size(b_eta, 1)), eta_3(size(b_eta, 1), size(b_eta, 1)), eta_4(size(b_eta, 1), size(b_eta, 1)))

        xi_1 = zero
        xi_2 = zero
        xi_3 = zero
        xi_4 = zero
        eta_1 = zero
        eta_2 = zero
        eta_3 = zero
        eta_4 = zero

        do i = 1, size(b_xi, 1) ! Loop over first basis functions
            do j = 1, size(b_xi, 1) ! Loop over second basis functions
                do k = 1, size(b_xi, 2) - 1! Loop over the polynomials
                    do l1 = 1, size(b_xi, 3) ! loop over order of first basis function
                        do l2 = 1, size(b_xi, 3) ! loop over order of second basis function
                            alpha = size(b_xi, 3) - l1
                            chi = size(b_xi, 3) - l2

                            !                            if (chi /= 0) then
                            xi_1(i, j) = xi_1(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * chi * (knotxi(k + 1)**(3 + alpha + chi) / (3 + alpha + chi) - knotxi(k + 1)**(1 + alpha + chi) / (1 + alpha + chi)) ! Finale
                            xi_1(i, j) = xi_1(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * chi * (knotxi(k)**(3 + alpha + chi) / (3 + alpha + chi) - knotxi(k)**(1 + alpha + chi) / (1 + alpha + chi)) ! Initial

                            xi_2(i, j) = xi_2(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * (knotxi(k + 1)**(3 + alpha + chi) / (3 + alpha + chi) - knotxi(k + 1)**(1 + alpha + chi) / (1 + alpha + chi)) ! Finale
                            xi_2(i, j) = xi_2(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * (knotxi(k)**(3 + alpha + chi) / (3 + alpha + chi) - knotxi(k)**(1 + alpha + chi) / (1 + alpha + chi)) ! Initial

                            !                            end if

                            xi_3(i, j) = xi_3(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * knotxi(k + 1)**(3 + alpha + chi) / (3 + alpha + chi) ! Finale
                            xi_3(i, j) = xi_3(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * knotxi(k)**(3 + alpha + chi) / (3 + alpha + chi) ! Initial

                            xi_4(i, j) = xi_4(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * (knotxi(k + 1)**(3 + alpha + chi) / (3 + alpha + chi) - 2 * knotxi(k + 1)**(1 + alpha + chi) / (1 + alpha + chi)) ! Finale
                            xi_4(i, j) = xi_4(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * (knotxi(k)**(3 + alpha + chi) / (3 + alpha + chi) - 2 * knotxi(k)**(1 + alpha + chi) / (1 + alpha + chi)) ! Initial
                        end do
                    end do
                end do

                do k = 1, size(b_eta, 2) - 1 ! Loop over the polynomials
                    do l1 = 1, size(b_eta, 3) ! loop over order
                        do l2 = 1, size(b_eta, 3) ! loop over order
                            beta = size(b_eta, 3) - l1
                            delta = size(b_eta, 3) - l2

                            !                            if (delta /= 0) then
                            eta_1(i, j) = eta_1(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * (knoteta(k + 1)**(1 + beta + delta) / (1 + beta + delta) - knoteta(k + 1)**(3 + beta + delta) / (3 + beta + delta)) ! Finale
                            eta_1(i, j) = eta_1(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * (knoteta(k)**(1 + beta + delta) / (1 + beta + delta) - knoteta(k)**(3 + beta + delta) / (3 + beta + delta)) ! Initial

                            eta_2(i, j) = eta_2(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * delta * (knoteta(k + 1)**(1 + beta + delta) / (1 + beta + delta) - knoteta(k + 1)**(3 + beta + delta) / (3 + beta + delta)) ! Finale
                            eta_2(i, j) = eta_2(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * delta * (knoteta(k)**(1 + beta + delta) / (1 + beta + delta) - knoteta(k)**(3 + beta + delta) / (3 + beta + delta)) ! Initial
                            !                            end if

                            eta_3(i, j) = eta_3(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * knoteta(k + 1)**(3 + beta + delta) / (3 + beta + delta) ! Finale
                            eta_3(i, j) = eta_3(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * knoteta(k + 1)**(3 + beta + delta) / (3 + beta + delta) ! Initial

                            eta_4(i, j) = eta_4(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * knoteta(k + 1)**(3 + beta + delta) / (3 + beta + delta) ! Finale
                            eta_4(i, j) = eta_4(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * knoteta(k)**(3 + beta + delta) / (3 + beta + delta) ! Initial
                        end do
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = 2 * c * mppi() * (R**2) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) - xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)) - one * xi_3(i2(2), j2(2)) * eta_3(i2(1), j2(1)) + one * xi_4(i2(1), j2(1)) * eta_4(i2(2), j2(2)))
                ! result(i, j) = 2*c*mppi()*(R**2)*((one*xi_1(i2(1), j2(1)) + xi_4(i2(1), j2(1)))*(eta_1(i2(2), j2(2))) + &
                !                                     xi_2(i2(1), j2(1))*(eta_3(i2(2), j2(2)) - eta_5(i2(2), j2(2))) + &
                !                                     xi_1(i2(1), j2(1))*eta_2(i2(2), j2(2)) - xi_3(i2(1), j2(1))*eta_3(i2(2), j2(2)))
            end do
        end do

    end subroutine int_C12three

    subroutine int_C21three(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the C21three integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        integer :: alpha, beta, chi, delta, i, j, k, l1, l2
        integer, dimension(2) :: i2, j2

        allocate(xi_1(size(b_xi, 1), size(b_xi, 1)), xi_2(size(b_xi, 1), size(b_xi, 1)), eta_1(size(b_eta, 1), size(b_eta, 1)), eta_2(size(b_eta, 1), size(b_eta, 1)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(b_xi, 1) ! Loop over first basis functions
            do j = 1, size(b_xi, 1) ! Loop over second basis functions
                do k = 1, size(b_xi, 2) - 1! Loop over the polynomials
                    do l1 = 1, size(b_xi, 3) ! loop over order of first basis function
                        do l2 = 1, size(b_xi, 3) ! loop over order of second basis function
                            alpha = size(b_xi, 3) - l1
                            chi = size(b_xi, 3) - l2
                            !                            if (chi /= 0) then
                            xi_1(i, j) = xi_1(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * chi * (knotxi(k + 1)**(3 + alpha + chi) / (3 + alpha + chi) - knotxi(k + 1)**(1 + alpha + chi) / (1 + alpha + chi))  ! Finale
                            xi_1(i, j) = xi_1(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * chi * (knotxi(k)**(3 + alpha + chi) / (3 + alpha + chi) - knotxi(k)**(1 + alpha + chi) / (1 + alpha + chi))  ! Initial

                            xi_2(i, j) = xi_2(i, j) + b_xi(i, k, l1) * b_xi(j, k, l2) * (knotxi(k + 1)**(3 + alpha + chi) / (3 + alpha + chi) - knotxi(k + 1)**(1 + alpha + chi) / (1 + alpha + chi)) ! Finale
                            xi_2(i, j) = xi_2(i, j) - b_xi(i, k, l1) * b_xi(j, k, l2) * (knotxi(k)**(3 + alpha + chi) / (3 + alpha + chi) - knotxi(k)**(1 + alpha + chi) / (1 + alpha + chi)) ! Initial
                            !                            end if

                        end do
                    end do
                end do

                do k = 1, size(b_eta, 2) - 1 ! Loop over the polynomials
                    do l1 = 1, size(b_eta, 3) ! loop over order
                        do l2 = 1, size(b_eta, 3) ! loop over order
                            beta = size(b_eta, 3) - l1
                            delta = size(b_eta, 3) - l2

                            !                            if (delta /= 0) then
                            eta_1(i, j) = eta_1(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * (knoteta(k + 1)**(1 + beta + delta) / (1 + beta + delta) - knoteta(k + 1)**(3 + beta + delta) / (3 + beta + delta)) ! Finale
                            eta_1(i, j) = eta_1(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * (knoteta(k)**(1 + beta + delta) / (1 + beta + delta) - knoteta(k)**(3 + beta + delta) / (3 + beta + delta)) ! Initial

                            eta_2(i, j) = eta_2(i, j) + b_eta(i, k, l1) * b_eta(j, k, l2) * delta * (knoteta(k + 1)**(1 + beta + delta) / (1 + beta + delta) - knoteta(k + 1)**(3 + beta + delta) / (3 + beta + delta)) ! Finale
                            eta_2(i, j) = eta_2(i, j) - b_eta(i, k, l1) * b_eta(j, k, l2) * delta * (knoteta(k)**(1 + beta + delta) / (1 + beta + delta) - knoteta(k)**(3 + beta + delta) / (3 + beta + delta)) ! Initial
                            !                            end if
                        end do
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = 2 * c * mppi() * (R**2) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) - xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_C21three

    subroutine int_S11one(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the S11one integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S11 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        type(mp_real), dimension(:, :, :, :), allocatable :: prod_xi, prod_eta
        integer :: alpha, beta, i, j, k, l
        integer, dimension(2) :: i2, j2

        allocate (prod_xi(size(b_xi, 1), size(b_xi, 1), size(b_xi, 2), size(b_xi, 3) + size(b_xi, 3) - 1), prod_eta(size(b_eta, 1), size(b_eta, 1), size(b_eta, 2), size(b_eta, 3) + size(b_eta, 3) - 1))

        do i = 1, size(b_xi, 1) ! Loop over the first basis functions
            do j = 1, size(b_xi, 1) ! Loop over the second basis functions
                do k = 1, size(b_xi, 2) ! Loop over the polynomials
                    ! Calculate the product of the B-spline coefficients
                    prod_xi(i, j, k, :) = fusion_coef(b_xi(i, k, :), b_xi(j, k, :))
                    prod_eta(i, j, k, :) = fusion_coef(b_eta(i, k, :), b_eta(j, k, :))
                end do
            end do
        end do

        allocate(xi_1(size(prod_xi, 1), size(prod_xi, 2)), xi_2(size(prod_xi, 1), size(prod_xi, 2)), eta_1(size(prod_eta, 1), size(prod_eta, 2)), eta_2(size(prod_eta, 1), size(prod_eta, 2)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(prod_xi, 1) ! Loop over first basis functions
            do j = 1, size(prod_xi, 2) ! Loop over second basis functions
                do k = 1, size(prod_xi, 3) - 1! Loop over the polynomials
                    do l = 1, size(prod_xi, 4) ! loop over order
                        alpha = size(prod_xi, 4) - l
                        xi_1(i, j) = xi_1(i, j) + prod_xi(i, j, k, l) * knotxi(k + 1)**(3 + alpha) / (3 + alpha) ! Finale
                        xi_1(i, j) = xi_1(i, j) - prod_xi(i, j, k, l) * knotxi(k)**(3 + alpha) / (3 + alpha) ! Initial

                        xi_2(i, j) = xi_2(i, j) + prod_xi(i, j, k, l) * knotxi(k + 1)**(1 + alpha) / (1 + alpha) ! Finale
                        xi_2(i, j) = xi_2(i, j) - prod_xi(i, j, k, l) * knotxi(k)**(1 + alpha) / (1 + alpha) ! Initial
                    end do
                end do

                do k = 1, size(prod_eta, 3) - 1 ! Loop over the polynomials
                    do l = 1, size(prod_eta, 4) ! loop over order
                        beta = size(prod_eta, 4) - l
                        eta_1(i, j) = eta_1(i, j) + prod_eta(i, j, k, l) * knoteta(k + 1)**(1 + beta) / (1 + beta) ! Finale
                        eta_1(i, j) = eta_1(i, j) - prod_eta(i, j, k, l) * knoteta(k)**(1 + beta) / (1 + beta) ! Initial

                        eta_2(i, j) = eta_2(i, j) + prod_eta(i, j, k, l) * knoteta(k + 1)**(3 + beta) / (3 + beta) ! Finale
                        eta_2(i, j) = eta_2(i, j) - prod_eta(i, j, k, l) * knoteta(k)**(3 + beta) / (3 + beta) ! Initial
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = 2 * mppi() * (R**3) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) - xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_S11one

    subroutine int_S22one(b_xi, b_eta, knotxi, knoteta, Z1, Z2, m, C, R, result)
        !> @brief This subroutine calculates the S22one integral for the H2+ molecule.
        !> @param b_xi : real(:, :, :) : the B-spline coefficients for the xi direction
        !> @param b_eta : real(:, :, :) : the B-spline coefficients for the eta direction
        !> @param knotxi : real(:) : the knot vector for the xi direction
        !> @param knoteta : real(:) : the knot vector for the eta direction
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei

        !> @return S22 : real : the value of the S11 integral
        type(mp_real), intent(in) :: Z1, Z2, m, C, R
        type(mp_real), intent(out), dimension(:, :) :: result
        type(mp_real), dimension(:), intent(in) :: knotxi, knoteta
        type(mp_real), dimension(:, :, :), intent(in) :: b_xi, b_eta

        type(mp_real), dimension(:, :), allocatable :: xi_1, xi_2, eta_1, eta_2
        type(mp_real), dimension(:, :, :, :), allocatable :: prod_xi, prod_eta
        integer :: alpha, beta, i, j, k, l
        integer, dimension(2) :: i2, j2

        allocate (prod_xi(size(b_xi, 1), size(b_xi, 1), size(b_xi, 2), size(b_xi, 3) + size(b_xi, 3) - 1), prod_eta(size(b_eta, 1), size(b_eta, 1), size(b_eta, 2), size(b_eta, 3) + size(b_eta, 3) - 1))

        do i = 1, size(b_xi, 1) ! Loop over the first basis functions
            do j = 1, size(b_xi, 1) ! Loop over the second basis functions
                do k = 1, size(b_xi, 2) ! Loop over the polynomials
                    ! Calculate the product of the B-spline coefficients
                    prod_xi(i, j, k, :) = fusion_coef(b_xi(i, k, :), b_xi(j, k, :))
                    prod_eta(i, j, k, :) = fusion_coef(b_eta(i, k, :), b_eta(j, k, :))
                end do
            end do
        end do

        allocate(xi_1(size(prod_xi, 1), size(prod_xi, 2)), xi_2(size(prod_xi, 1), size(prod_xi, 2)), eta_1(size(prod_eta, 1), size(prod_eta, 2)), eta_2(size(prod_eta, 1), size(prod_eta, 2)))

        xi_1 = zero
        xi_2 = zero
        eta_1 = zero
        eta_2 = zero

        do i = 1, size(prod_xi, 1) ! Loop over first basis functions
            do j = 1, size(prod_xi, 2) ! Loop over second basis functions
                do k = 1, size(prod_xi, 3) - 1! Loop over the polynomials
                    do l = 1, size(prod_xi, 4) ! loop over order
                        alpha = size(prod_xi, 4) - l
                        xi_1(i, j) = xi_1(i, j) + prod_xi(i, j, k, l) * (knotxi(k + 1)**(5 + alpha) / (5 + alpha) - knotxi(k + 1)**(3 + alpha) / (3 + alpha)) ! Finale
                        xi_1(i, j) = xi_1(i, j) - prod_xi(i, j, k, l) * (knotxi(k)**(5 + alpha) / (5 + alpha) - knotxi(k)**(3 + alpha) / (3 + alpha)) ! Initial

                        xi_2(i, j) = xi_2(i, j) + prod_xi(i, j, k, l) * (knotxi(k + 1)**(3 + alpha) / (3 + alpha) - knotxi(k + 1)**(1 + alpha) / (1 + alpha)) ! Finale
                        xi_2(i, j) = xi_2(i, j) - prod_xi(i, j, k, l) * (knotxi(k)**(3 + alpha) / (3 + alpha) - knotxi(k)**(1 + alpha) / (1 + alpha)) ! Initial
                    end do
                end do

                do k = 1, size(prod_eta, 3) - 1 ! Loop over the polynomials
                    do l = 1, size(prod_eta, 4) ! loop over order
                        beta = size(prod_eta, 4) - l
                        eta_1(i, j) = eta_1(i, j) + prod_eta(i, j, k, l) * (knoteta(k + 1)**(1 + beta) / (1 + beta) - knoteta(k + 1)**(3 + beta) / (3 + beta)) ! Finale
                        eta_1(i, j) = eta_1(i, j) - prod_eta(i, j, k, l) * (knoteta(k)**(1 + beta) / (1 + beta) - knoteta(k)**(3 + beta) / (3 + beta)) ! Initial

                        eta_2(i, j) = eta_2(i, j) + prod_eta(i, j, k, l) * (knoteta(k + 1)**(3 + beta) / (3 + beta) - knoteta(k + 1)**(5 + beta) / (5 + beta)) ! Finale
                        eta_2(i, j) = eta_2(i, j) - prod_eta(i, j, k, l) * (knoteta(k)**(3 + beta) / (3 + beta) - knoteta(k)**(5 + beta) / (5 + beta)) ! Initial
                    end do
                end do
            end do
        end do

        result = zero

        do i = 1, size(b_xi, 1)**2 ! Loop over xi basis functions
            do j = 1, size(b_eta, 1)**2 ! Loop over eta basis functions
                i2 = indexToPair(i, size(b_xi, 1))
                j2 = indexToPair(j, size(b_eta, 1))

                result(i, j) = 2 * mppi() * (R**3) * (xi_1(i2(1), j2(1)) * eta_1(i2(2), j2(2)) - xi_2(i2(1), j2(1)) * eta_2(i2(2), j2(2)))
            end do
        end do

    end subroutine int_S22one

    subroutine init_h2plus_sep(d, n, n_remove, Z1, Z2, m, C, R, ximax, ximin, jz2, epsilon, eta_slp, save_step)
        !> @brief This subroutine initializes the B-spline coefficients for the H2+ molecule.
        !> @param d : integer : the degree of the B-spline
        !> @param n : integer : the number of usable B-splines
        !> @param n_remove : integer : the number of knots to remove from each end
        !> @param Z1 : real : the number of protons for the first atom
        !> @param Z2 : real : the number of protons for the second atom
        !> @param m : real : the mass of the electron
        !> @param C : real : the speed of light
        !> @param R : real : the distance between the two nuclei
        !> @param ximax : real : the maximum position of the B-spline on xi-axis
        !> @param ximin : real : the minimum position of the B-spline on xi-axis
        !> @param jz2 : real : the quantum number (2*jz)
        !> @param epsilon : real : the machine epsilon
        !> @param eta_slp : real : the parameter for the generation of the knot vector on eta
        !> @param save_step : boolean : whether to save matrices to a file
        type(mp_real), intent(in) :: Z1, Z2, m, C, R, ximax, ximin, epsilon, eta_slp
        integer, intent(in) :: d, n, n_remove, jz2
        logical, intent(in) :: save_step

        double precision :: tm0, tm1

        type(mp_real), dimension(:), allocatable :: knotxi, knoteta, knotxi_eps, knoteta_eps, knotxi_tmp
        type(mp_real), dimension(:, :), allocatable :: S11one, S22one, C11one, C11two, C22one, C22two, C11three, C22three, C12three, C21three, C_mat, S_mat, vect
        type(mp_real), dimension(:), allocatable :: w, fv1, fv2
        type(mp_real), dimension(:, :, :), allocatable :: bspline_xi, bspline_xi_tmp, bspline_eta
        logical :: debug_bool = .false.

        integer :: ntot, i, j, ierr
        integer :: i2(2)

        ! ----------------------
        ! Development : Partial diag
        logical :: tot_diag = .false.
        integer :: maxit = 10
        type(mp_real) :: eig

        eig = '1.87777625d4'
        ! ----------------------

        ntot = n + d + 2 * n_remove
        zero = '0.d0'
        one = '1.d0'

        write (6, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d, " and Number of BSplines to remove: ", n_remove

        print *, "Generating B-spline knots..."
        tm0 = second()
        ! Generate the knot vectors for xi and eta
        allocate (knotxi_tmp(ntot + 1), knotxi(ntot), knoteta(ntot))

        knotxi_tmp = knot_xi(d, n + 1, n_remove, ximin, ximax)
        knotxi = knotxi_tmp(1:ntot) ! Remove the last knot to avoid singularities
        knoteta = knot_eta(d, n, n_remove, eta_slp)

        tm1 = second()
        print *, "Time taken to generate knots: ", tm1 - tm0, " seconds"

        print *, "Generating B-spline coefficients..."
        tm0 = second()
        ! Generate the B-spline coefficients for xi and eta
        allocate (bspline_xi(n, ntot, d), bspline_xi_tmp(n, ntot + 1, d), bspline_eta(n, ntot, d))

        do i = 1 + n_remove, n + n_remove ! Loop over the number of B-splines
            if (debug_bool) then
                print *, "Generating B-spline coefficients for B-spline ", i, "on xi-axis..."
            end if
            call init_bspine(d, i, knotxi_tmp, ntot + 1, bspline_xi_tmp(i - n_remove, :, :), debug_bool)
            if (debug_bool) then
                print *, "Generating B-spline coefficients for B-spline ", i, "on eta-axis..."
            end if
            call init_bspine(d, i, knoteta, ntot, bspline_eta(i - n_remove, :, :), debug_bool)
        end do

        bspline_xi(:, :, :) = bspline_xi_tmp(:, 1:ntot, :)

        ! Modify the knots to avoid singularities.
        allocate (knotxi_eps(ntot), knoteta_eps(ntot))
        knotxi_eps = knotxi
        knoteta_eps = knoteta
        do i = 1, d
            knotxi_eps(i) = knotxi_eps(i) + epsilon
            knoteta_eps(i) = knoteta_eps(i) + epsilon
            knoteta_eps(ntot - i + 1) = knoteta_eps(ntot - i + 1) - epsilon
        end do

        if (debug_bool) then
            print *, "Adjusted B-spline knots for xi-axis: "
            call write_lists(knotxi_eps, 6, 25, 10)
            print *, "Adjusted B-spline knots for eta-axis: "
            call write_lists(knoteta_eps, 6, 25, 10)
        end if

        tm1 = second()
        print *, "Time taken to generate B-spline coefficients: ", tm1 - tm0, " seconds"

        print *, "Calculating the C11one integral..."
        tm0 = second()
        ! Calculate the C11one integral
        allocate (C11one(n**2, n**2))

        C11one = zero

        call int_C11one(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C11one)

        tm1 = second()
        print *, "Time taken to calculate C11one integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "C11one_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(C11one(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the C11two integral..."
        tm0 = second()
        ! Calculate the C11two integral
        allocate (C11two(n**2, n**2))

        call int_C11two(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C11two)
        tm1 = second()
        print *, "Time taken to calculate C11two integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "C11two_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(C11two(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the C22one integral..."
        tm0 = second()
        ! Calculate the C22one integral
        allocate (C22one(n**2, n**2))

        call int_C22one(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C22one)
        tm1 = second()
        print *, "Time taken to calculate C22one integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "C22one_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(C22one(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the C22two integral..."
        tm0 = second()
        ! Calculate the C22two integral
        allocate (C22two(n**2, n**2))

        call int_C22two(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C22two)
        tm1 = second()
        print *, "Time taken to calculate C22two integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "C22two_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(C22two(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the C11three integral..."
        tm0 = second()
        ! Calculate the C11three integral
        allocate (C11three(n**2, n**2))

        call int_C11three(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C11three)
        tm1 = second()
        print *, "Time taken to calculate C11three integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "C11three_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(C11three(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the C22three integral..."
        tm0 = second()
        ! Calculate the C22three integral
        allocate (C22three(n**2, n**2))

        call int_C22three(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C22three)
        tm1 = second()
        print *, "Time taken to calculate C22three integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "C22three_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(C22three(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the C12three integral..."
        tm0 = second()
        ! Calculate the C12three integral
        allocate (C12three(n**2, n**2))

        call int_C12three(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C12three)
        tm1 = second()
        print *, "Time taken to calculate C12three integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "C12three_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(C12three(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the C21three integral..."
        tm0 = second()
        ! Calculate the C21three integral
        allocate (C21three(n**2, n**2))

        call int_C21three(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, C21three)
        tm1 = second()
        print *, "Time taken to calculate C21three integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "C21three_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(C21three(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the S11one integral..."
        tm0 = second()
        ! Calculate the S11one integral
        allocate (S11one(n**2, n**2))

        call int_S11one(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, S11one)
        tm1 = second()
        print *, "Time taken to calculate S11one integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "S11one_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(S11one(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Calculating the S22one integral..."
        tm0 = second()
        ! Calculate the S22one integral
        allocate (S22one(n**2, n**2))

        call int_S22one(bspline_xi, bspline_eta, knotxi_eps, knoteta_eps, Z1, Z2, m, C, R, S22one)
        tm1 = second()
        print *, "Time taken to calculate S22one integral: ", tm1 - tm0

        if (save_step) then
            open (10, file = "S22one_sep.csv", status = 'replace')
            do i = 1, n**2
                call write_csv(S22one(i, :), 10, 30, 10)
            end do
            close (10)
        end if

        print *, "Generating the C Matrix..."
        tm0 = second()
        ! Generate the C matrix
        allocate (C_mat(4 * n**2, 4 * n**2))

        C_mat = zero
        do i = 1, n**2
            do j = 1, n**2
                ! Diagonal blocks
                C_mat(i, j) = C11one(i, j)
                C_mat(i + n**2, j + n**2) = C22one(i, j)
                C_mat(i + 2 * n**2, j + 2 * n**2) = C11two(i, j)
                C_mat(i + 3 * n**2, j + 3 * n**2) = C22two(i, j)
                ! Off-diagonal blocks
                C_mat(i, j + 2 * n**2) = C11three(i, j)
                C_mat(i, j + 3 * n**2) = C12three(i, j)
                C_mat(i + n**2, j + 2 * n**2) = C21three(i, j)
                C_mat(i + n**2, j + 3 * n**2) = C22three(i, j)
                C_mat(i + 2 * n**2, j) = C11three(j, i)
                C_mat(i + 2 * n**2, j + n**2) = C21three(j, i)
                C_mat(i + 3 * n**2, j) = C12three(j, i)
                C_mat(i + 3 * n**2, j + n**2) = C22three(j, i)
            end do
        end do

        if (save_step) then
            print *, "Saving C_mat matrix to file..."
            open (unit = 13, file = 'C_mat.csv', status = 'replace')
            do i = 1, 4 * n**2
                call write_csv(C_mat(i, :), 13, 50, 30)
            end do
            close (13)
        end if

        tm1 = second()
        print *, "Time taken to generate C matrix: ", tm1 - tm0, " seconds"

        print *, "Generating the S Matrix..."
        tm0 = second()
        ! Generate the S matrix
        allocate (S_mat(4 * n**2, 4 * n**2))

        S_mat = zero
        do i = 1, n**2
            do j = 1, n**2
                ! Diagonal blocks
                S_mat(i, j) = S11one(i, j)
                S_mat(i + n**2, j + n**2) = S22one(i, j)
                S_mat(i + 2 * n**2, j + 2 * n**2) = S11one(i, j)
                S_mat(i + 3 * n**2, j + 3 * n**2) = S22one(i, j)
            end do
        end do

        if (save_step) then
            print *, "Saving S_mat matrix to file..."
            open (unit = 14, file = 'S_mat.csv', status = 'replace')
            do i = 1, 4 * n**2
                call write_csv(S_mat(i, :), 14, 50, 30)
            end do
            close (14)
        end if

        tm1 = second()
        print *, "Time taken to generate S matrix: ", tm1 - tm0, " seconds"
        print *, "C matrix and S matrix generated successfully."

        if (debug_bool) then
            print *, "Check if C is hermitian..."
            do i = 1, 4 * n**2
                do j = 1, 4 * n**2
                    if (C_mat(i, j) - C_mat(j, i) > epsilonn(one)) then
                        print *, "C is not hermitian at (", i, ",", j, ") :"
                        call mpwrite(6, 35, 15, C_mat(i, j) - C_mat(j, i))
                    end if
                end do
            end do

            print *, "Check if S is hermitian..."
            do i = 1, 4 * n**2
                do j = 1, 4 * n**2
                    if (S_mat(i, j) - S_mat(j, i) > epsilonn(one)) then
                        print *, "S is not hermitian at (", i, ",", j, ") :"
                        call mpwrite(6, 35, 15, S_mat(i, j) - S_mat(j, i))
                    end if
                end do
            end do
        end if

        print *, "Calculating the eigenvalues..."
        tm0 = second()
        ! Calculate the eigenvalues and eigenvectors
        if (tot_diag) then
            allocate (w(4 * n**2), fv1(4 * n**2), fv2(4 * n**2), vect(4 * n**2, 4 * n**2))

            call rsg(4 * n**2, 4 * n**2, C_mat, S_mat, w, 0, vect, fv1, fv2, ierr)

            print *, "Error code: ", ierr
        else
            call pdiag(4 * n**2, C_mat, S_mat, eig, maxit)
        end if

        tm1 = second()
        print *, "Time taken to calculate eigenvalues: ", tm1 - tm0, " seconds"

        print *, "Saving logs to file..."
        ! Save logs
        open (unit = 1, file = 'log_file', status = 'replace')
        write (1, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d, " and Number of BSplines to remove: ", n_remove
        write (1, '(a)') "Speed of light: "
        call mpwrite(1, 35, 15, C)
        write (1, '(a)') "Mass of the electron: "
        call mpwrite(1, 35, 15, m)
        write (1, '(a)') "Slope for eta: "
        call mpwrite(1, 35, 15, eta_slp)
        write (1, '(a)') "Non-adjusted Xi knot vector: "
        call write_lists(knotxi, 1, 35, 15)
        write (1, '(a)') "Non-adjusted Eta knot vector: "
        call write_lists(knoteta, 1, 35, 15)
        write (1, '(a)') "Internuclear distance: "
        call mpwrite(1, 35, 15, R)
        ! write (1, '(a)') "B-Spline basis function Order: "
        ! do i = 1, n**2
        !    i2 = indexToPair(i, n)
        !    write (1, '(i4, a, i4, a, i4)', advance='no') i, " ", i2(1), " ", i2(2)
        !    write (1, '(a)') " "
        ! end do
        write (1, '(a)') "Adjusted Xi knot vector: "
        call write_lists(knotxi_eps, 1, 35, 15)
        write (1, '(a)') "Adjusted Eta knot vector: "
        call write_lists(knoteta_eps, 1, 35, 15)
        if (debug_bool) then
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "S11 integral: "
            do i = 1, n**2
                call write_lists(S11one(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "S22 integral: "
            do i = 1, n**2
                call write_lists(S22one(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "C11one integral: "
            do i = 1, n**2
                call write_lists(C11one(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "C11two integral: "
            do i = 1, n**2
                call write_lists(C11two(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "C22one integral: "
            do i = 1, n**2
                call write_lists(C22one(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "C22two integral: "
            do i = 1, n**2
                call write_lists(C22two(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "C11three integral: "
            do i = 1, n**2
                call write_lists(C11three(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "C22three integral: "
            do i = 1, n**2
                call write_lists(C22three(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "C12three integral: "
            do i = 1, n**2
                call write_lists(C12three(i, :), 1, 35, 15)
            end do
            write (1, '(a)') "--------------------------------------------------------------"
            write (1, '(a)') "C21three integral: "
            do i = 1, n**2
                call write_lists(C21three(i, :), 1, 35, 15)
            end do
        end if
        write (1, '(a)') "--------------------------------------------------------------"
        write (1, '(a)') "Value of mc²:"
        call mpwrite(1, 35, 15, m * c * c)
        write (1, '(a)') "--------------------------------------------------------------"
        if (tot_diag) then
            write (1, '(a)') "Problem totally diagonalized using the rsg subroutine."
            write (1, '(a)') "Eigenvalues: "
            do i = 1, 4 * n**2
                write (1, '(i4, a, i4)', advance = 'no') i, " "
                call mpwrite(1, 35, 15, w(4 * n**2 - i + 1))
                write (1, '(i4, a, i4)', advance = 'no') i, " Translated by -mc^2 : "
                call mpwrite(1, 35, 15, w(4 * n**2 - i + 1) - m * c * c)
                write (1, '(a)') " "
            end do
        else
            write (1, '(a)') "Problem partially diagonalized using the invsg subroutine."
            write (1, '(a)') "Eigenvalue :"
            call mpwrite(1, 35, 15, eig)
            write (1, '(a)') "Eigenvalue translated by -mc^2 :"
            call mpwrite(1, 35, 15, eig - m * c * c)
        end if
        close (1)

        print *, "Logs saved to log_file."
        print *, "Saving eigenvalues to file..."

        open (unit = 12, file = 'eigenvalues.txt', status = 'replace')
        write (12, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, " and Order of BSplines: ", d, " and Number of BSplines to remove: ", n_remove
        write (12, '(a)') "Speed of light: "
        call mpwrite(12, 35, 15, C)
        write (12, '(a)') "Mass of the electron: "
        call mpwrite(12, 35, 15, m)
        write (12, '(a)') "Slope for eta: "
        call mpwrite(12, 35, 15, eta_slp)
        write (12, '(a)') "Internuclear distance: "
        call mpwrite(12, 35, 15, R)
        write (12, '(a)') "Non-adjusted Xi knot vector: "
        call write_lists(knotxi, 12, 35, 15)
        write (12, '(a)') "Non-adjusted Eta knot vector: "
        call write_lists(knoteta, 12, 35, 15)
        write (12, '(a)') "Epsilon: "
        call mpwrite(12, 35, 15, epsilon)
        if (tot_diag) then
            do i = 1, 4 * n**2
                if (w(4 * n**2 - i + 1) > -m * c * c .and. w(4 * n**2 - i + 1) < m * c * c) then
                    write (12, '(i4, a, i4)', advance = 'no') i, " "
                    call mpwrite(12, 35, 15, w(4 * n**2 - i + 1))
                    write (12, '(i4, a, i4)', advance = 'no') i, " Translated by -mc^2 : "
                    call mpwrite(12, 35, 15, w(4 * n**2 - i + 1) - m * c * c)
                    write (12, '(a)') " "
                end if
            end do
        else
            write (12, '(i4, a, i4)', advance = 'no') 1, " "
            call mpwrite(12, 35, 15, eig)
            write (12, '(i4, a, i4)', advance = 'no') 1, " Translated by -mc^2 : "
            call mpwrite(12, 35, 15, eig - m * c * c)
            write (12, '(a)') " "
        end if
        close (12)

        print *, "Eigenvalues saved to eigenvalues.txt."

    end subroutine init_h2plus_sep

    function epsilonn(alpha)
        !> @brief Calculate the machine epsilon
        !> @param alpha The value to calculate the machine epsilon
        USE mpmodule
        implicit type(mp_real) (a - h, o - z)

        ten = '10.d0'
        epsilonn = ten**(-mpipl)

        return
    end function epsilonn

end module h2plus_sep
