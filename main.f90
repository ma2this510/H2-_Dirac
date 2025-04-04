program main
    use mpmodule
    use bspline_gen
    use h2plus
    implicit none

    integer :: d, n_remove, n, i
    type(mp_real) :: Z1, Z2, Z, Zmax, Zmin, Rmin, Rmax, c, zero, one

    zero = '0.d0'
    one = '1.d0'

    !-------------------------------------------------
    ! Define Important Constants
    d = 3 ! Order of the B-Spline (order Mathematica + 1)
    n = 20 ! Number of B-Spline knots per region
    n_remove = 2 ! Number of knots to remove from each end
    Z1 = '2.0d0' ! number of protons for the first atom
    Z2 = '3.0d0' ! number of protons for the second atom
    C = '137.0359895d0' ! check CODATA 1986
    Z = '1.0d1' ! position of both nuclei
    Zmax = '50.0d0' ! maximum position of the B-spline on z-axis
    Zmin = '1.0d-4' ! minimum position of the B-spline on z-axis
    Rmin = '1.0d-4' ! minimum position of the B-spline on r-axis
    Rmax = '1.0d1' ! maximum position of the B-spline on r-axis
    !-------------------------------------------------

    call init_h2plus(d, n, n_remove, Z1, Z2, Z, Zmax, Zmin, Rmin, Rmax)
end program main

function epsilonn(alpha)
    !> @brief Calculate the machine epsilon
    !> @param alpha The value to calculate the machine epsilon
    USE mpmodule
    implicit type(mp_real) (a - h, o - z)
 
    ten = '10.d0'
    epsilonn = ten**(-mpipl)
 
    return
 end function epsilonn
 