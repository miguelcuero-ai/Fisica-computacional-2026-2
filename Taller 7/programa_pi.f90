program programa_pi
    use modulo_romberg
    implicit none

    real(dp) :: a, b, tol
    integer :: kmax, k_final, n_eval
    real(dp), allocatable :: S(:,:)
    real(dp) :: pi_calc, pi_exact, error_pi

    a = 0.0_dp
    b = 1.0_dp
    tol = 1.0e-5_dp
    kmax = 20

    call romberg(a, b, tol, kmax, 2, S, k_final, n_eval)

    pi_calc = S(k_final, k_final)
    pi_exact = 4.0_dp * atan(1.0_dp)
    error_pi = abs(pi_calc - pi_exact)

    open(unit=30, file='resultado_pi.dat', status='replace')
    write(30, *) k_final, n_eval, pi_calc, error_pi
    close(30)

end program programa_pi
