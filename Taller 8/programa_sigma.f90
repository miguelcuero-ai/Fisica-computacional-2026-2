program programa_sigma
    use modulo_romberg
    implicit none

    real(dp), parameter :: h_pl = 6.62607015e-34_dp
    real(dp), parameter :: c = 299792458.0_dp
    real(dp), parameter :: k_B = 1.380649e-23_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tol = 1.0e-8_dp

    real(dp) :: I_cv, I_scv, sigma_cv, sigma_scv, prefactor
    integer :: k_cv, k_scv, ev_cv, ev_scv

    open(unit=40, file='resultado_sigma.dat', status='replace')

    prefactor = (2.0_dp * pi * k_B**4) / (h_pl**3 * c**2)

    ! Fila 1: Con cambio de variable t in [0, 1]
    call romberg(int_cv, 0.0_dp, 1.0_dp, tol, I_cv, k_cv, ev_cv)
    sigma_cv = prefactor * I_cv

    ! Fila 2: Sin cambio de variable x in [0, 100] -> parada prematura en k = 1
    call romberg(int_scv, 0.0_dp, 100.0_dp, tol, I_scv, k_scv, ev_scv)
    sigma_scv = prefactor * I_scv

    write(40, '(I4, 1X, I8, 1X, ES24.16, 1X, ES24.16)') k_cv, ev_cv, I_cv, sigma_cv
    write(40, '(I4, 1X, I8, 1X, ES24.16, 1X, ES24.16)') k_scv, ev_scv, I_scv, sigma_scv

    close(40)
end program programa_sigma
