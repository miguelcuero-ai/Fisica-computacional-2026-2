program programa_energia
    use modulo_romberg
    implicit none

    real(dp) :: a, b, tol
    integer :: kmax, k_final, n_eval
    real(dp), allocatable :: S(:,:)
    real(dp) :: e_dis, e_ref, diff
    real(dp) :: gamma, omega, w0_sq, x0, v0, xT, vT, E0, ET
    integer :: k, j
    
    a = 0.0_dp
    b = 3.0_dp
    tol = 1.0e-6_dp
    kmax = 20

    gamma = 0.5_dp
    omega = 2.0_dp * 4.0_dp * atan(1.0_dp)
    w0_sq = omega**2 + gamma**2

    call romberg(a, b, tol, kmax, 1, S, k_final, n_eval)
    e_dis = S(k_final, k_final)

    x0 = exp(-gamma * a) * cos(omega * a)
    v0 = -exp(-gamma * a) * (gamma * cos(omega * a) + omega * sin(omega * a))
    E0 = 0.5_dp * v0**2 + 0.5_dp * w0_sq * x0**2

    xT = exp(-gamma * b) * cos(omega * b)
    vT = -exp(-gamma * b) * (gamma * cos(omega * b) + omega * sin(omega * b))
    ET = 0.5_dp * vT**2 + 0.5_dp * w0_sq * xT**2

    e_ref = E0 - ET
    diff = abs(e_dis - e_ref)

    open(unit=10, file='tabla_energia.dat', status='replace')
    do k = 0, k_final
        write(10, *) (S(k,j), j=0, k)
    end do
    close(10)

    open(unit=20, file='resultado_energia.dat', status='replace')
    write(20, *) k_final, n_eval, e_dis, e_ref, diff
    close(20)

end program programa_energia
