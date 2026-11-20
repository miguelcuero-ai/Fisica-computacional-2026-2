program programa_wien
    use modulo_ceros
    implicit none

    real(dp), parameter :: h = 6.62607015e-34_dp
    real(dp), parameter :: c = 299792458.0_dp
    real(dp), parameter :: k_B = 1.380649e-23_dp
    real(dp), parameter :: tol = 1.0e-10_dp

    real(dp) :: r_b5, r_n5, r_s5, r_b3, r_n3, r_s3, r_n5_05, b_const
    integer :: it_b5, it_n5, it_s5, it_b3, it_n3, it_s3, it_n5_05

    open(unit=10, file='iteraciones_wien.dat', status='replace')
    open(unit=20, file='resultado_wien.dat', status='replace')
    open(unit=30, file='constante_wien.dat', status='replace')

    ! --- n = 5 (registra iteraciones en iteraciones_wien.dat) ---
    call biseccion(5, 1.0_dp, 10.0_dp, tol, r_b5, it_b5, 10, 1)
    call newton(5, 10.0_dp, tol, r_n5, it_n5, 10, 2)
    call secante(5, 10.0_dp, 9.0_dp, tol, r_s5, it_s5, 10, 3)

    ! --- n = 3 ---
    call biseccion(3, 1.0_dp, 10.0_dp, tol, r_b3, it_b3, 10, 0)
    call newton(3, 10.0_dp, tol, r_n3, it_n3, 10, 0)
    call secante(3, 10.0_dp, 9.0_dp, tol, r_s3, it_s3, 10, 0)

    ! --- n = 5 (Newton con x0 = 0.5 -> raíz trivial 0) ---
    call newton(5, 0.5_dp, tol, r_n5_05, it_n5_05, 10, 0)

    ! 7 filas en resultado_wien.dat: n, metodo, raiz, iteraciones
    write(20, '(I2, 1X, I2, 1X, ES24.16, 1X, I4)') 5, 1, r_b5, it_b5
    write(20, '(I2, 1X, I2, 1X, ES24.16, 1X, I4)') 5, 2, r_n5, it_n5
    write(20, '(I2, 1X, I2, 1X, ES24.16, 1X, I4)') 5, 3, r_s5, it_s5
    write(20, '(I2, 1X, I2, 1X, ES24.16, 1X, I4)') 3, 1, r_b3, it_b3
    write(20, '(I2, 1X, I2, 1X, ES24.16, 1X, I4)') 3, 2, r_n3, it_n3
    write(20, '(I2, 1X, I2, 1X, ES24.16, 1X, I4)') 3, 3, r_s3, it_s3
    write(20, '(I2, 1X, I2, 1X, ES24.16, 1X, I4)') 5, 4, r_n5_05, it_n5_05

    ! constante_wien.dat: x5, x3, b [m K]
    b_const = (h * c) / (k_B * r_n5)
    write(30, '(ES24.16, 1X, ES24.16, 1X, ES24.16)') r_n5, r_n3, b_const

    close(10)
    close(20)
    close(30)
end program programa_wien
