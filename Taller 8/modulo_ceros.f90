module modulo_ceros
    implicit none
    integer, parameter :: dp = kind(0.d0)

contains

    real(dp) function f_wien(x, n)
        real(dp), intent(in) :: x
        integer, intent(in) :: n
        f_wien = x - real(n, dp) * (1.0_dp - exp(-x))
    end function f_wien

    real(dp) function df_wien(x, n)
        real(dp), intent(in) :: x
        integer, intent(in) :: n
        df_wien = 1.0_dp - real(n, dp) * exp(-x)
    end function df_wien

    subroutine biseccion(n_val, a, b, tol, raiz, iter, iunit, id_metodo)
        integer, intent(in) :: n_val, iunit, id_metodo
        real(dp), intent(in) :: a, b, tol
        real(dp), intent(out) :: raiz
        integer, intent(out) :: iter

        real(dp) :: fa, fb, fc, a_local, b_local, c
        integer :: k

        a_local = a
        b_local = b
        fa = f_wien(a_local, n_val)
        fb = f_wien(b_local, n_val)
        iter = 0

        if (fa * fb > 0.0_dp) then
            print *, "Error en Biseccion: No hay cambio de signo."
            return
        end if

        do k = 1, 100
            c = 0.5_dp * (a_local + b_local)
            fc = f_wien(c, n_val)
            iter = k

            if (id_metodo > 0) write(iunit, '(I2, 1X, I4, 1X, ES24.16)') id_metodo, iter, c

            ! Criterio de parada exacto: b - c <= tol
            if ((b_local - c) <= tol) then
                raiz = c
                return
            end if

            if (fa * fc < 0.0_dp) then
                b_local = c
                fb = fc
            else
                a_local = c
                fa = fc
            end if
        end do
        print *, "Error en Biseccion: Limite de 100 iteraciones."
    end subroutine biseccion

    subroutine newton(n_val, x0, tol, raiz, iter, iunit, id_metodo)
        integer, intent(in) :: n_val, iunit, id_metodo
        real(dp), intent(in) :: x0, tol
        real(dp), intent(out) :: raiz
        integer, intent(out) :: iter

        real(dp) :: x_act, x_new, deriv
        integer :: k

        x_act = x0
        iter = 0

        do k = 1, 100
            deriv = df_wien(x_act, n_val)
            if (abs(deriv) < epsilon(1.0_dp)) then
                print *, "Error en Newton: Denominador nulo."
                return
            end if

            x_new = x_act - f_wien(x_act, n_val) / deriv
            iter = k

            if (id_metodo > 0) write(iunit, '(I2, 1X, I4, 1X, ES24.16)') id_metodo, iter, x_new

            ! Criterio de parada exacto: |x_{k+1} - x_k| < tol
            if (abs(x_new - x_act) < tol) then
                raiz = x_new
                return
            end if
            x_act = x_new
        end do
        print *, "Error en Newton: Limite de 100 iteraciones."
    end subroutine newton

    subroutine secante(n_val, x0, x1, tol, raiz, iter, iunit, id_metodo)
        integer, intent(in) :: n_val, iunit, id_metodo
        real(dp), intent(in) :: x0, x1, tol
        real(dp), intent(out) :: raiz
        integer, intent(out) :: iter

        real(dp) :: xa, xb, x_new, fa, fb, den
        integer :: k

        xa = x0
        xb = x1
        fa = f_wien(xa, n_val)
        fb = f_wien(xb, n_val)
        iter = 0

        do k = 1, 100
            den = fb - fa
            if (abs(den) < epsilon(1.0_dp)) then
                print *, "Error en Secante: Denominador nulo."
                return
            end if

            x_new = xb - fb * (xb - xa) / den
            iter = k

            if (id_metodo > 0) write(iunit, '(I2, 1X, I4, 1X, ES24.16)') id_metodo, iter, x_new

            ! Criterio de parada exacto: |x_{k+1} - x_k| < tol
            if (abs(x_new - xb) < tol) then
                raiz = x_new
                return
            end if

            xa = xb
            fa = fb
            xb = x_new
            fb = f_wien(x_new, n_val)
        end do
        print *, "Error en Secante: Limite de 100 iteraciones."
    end subroutine secante

end module modulo_ceros
