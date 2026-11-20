module modulo_romberg
    implicit none
    integer, parameter :: dp = kind(0.d0)

contains

    real(dp) function int_cv(t)
        real(dp), intent(in) :: t
        real(dp) :: x
        if (t <= 0.0_dp .or. t >= 1.0_dp) then
            int_cv = 0.0_dp
        else
            x = t / (1.0_dp - t)
            if (x > 500.0_dp) then
                int_cv = 0.0_dp
            else
                int_cv = (x**3 / (exp(x) - 1.0_dp)) / ((1.0_dp - t)**2)
            end if
        end if
    end function int_cv

    real(dp) function int_scv(x)
        real(dp), intent(in) :: x
        if (x <= 0.0_dp .or. x > 500.0_dp) then
            int_scv = 0.0_dp
        else
            int_scv = x**3 / (exp(x) - 1.0_dp)
        end if
    end function int_scv

    subroutine romberg(f, a, b, tol, integral, k_final, evals)
        interface
            real(dp) function f(x)
                import dp
                real(dp), intent(in) :: x
            end function f
        end interface
        real(dp), intent(in) :: a, b, tol
        real(dp), intent(out) :: integral
        integer, intent(out) :: k_final, evals

        integer, parameter :: max_k = 20
        real(dp) :: R(0:max_k, 0:max_k)
        real(dp) :: h, suma
        integer :: k, j, i, n_puntos

        h = b - a
        R(0,0) = 0.5_dp * h * (f(a) + f(b))
        evals = 2

        do k = 1, max_k
            n_puntos = 2**(k-1)
            h = (b - a) / real(2**k, dp)
            suma = 0.0_dp

            do i = 1, n_puntos
                suma = suma + f(a + real(2*i - 1, dp) * h)
            end do

            R(k,0) = 0.5_dp * R(k-1,0) + h * suma
            evals = evals + n_puntos

            do j = 1, k
                R(k,j) = R(k,j-1) + (R(k,j-1) - R(k-1,j-1)) / (4.0_dp**j - 1.0_dp)
            end do

            if (abs(R(k,k) - R(k-1,k-1)) < tol) then
                integral = R(k,k)
                k_final = k
                return
            end if
        end do

        integral = R(max_k, max_k)
        k_final = max_k
    end subroutine romberg

end module modulo_romberg
