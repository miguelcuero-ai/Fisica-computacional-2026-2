module modulo_romberg
    implicit none
    integer, parameter :: dp = kind(1.0d0)

contains

    real(dp) function f_energia(t) result(f)
        real(dp), intent(in) :: t
        real(dp) :: gamma, omega, v
        gamma = 0.5_dp
        omega = 2.0_dp * 4.0_dp * atan(1.0_dp)
        v = -exp(-gamma * t) * (gamma * cos(omega * t) + omega * sin(omega * t))
        f = 2.0_dp * gamma * v**2
    end function f_energia

    real(dp) function f_pi(x) result(f)
        real(dp), intent(in) :: x
        f = 4.0_dp / (1.0_dp + x**2)
    end function f_pi

    real(dp) function f_integrando(x, opcion) result(f)
        real(dp), intent(in) :: x
        integer, intent(in) :: opcion
        if (opcion == 1) then
            f = f_energia(x)
        else if (opcion == 2) then
            f = f_pi(x)
        else
            f = 0.0_dp
        end if
    end function f_integrando

    subroutine romberg(a, b, tol, kmax, opcion, S, k_final, n_eval)
        real(dp), intent(in) :: a, b, tol
        integer, intent(in) :: kmax, opcion
        real(dp), allocatable, intent(out) :: S(:,:)
        integer, intent(out) :: k_final, n_eval
        
        real(dp) :: h, suma
        integer :: k, i, j, puntos_nuevos

        if (allocated(S)) deallocate(S)
        allocate(S(0:kmax, 0:kmax))
        S = 0.0_dp

        h = b - a
        S(0,0) = (h / 2.0_dp) * (f_integrando(a, opcion) + f_integrando(b, opcion))
        n_eval = 2

        do k = 1, kmax
            h = h / 2.0_dp
            puntos_nuevos = 2**(k-1)
            suma = 0.0_dp
            
            do i = 1, puntos_nuevos
                suma = suma + f_integrando(a + real(2*i - 1, dp) * h, opcion)
            end do
            
            n_eval = n_eval + puntos_nuevos
            S(k,0) = S(k-1,0) / 2.0_dp + h * suma
            
            do j = 1, k
                S(k,j) = S(k,j-1) + (S(k,j-1) - S(k-1,j-1)) / (4.0_dp**j - 1.0_dp)
            end do
            
            if (abs(S(k,k) - S(k-1,k-1)) < tol) then
                k_final = k
                return
            end if
        end do
        k_final = kmax
    end subroutine romberg

end module modulo_romberg
