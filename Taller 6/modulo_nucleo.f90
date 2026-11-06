! Taller 6 -- Módulo de la fórmula semiempírica (escriba aquí su módulo)
module modulo_nucleo
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! Constantes del modelo (Zubairi & Weber)
  real(dp), parameter :: a_V = 15.8_dp
  real(dp), parameter :: a_S = 18.3_dp
  real(dp), parameter :: a_C = 0.71_dp
  real(dp), parameter :: a_A = 92.7_dp
  real(dp), parameter :: M_n = 939.57_dp
  real(dp), parameter :: M_H = 938.78_dp

contains

  ! Término de Volumen: -a_V
  real(dp) function termino_volumen(A)
    integer, intent(in) :: A
    termino_volumen = -a_V
  end function termino_volumen

  ! Término de Superficie: a_S * A^(-1/3)
  real(dp) function termino_superficie(A)
    integer, intent(in) :: A
    termino_superficie = a_S * (real(A, dp)**(-1.0_dp / 3.0_dp))
  end function termino_superficie

  ! Término de Coulomb: a_C * Z^2 * A^(-4/3)
  real(dp) function termino_coulomb(Z, A)
    integer, intent(in) :: Z, A
    termino_coulomb = a_C * (real(Z, dp)**2) * (real(A, dp)**(-4.0_dp / 3.0_dp))
  end function termino_coulomb

  ! Término de Asimetría: a_A * (Z - A/2)^2 / A^2
  real(dp) function termino_asimetria(Z, A)
    integer, intent(in) :: Z, A
    termino_asimetria = a_A * ((real(Z, dp) - real(A, dp)/2.0_dp)**2) / (real(A, dp)**2)
  end function termino_asimetria

  ! Energía de enlace por nucleón E_B(Z, A)
  real(dp) function E_B(Z, A)
    integer, intent(in) :: Z, A
    E_B = termino_volumen(A) + termino_superficie(A) + &
          termino_coulomb(Z, A) + termino_asimetria(Z, A)
  end function E_B

  ! Masa atómica aproximada M(Z, A)
  real(dp) function M(Z, A)
    integer, intent(in) :: Z, A
    M = real(Z, dp)*M_H + real(A - Z, dp)*M_n + real(A, dp)*E_B(Z, A)
  end function M

  ! Subrutina para buscar el Z de máximo enlace y su energía de enlace E_B para un A dado
  subroutine Z_estable(A, Z_estable_valle, E_B_minimo)
    integer, intent(in) :: A
    integer, intent(out) :: Z_estable_valle
    real(dp), intent(out) :: E_B_minimo
    integer :: Z_prueba
    real(dp) :: E_B_actual

    ! Inicializamos con el primer protón posible (Z=1)
    Z_estable_valle = 1
    E_B_minimo = E_B(1, A)

    do Z_prueba = 2, A
       E_B_actual = E_B(Z_prueba, A)
       ! Buscamos el núcleo más ligado, que corresponde al valor de E_B más negativo (mínimo algebraico)
       if (E_B_actual < E_B_minimo) then
          E_B_minimo = E_B_actual
          Z_estable_valle = Z_prueba
       end if
    end do
  end subroutine Z_estable

end module modulo_nucleo
