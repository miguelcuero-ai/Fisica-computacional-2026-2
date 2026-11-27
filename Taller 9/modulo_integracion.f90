! ============================================================
! Módulo:    integracion
! Propósito: reglas compuestas del trapecio y de Simpson para
!            una función muestreada en n puntos equiespaciados.
!            Se compila como librería de Python con F2PY:
!            python -m numpy.f2py -c modulo_integracion.f90 -m integ
! ============================================================
module integracion
  implicit none
contains

  ! trapecio: regla compuesta del trapecio
  ! Entradas: n, número de puntos; y(n), valores de la función; h, paso
  ! Salida:   s, integral aproximada
  ! Caso inválido (n < 2): s = 0 (sin stop, para no tumbar el intérprete)
  subroutine trapecio(n, y, h, s)
    integer, intent(in)  :: n
    real(8), intent(in)  :: y(n), h
    real(8), intent(out) :: s

    if (n < 2) then
      s = 0.0d0
      return
    end if
    ! T = h*[(y1 + yn)/2 + suma de los puntos interiores]
    s = h*(0.5d0*(y(1) + y(n)) + sum(y(2:n-1)))
  end subroutine trapecio

  ! simpson: regla compuesta de Simpson (n impar, es decir, intervalos pares)
  ! Entradas: n, número de puntos (impar); y(n), valores; h, paso
  ! Salida:   s, integral aproximada
  ! Caso inválido (n < 3 o n par): s = 0 (sin stop)
  subroutine simpson(n, y, h, s)
    integer, intent(in)  :: n
    real(8), intent(in)  :: y(n), h
    real(8), intent(out) :: s

    if (n < 3 .or. mod(n, 2) == 0) then
      s = 0.0d0
      return
    end if
    ! S = h/3*[y1 + yn + 4*(i par: 2,4,...,n-1) + 2*(i impar: 3,5,...,n-2)]
    s = h/3.0d0*(y(1) + y(n) + 4.0d0*sum(y(2:n-1:2)) + 2.0d0*sum(y(3:n-2:2)))
  end subroutine simpson

end module integracion
