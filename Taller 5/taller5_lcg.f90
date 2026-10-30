! Taller 5 -- pi por acierto y error con el LCG (57, 1, 256) en lugar de random_number.
! Salida: lcg.dat (N, pi_est).
program taller5_lcg
  implicit none
  ! Definimos precisión doble para realizar los cálculos decimales.
  integer, parameter :: dp = kind(1.0d0)
  ! Parámetros del Generador Lineal Congruente (LCG): a = 57, c = 1, M = 256.
  ! M = 256 significa que el estado del generador solo toma 256 valores enteros distintos.
  integer, parameter :: a_lcg = 57, c_lcg = 1, M_lcg = 256
  integer :: k, r, u_dat
  ! Usamos enteros de 64 bits (8 bytes) para permitir conteos hasta N = 10^7.
  integer(8) :: i, n_total, n_aciertos
  real(dp) :: x, y, pi_est

  ! Abrimos el archivo lcg.dat para registrar los resultados de la simulación.
  open(newunit=u_dat, file='lcg.dat', status='replace')

  ! --- Ciclo principal: evaluamos la estimación para N = 10^2, 10^3, ..., 10^7 ---
  do k = 2, 7
     n_total = 10_8**k                           ! N total de puntos a generar
     r = 10                                     ! Semilla inicial del LCG, la misma para cada N
     n_aciertos = 0_8                        ! Reiniciamos el contador de puntos dentro del círculo

     ! --- Generación de N puntos (x, y) usando parejas consecutivas del LCG ---
     do i = 1_8, n_total
        ! Primer avance del estado del LCG para obtener la coordenada x
        r = mod(a_lcg*r + c_lcg, M_lcg)
        ! Normalizamos r/M a [0,1) y luego aplicamos el cambio de intervalo a [-1,1)
        x = 2.0_dp * real(r, dp) / real(M_lcg, dp) - 1.0_dp

        ! Segundo avance del estado del LCG para obtener la coordenada y
        r = mod(a_lcg*r + c_lcg, M_lcg)
        y = 2.0_dp * real(r, dp) / real(M_lcg, dp) - 1.0_dp

        ! Verificamos si el punto cae dentro del círculo de radio 1 (x^2 + y^2 < 1)
        if (x*x + y*y < 1.0_dp) n_aciertos = n_aciertos + 1_8
     end do

     ! Estimación de pi con los puntos generados por el LCG
     pi_est = 4.0_dp * real(n_aciertos, dp) / real(n_total, dp)

     ! Guardamos N y pi_est en el archivo e imprimimos por pantalla
     write(u_dat,'(i10,es25.16e3)') n_total, pi_est
     print '(i10,f12.6)', n_total, pi_est
  end do
  close(u_dat)                                  ! Cerramos el archivo lcg.dat
end program taller5_lcg
