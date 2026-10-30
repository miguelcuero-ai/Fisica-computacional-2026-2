! Taller 5 -- Convergencia del estimador de pi por acierto y error.
! Salidas: convergencia.dat (N, media, desviación) y pendiente.dat (pendiente m).
program taller5_pi
  implicit none
  ! Seleccionamos la precisión de doble precisión para los números decimales.
  integer, parameter :: dp = kind(1.0d0)
  ! Número de corridas independientes por cada tamaño N para medir la variabilidad.
  integer, parameter :: M = 20
  integer :: k, j, n_sem, u_dat, n_puntos
  ! Usamos enteros de 64 bits (8 bytes) para contar hasta N = 10^7 sin desbordar memoria.
  integer(8) :: i, n_total, n_aciertos
  integer, allocatable :: semilla(:)
  real(dp) :: u, v, x, y, pi_est, suma, suma2, media, desv
  real(dp) :: X_k, Y_k, sx, sy, sxx, sxy, pendiente

  ! --- Configuración inicial del generador ---
  ! Consultamos cuántos enteros necesita Fortran para guardar su estado aleatorio.
  call random_seed(size=n_sem)
  allocate(semilla(n_sem))

  ! --- Acumuladores para el ajuste lineal de mínimos cuadrados ---
  ! Guardarán sum(X), sum(Y), sum(X^2) y sum(XY) donde X = log10(N) y Y = log10(desv).
  sx = 0.0_dp; sy = 0.0_dp; sxx = 0.0_dp; sxy = 0.0_dp; n_puntos = 0

  ! Abrimos el archivo de salida; 'replace' sobreescribe datos de ejecuciones anteriores.
  open(newunit=u_dat, file='convergencia.dat', status='replace')

  ! --- Ciclo externo: recorre los tamaños de muestra N = 10^2, 10^3, ..., 10^7 ---
  do k = 2, 7
     n_total = 10_8**k                           ! N total de puntos para este nivel

     ! Reiniciamos las sumas acumuladas entre las M corridas del mismo N.
     suma = 0.0_dp; suma2 = 0.0_dp

     ! --- Repetimos el experimento M = 20 veces para este tamaño N ---
     do j = 1, M
        ! Cambiamos la semilla en cada corrida para garantizar que sean independientes.
        semilla = 1000*k + j
        call random_seed(put=semilla)

        n_aciertos = 0_8                        ! Contador de puntos dentro del círculo

        ! --- Generación de N puntos en el cuadrado [-1,1) x [-1,1) ---
        do i = 1_8, n_total
           call random_number(u)                ! u en [0,1)
           call random_number(v)                ! v en [0,1)
           ! Cambio de intervalo: transformamos de [0,1) a [-1,1)
           x = 2.0_dp*u - 1.0_dp
           y = 2.0_dp*v - 1.0_dp

           ! Si cumple x^2 + y^2 < 1, el punto está dentro del círculo unitario.
           if (x*x + y*y < 1.0_dp) n_aciertos = n_aciertos + 1_8
        end do

        ! Estimación de pi: pi ~ 4 * (aciertos / N)
        ! real(..., dp) evita una división entera que truncaría el resultado a cero.
        pi_est = 4.0_dp * real(n_aciertos, dp) / real(n_total, dp)

        ! Acumulamos la estimación y su cuadrado para hallar media y desviación muestral.
        suma  = suma + pi_est
        suma2 = suma2 + pi_est*pi_est
     end do

     ! --- Cálculo estadístico de las M = 20 corridas ---
     media = suma / real(M, dp)                 ! Media muestral (pi_bar)
     ! Desviación estándar muestral (s) entre corridas con corrección M / (M - 1).
     desv  = sqrt((real(M, dp) / real(M - 1, dp)) * (suma2 / real(M, dp) - media*media))

     ! Guardamos N, media y desviación en el archivo para la gráfica.
     write(u_dat,'(i10,2es25.16e3)') n_total, media, desv
     print '(i10,2f12.6)', n_total, media, desv

     ! --- Acumulamos datos para calcular la pendiente m ---
     ! Evaluamos la escala logarítmica: log10(s) vs log10(N)
     X_k = log10(real(n_total, dp))
     Y_k = log10(desv)
     sx  = sx + X_k
     sy  = sy + Y_k
     sxx = sxx + X_k*X_k
     sxy = sxy + X_k*Y_k
     n_puntos = n_puntos + 1
  end do
  close(u_dat)

  ! --- Ajuste por mínimos cuadrados para obtener la pendiente m ---
  ! La teoría predice s ~ 1.642 / sqrt(N), por lo que la pendiente debe ser m ~ -0.5.
  pendiente = (real(n_puntos, dp)*sxy - sx*sy) / (real(n_puntos, dp)*sxx - sx*sx)

  open(newunit=u_dat, file='pendiente.dat', status='replace')
  write(u_dat,'(es25.16e3)') pendiente
  close(u_dat)
  print '(a,f9.4)', 'Pendiente de log10(s) contra log10(N): ', pendiente
end program taller5_pi
