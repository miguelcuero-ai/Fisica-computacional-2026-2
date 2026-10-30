! Taller 5 -- Integral de sin(x) en [0,pi/2) por valor medio y acierto y error.
! Salida: seno.dat (media_vm, desv_vm, media_ae, desv_ae).
program taller5_seno
  implicit none
  ! Seleccionamos doble precisión para los reales.
  integer, parameter :: dp = kind(1.0d0)
  ! N = 1000 puntos por cada simulación; R = 1000 repeticiones para evaluar la dispersión.
  integer, parameter :: N = 1000, R = 1000
  ! Límites de integración: a = 0, b = pi/2 (usando 2*atan(1) para precisión exactitud de pi/2) y H = 1 de altura.
  real(dp), parameter :: a = 0.0_dp, b = 2.0_dp*atan(1.0_dp), H = 1.0_dp
  integer :: rep, i, n_bajo, n_sem, u_dat
  integer, allocatable :: semilla(:)
  real(dp) :: u, v, x, y, suma_f, I_vm, I_ae
  real(dp) :: s_vm, s2_vm, s_ae, s2_ae, media_vm, desv_vm, media_ae, desv_ae

  ! --- Inicialización de la semilla ---
  ! Consultamos el número de enteros requeridos para el estado del generador.
  call random_seed(size=n_sem)
  allocate(semilla(n_sem))
  semilla = 314
  ! Fijamos la semilla una sola vez antes de iniciar las R corridas repetidas.
  call random_seed(put=semilla)

  ! Acumuladores de las estimaciones y sus cuadrados para ambos métodos.
  s_vm = 0.0_dp; s2_vm = 0.0_dp; s_ae = 0.0_dp; s2_ae = 0.0_dp

  ! --- Ciclo principal: realizamos R repeticiones de la simulación ---
  do rep = 1, R

     ! =================================================================
     ! MÉTODO 1: Valor medio
     ! Evaluamos la función continua sin(x) en N puntos uniformes en [a,b).
     ! =================================================================
     suma_f = 0.0_dp
     do i = 1, N
        call random_number(u)                   ! u en [0,1)
        x = a + (b - a)*u                       ! Cambio de intervalo: [0,1) -> [0, pi/2)
        suma_f = suma_f + sin(x)                ! Acumulamos f(x) = sin(x)
     end do
     ! Integral = ancho del intervalo * promedio de la función f(x)
     I_vm = (b - a) * suma_f / real(N, dp)

     ! =================================================================
     ! MÉTODO 2: Acierto y error
     ! Generamos puntos (x,y) uniformes en el rectángulo [a,b) x [0,H).
     ! =================================================================
     n_bajo = 0                                 ! Contador de aciertos (puntos bajo la curva)
     do i = 1, N
        call random_number(u)
        call random_number(v)
        x = a + (b - a)*u                       ! Coordenada x en [0, pi/2)
        y = H * v                               ! Coordenada y en [0, 1)
        ! Si y < sin(x), el punto se encuentra por debajo de la curva.
        if (y < sin(x)) n_bajo = n_bajo + 1
     end do
     ! Integral = área del rectángulo envolvente * fracción de aciertos
     I_ae = (b - a) * H * real(n_bajo, dp) / real(N, dp)

     ! --- Acumulamos resultados para calcular media y varianza muestral ---
     s_vm  = s_vm + I_vm
     s2_vm = s2_vm + I_vm*I_vm
     s_ae  = s_ae + I_ae
     s2_ae = s2_ae + I_ae*I_ae
  end do

  ! --- Estadística final de las R = 1000 repeticiones ---
  ! Calculamos la media muestral de cada estimador (ambas deben estar cerca de 1.0)
  media_vm = s_vm / real(R, dp)
  media_ae = s_ae / real(R, dp)

  ! Calculamos la desviación estándar muestral entre corridas con la corrección R / (R - 1)
  desv_vm  = sqrt((real(R, dp) / real(R - 1, dp)) * (s2_vm / real(R, dp) - media_vm*media_vm))
  desv_ae  = sqrt((real(R, dp) / real(R - 1, dp)) * (s2_ae / real(R, dp) - media_ae*media_ae))

  ! Guardamos los cuatro resultados numéricos ordenados para la lectura automática del script de verificación.
  open(newunit=u_dat, file='seno.dat', status='replace')
  write(u_dat,'(4es25.16e3)') media_vm, desv_vm, media_ae, desv_ae
  close(u_dat)

  ! Mostramos las estadísticas en consola
  print '(a,2f10.5)', 'Valor medio (media, desviacion): ', media_vm, desv_vm
  print '(a,2f10.5)', 'Acierto y error (media, desviacion): ', media_ae, desv_ae
end program taller5_seno
