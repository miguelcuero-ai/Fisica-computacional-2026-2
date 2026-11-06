! Taller 6 -- Programa principal
program programa_nucleo
  use modulo_nucleo
  implicit none

  integer :: u_nuc, u_val, u_max
  integer :: i, A, Z
  real(dp) :: E_B_valle, E_B_global_min, E_B_nuc, M_nuc
  integer :: A_mas_ligado, Z_mas_ligado

  ! Arreglos para los 8 núcleos del enunciado en el orden exacto
  integer, dimension(8), parameter :: Z_8 = (/ 94, 32, 72, 18, 26, 20, 92, 92 /)
  integer, dimension(8), parameter :: A_8 = (/ 244, 70, 179, 38, 56, 40, 235, 238 /)

  ! 1. Escribir nucleos.dat
  open(newunit=u_nuc, file='nucleos.dat', status='replace')
  do i = 1, 8
     Z = Z_8(i)
     A = A_8(i)
     E_B_nuc = E_B(Z, A)
     M_nuc = M(Z, A)
     write(u_nuc, '(2i6, 5f16.8, f18.4)') Z, A, &
          termino_volumen(A), &
          termino_superficie(A), &
          termino_coulomb(Z, A), &
          termino_asimetria(Z, A), &
          E_B_nuc, M_nuc
  end do
  close(u_nuc)

  ! 2. Escribir valle.dat y encontrar el máximo global
  open(newunit=u_val, file='valle.dat', status='replace')
  
  ! Inicialización para buscar el mínimo global de E_B
  A_mas_ligado = 4
  call Z_estable(4, Z_mas_ligado, E_B_global_min)
  
  ! Escribimos para A = 4
  write(u_val, '(2i6, f16.8)') 4, Z_mas_ligado, E_B_global_min

  do A = 5, 260
     call Z_estable(A, Z, E_B_valle)
     write(u_val, '(2i6, f16.8)') A, Z, E_B_valle
     
     if (E_B_valle < E_B_global_min) then
        E_B_global_min = E_B_valle
        A_mas_ligado = A
        Z_mas_ligado = Z
     end if
  end do
  close(u_val)

  ! 3. Escribir maximo.dat
  open(newunit=u_max, file='maximo.dat', status='replace')
  write(u_max, '(2i6, f16.8)') A_mas_ligado, Z_mas_ligado, E_B_global_min
  close(u_max)

end program programa_nucleo
