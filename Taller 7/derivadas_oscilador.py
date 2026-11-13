import numpy as np
import matplotlib.pyplot as plt

gamma = np.float64(0.5)
omega = np.float64(2.0 * np.pi)
t_fin = np.float64(3.0)

def exact_v(t):
    return -np.exp(-gamma * t) * (gamma * np.cos(omega * t) + omega * np.sin(omega * t))

def exact_a(t):
    return np.exp(-gamma * t) * ((gamma**2 - omega**2) * np.cos(omega * t) + 2.0 * gamma * omega * np.sin(omega * t))

def calcular_derivadas(dt):
    dt = np.float64(dt)
    N = int(np.round(t_fin / dt)) + 1
    t = np.linspace(0.0, t_fin, N, dtype=np.float64)
    x = np.exp(-gamma * t) * np.cos(omega * t)
    
    v_adelantada = (x[1:] - x[:-1]) / dt
    err_fw = np.max(np.abs(v_adelantada - exact_v(t[:-1])))
    
    v_atrasada = (x[1:] - x[:-1]) / dt
    err_bw = np.max(np.abs(v_atrasada - exact_v(t[1:])))
    
    v_central = (x[2:] - x[:-2]) / (2.0 * dt)
    err_cen_v = np.max(np.abs(v_central - exact_v(t[1:-1])))
    
    v_5puntos = (-x[4:] + 8.0*x[3:-1] - 8.0*x[1:-3] + x[:-4]) / (12.0 * dt)
    err_5p = np.max(np.abs(v_5puntos - exact_v(t[2:-2])))
    
    a_central = (x[2:] - 2.0*x[1:-1] + x[:-2]) / (dt**2)
    err_cen_a = np.max(np.abs(a_central - exact_a(t[1:-1])))
    
    return t, x, v_central, a_central, err_fw, err_bw, err_cen_v, err_5p, err_cen_a

t1, x1, v_cen1, a_cen1, e_fw1, e_bw1, e_cv1, e_5p1, e_ca1 = calcular_derivadas(0.01)
_, _, _, _, e_fw2, e_bw2, e_cv2, e_5p2, e_ca2 = calcular_derivadas(0.005)

datos_derivadas = np.column_stack((t1[1:-1], x1[1:-1], v_cen1, a_cen1))
np.savetxt('derivadas_oscilador.dat', datos_derivadas, fmt='%.15e')

datos_errores = np.array([
    [0.01, e_fw1, e_bw1, e_cv1, e_5p1, e_ca1],
    [0.005, e_fw2, e_bw2, e_cv2, e_5p2, e_ca2]
])
np.savetxt('errores_derivadas.dat', datos_errores, fmt='%.15e')

plt.figure(figsize=(10, 6))
t_exact = np.linspace(0, t_fin, 1000)
plt.plot(t_exact, exact_v(t_exact), 'b-', label='$v$ exacta')
plt.plot(t_exact, exact_a(t_exact), 'r-', label='$a$ exacta')
plt.plot(t1[1:-1], v_cen1, 'b--', label='$v$ central (numérica)')
plt.plot(t1[1:-1], a_cen1, 'r--', label='$a$ central (numérica)')
plt.title('Velocidad y Aceleración del Oscilador Amortiguado')
plt.xlabel('Tiempo (s)')
plt.ylabel('Amplitud')
plt.legend()
plt.grid(True)
plt.savefig('derivadas_oscilador.png', dpi=300)
