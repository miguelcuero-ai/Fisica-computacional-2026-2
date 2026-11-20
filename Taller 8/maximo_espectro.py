import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# --- Constantes ---
h = 6.62607015e-34
c = 299792458.0
kB = 1.380649e-23
T = 5772.0

# --- 1. Leer constante_wien.dat ---
data_wien = np.loadtxt('constante_wien.dat')
x5 = data_wien[0]
x3 = data_wien[1]
b_wien = data_wien[2]

# Valores teóricos esperados (en nm)
lam_max_teo = (b_wien / T) * 1e9
lam_max_freq = (h * c / (kB * T * x3)) * 1e9

# --- Función de cuerpo negro (Planck) ---
def planck(lam_nm):
    lam = lam_nm * 1e-9
    return (2.0 * h * c**2 / lam**5) / (np.exp(h * c / (lam * kB * T)) - 1.0)

dlams = [100.0, 50.0, 25.0, 12.5]
resultados = []
lam_plot, B_plot, spline_plot, lam_max_plot = None, None, None, None

for dlam in dlams:
    lam = np.arange(100.0, 3000.0 + dlam/2, dlam)
    B = planck(lam)

    # (a) Mayor valor del arreglo
    idx_a = np.argmax(B)
    lam_a = lam[idx_a]
    err_a = abs(lam_a - lam_max_teo)

    # (b) Cero del spline de la derivada central
    lam_int = lam[1:-1]
    dB_int = (B[2:] - B[:-2]) / (2.0 * dlam)

    spline_b = CubicSpline(lam_int, dB_int)
    d2_spline_b = spline_b.derivative()
    raices_b = spline_b.roots()
    valid_b = [r for r in raices_b if 100 <= r <= 3000 and d2_spline_b(r) < 0]
    lam_b = min(valid_b, key=lambda x: abs(x - lam_a)) if len(valid_b) > 0 else np.nan
    err_b = abs(lam_b - lam_max_teo) if not np.isnan(lam_b) else np.nan

    # (c) Cero de la derivada del spline de B_lam
    spline_c = CubicSpline(lam, B)
    dspline_c = spline_c.derivative()
    d2_spline_c = dspline_c.derivative()
    raices_c = dspline_c.roots()
    valid_c = [r for r in raices_c if 100 <= r <= 3000 and d2_spline_c(r) < 0]
    lam_c = min(valid_c, key=lambda x: abs(x - lam_a)) if len(valid_c) > 0 else np.nan
    err_c = abs(lam_c - lam_max_teo) if not np.isnan(lam_c) else np.nan

    resultados.append([dlam, err_a, err_b, err_c])

    if dlam == 25.0:
        lam_plot = lam
        B_plot = B
        spline_plot = spline_c
        lam_max_plot = lam_c

# --- Escribir maximo_espectro.dat ---
with open('maximo_espectro.dat', 'w') as f:
    for res in resultados:
        f.write(f"{res[0]:.1f} {res[1]:.16e} {res[2]:.16e} {res[3]:.16e}\n")

# --- Gráfica para malla de 25 nm ---
lam_fino = np.linspace(100, 3000, 1000)
plt.figure(figsize=(9, 5))
plt.plot(lam_plot, B_plot, 'k.', label=r'Datos (Malla $\Delta\lambda = 25$ nm)')
plt.plot(lam_fino, spline_plot(lam_fino), 'C0-', alpha=0.7, label=r'Spline $B_\lambda$')
plt.axvline(lam_max_plot, color='C3', linestyle='--', label=fr'$\lambda_{{max}}$ = {lam_max_plot:.2f} nm')
plt.axvline(lam_max_freq, color='C2', linestyle='-.', label=fr'$\lambda_{{max, frec}}$ = {lam_max_freq:.2f} nm')

plt.xlabel(r'Longitud de onda $\lambda$ (nm)')
plt.ylabel(r'Radiancia espectral $B_\lambda$')
plt.title('Espectro Solar y Máximos de Emisión (T = 5772 K)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.savefig('espectro_solar.png', dpi=300, bbox_inches='tight')
