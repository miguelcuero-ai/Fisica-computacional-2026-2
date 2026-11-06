"""Taller 6 -- Módulo de análisis (escriba aquí sus cuatro funciones)."""
import numpy as np
import matplotlib.pyplot as plt

def leer_valle(archivo):
    """Lee el archivo valle.dat escrito por el programa de Fortran.

    Entrada:
        archivo: ruta o nombre del archivo de entrada
    Salidas:
        A: arreglo numpy con los valores de A (número de nucleones)
        Z: arreglo numpy con los valores de Z (número de protones)
        E_B: arreglo numpy con las energías de enlace por nucleón (MeV)
    """
    datos = np.loadtxt(archivo)
    A = datos[:, 0]
    Z = datos[:, 1]
    E_B = datos[:, 2]
    return A, Z, E_B

def encontrar_maximo_ligado(A, Z, E_B):
    """Encuentra el núcleo con el mayor enlace por nucleón (E_B más negativo)
    recorriendo los datos con un ciclo y un condicional (sin min, np.min ni np.argmin).

    Entradas:
        A: arreglo numpy con el número de nucleones
        Z: arreglo numpy con el número de protones
        E_B: arreglo numpy con las energías de enlace por nucleón
    Salidas:
        A_max: número de nucleones del núcleo más ligado (float)
        Z_max: número de protones del núcleo más ligado (float)
        E_B_max: energía de enlace por nucleón de dicho núcleo (float)
    """
    A_max = A[0]
    Z_max = Z[0]
    E_B_max = E_B[0]
    
    for i in range(1, len(E_B)):
        # Buscamos el E_B más negativo (el mínimo algebraico)
        if E_B[i] < E_B_max:
            E_B_max = E_B[i]
            A_max = A[i]
            Z_max = Z[i]
            
    return A_max, Z_max, E_B_max

def escribir_maximo(archivo, A, Z, E_B):
    """Escribe el núcleo más ligado en el archivo de salida con el mismo formato de maximo.dat.

    Entradas:
        archivo: ruta del archivo de salida (maximo_python.dat)
        A: número de nucleones (float o int)
        Z: número de protones (float o int)
        E_B: energía de enlace por nucleón
    """
    with open(archivo, 'w') as f:
        f.write(f"{int(A):6d}{int(Z):6d}{E_B:16.8f}\n")

def graficar_curva(A, Z, E_B, experimentales, archivo_salida):
    """Grafica B/A = -E_B contra A a lo largo de la curva de máximo enlace junto con
    los puntos experimentales del AME2020, guardando el resultado como imagen PNG.

    Entradas:
        A: arreglo numpy de nucleones
        Z: arreglo numpy de protones
        E_B: arreglo numpy de energías de enlace por nucleón (negativas)
        experimentales: lista de tuplas (símbolo, Z, A, B/A [MeV])
        archivo_salida: ruta de la imagen para guardar el archivo PNG
    """
    # El modelo entrega E_B < 0 para núcleos ligados. B/A representa -E_B (positivo)
    B_A_modelo = -E_B
    
    plt.figure(figsize=(8, 5))
    
    # Graficar la curva del modelo
    plt.plot(A, B_A_modelo, '-', color='#162A54', label='Modelo (Bethe-Weizsäcker)', linewidth=2)
    
    # Graficar los puntos experimentales
    A_exp = [item[2] for item in experimentales]
    B_A_exp = [item[3] for item in experimentales]
    plt.scatter(A_exp, B_A_exp, color='#B02A20', marker='o', s=50, zorder=5, label='Datos Experimentales (AME2020)')
    
    # Etiquetas para cada punto experimental para facilitar el análisis
    for item in experimentales:
        plt.annotate(item[0], (item[2], item[3]), textcoords="offset points", xytext=(5,5), ha='left', fontsize=8)
        
    plt.title('Energía de enlace por nucleón ($B/A$) vs Número de masa ($A$)')
    plt.xlabel('$A$ (Número de nucleones)')
    plt.ylabel('$B/A$ (Energía de enlace por nucleón) [MeV]')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    
    # Guardar gráfica
    plt.savefig(archivo_salida, dpi=300)
    plt.close()
