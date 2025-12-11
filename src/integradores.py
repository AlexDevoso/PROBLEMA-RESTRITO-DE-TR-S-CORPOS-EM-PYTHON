# Integradores numéricos
# Simulação problema restrito
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from funcoes_orbitais import *

nMax = 10000  # número máximo de passos
dt = 30 / dias_por_ano           # passo de tempo em anos (30 dias)
times = np.arange(0, dt * nMax, dt)



trajetoria = np.zeros((nMax, 6))
trajetoria[0] = Y0
k = 1

t_eval = np.arange(0, 2 * dt, dt)


energia_vals = np.zeros(nMax)
jacobi_vals = np.zeros(nMax)
jacobi_vals[0] = constante_jacobi(Y0)
df_energia = np.zeros(nMax) 

while k < nMax:
    
    sol = solve_ivp(derivadas_estado, method="DOP853", t_span=[0, 2*dt],
        y0=trajetoria[k-1, :], t_eval=t_eval)
    ri = sol.y[:, -1]
    trajetoria[k, :] = ri
    
    jacobi_vals[k] = constante_jacobi(trajetoria[k])
    df_energia[k] = np.abs((jacobi_vals[k] - jacobi_vals[0]) / jacobi_vals[0])
    k += 1