# Simulação problema restrito
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from funcoes_orbitais import *
from integradores import *



plt.rcParams['figure.dpi'] = 300
golden_ratio = (1 + np.sqrt(5)) / 2
figure_width = 3
figure_height = figure_width / golden_ratio
plt.rcParams['figure.figsize'] = [figure_width, figure_height]


plt.plot(trajetoria[0:k, 0], trajetoria[0:k, 1])
plt.xlabel("x (UA)")
plt.ylabel("y (UA)")
plt.title("Trajetória projetada (plano xy)")
plt.gca().set_aspect('equal', adjustable='box')
plt.show()


