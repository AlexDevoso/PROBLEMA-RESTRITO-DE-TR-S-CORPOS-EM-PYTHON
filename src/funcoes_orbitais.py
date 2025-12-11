# Funções orbitais
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --------------------------
G = 6.67430E-11                # N m^2 / kg^2
ua_m = 149597870700            # 1 UA em metros
seg_por_dia = 86400           # segundos por dia
dias_por_ano = 365.25         # dias por ano

# parâmetros gravitacionais
GMsun = 1.32712440041279419E20
GMneptune = 6836527.100580E9

#constantes para unidades (UA e anos)
mu_sol = GMsun * ((dias_por_ano * seg_por_dia)**2 / ua_m**3)      # mu1 em UA^3 / ano^2
mu_netuno = GMneptune * ((dias_por_ano * seg_por_dia)**2 / ua_m**3)  # mu2 em UA^3 / ano^2

#azão conservada do problema reduzido (Dúvidas)
mu = GMneptune / (GMsun + GMneptune)

# velocidade angular do referencial rotativo (ano^-1)
omega = 2 * np.pi / 166.53457869719753

# distância característica (aphelium) usada para posicionamento inicial dos primários
r_aphelio = ((mu_sol + mu_netuno) / (omega**2))**(1/3)
# posições dos dois corpos primários em x (no referencial rotativo)
x_sol = -mu * r_aphelio
x_netuno = (1 - mu) * r_aphelio

# condição inicial da partícula (unidades: UA e anos)
# --------------------------
x0, y0, z0 = 18.955301453285813, -28.505415095903587, -0.8022897227312524
u0, v0, w0 = -0.05371381524210492, -0.3168520373583842, -0.33413526140611555

R0 = np.array((x0, y0, z0))
V0 = np.array((u0, v0, w0))
Y0 = np.hstack((R0, V0))   # vetor de estado [x,y,z,xdot,ydot,zdot]


# função que retorna as derivadas (Dúvidas sobre a derivadas )
def derivadas_estado(t, Y):
  x, y, z = Y[:3]
  xdot, ydot, zdot = Y[3:]
  Ydot = np.zeros_like(Y)
  # primeiras três derivadas são as velocidades
  Ydot[:3] = Y[3:]
  distancia_ao_sol = np.sqrt((x - x_sol)**2 + y**2 + z**2)
  distancia_a_netuno = np.sqrt((x - x_netuno)**2 + y**2 + z**2)
  Ydot[3] = (
      2 * ydot * omega
      + x * omega**2
      - mu_sol * (x - x_sol) / distancia_ao_sol**3
      - mu_netuno * (x - x_netuno) / distancia_a_netuno**3
  )
  Ydot[4] = (
      -2 * xdot * omega
      + y * omega**2
      - mu_sol * y / distancia_ao_sol**3
      - mu_netuno * y / distancia_a_netuno**3
  )
  Ydot[5] = z * (- mu_sol / distancia_ao_sol**3 - mu_netuno / distancia_a_netuno**3)
  return Ydot

#função para definir a constante de jacobi
def constante_jacobi(Y):
    x, y, z = Y[:3]
    xdot, ydot, zdot = Y[3:]

    r1 = np.sqrt((x - x_sol)**2 + y**2 + z**2)
    r2 = np.sqrt((x - x_netuno)**2 + y**2 + z**2)

    # potencial efetivo U
    U = 0.5 * omega**2 * (x**2 + y**2) + mu_sol / r1 + mu_netuno / r2

    v2 = xdot**2 + ydot**2 + zdot**2

    C = 2 * U - v2
    return C