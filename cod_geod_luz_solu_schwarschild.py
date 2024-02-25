import numpy as np
import einsteinpy as eins
import math


c = 1 # velocidade da luz
G = 1 # constante gravitacional
M = 1 # massa do buraco negro
R = 2*G*M/c**2

# conversão de coordenadas

def Boyer2Cart(r, theta, phi):
     x = sqrt(r**2 + a**2) * sin(theta) * cos(phi)
     y = sqrt(r**2 + a**2) * sin(theta) * sin(phi)
     z = r * cos(theta)
     return np.column_stack((x, y, z))

# resoluções da tela, em pixels

resolution_height = 4
resolution_width = 4

# dimensões da janela observacional, não sei se é aplicável em python

window_height = 0.00001 #h_P
window_width = (resolution_width / resolution_height) * window_height #w_P
distance_from_window = -1.4e-4 #d_P

# stepsize runge-kutta

stepsize = 0.1

# pi

pi = math.pi

# lista vazia pra guardar as curvas

curvas = np.array([])

for j in range(resolution_width):
    for i in range(resolution_height):
        #(i,j) é a localização de um pixel na imagem dada
        #localização dos pixels na janela observacional
        h = window_height / 2 - (i - 1) * window_height / (resolution_height - 1)
        w = -window_width / 2 + (j - 1) * window_width / (resolution_width - 1)

        def geodesicas:
            def f(r, t_dot, v, omega, phi_dot, theta):
                f1 = (R / (R - r)) * t_dot * v
                f2 = ((R - r) * R / (2 * r * (r - R)**2)) * v**2 - ((R (R - r) / (2 * r**3))) * t_dot**2 - 