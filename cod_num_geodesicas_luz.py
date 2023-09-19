import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


# condições iniciais

G = 1
M = 1 #massa do buraco negro
a = 0.6 #momento angular
R = 2*G*M 
radius_celestial_sphere = 80
aDiskMin = 2*R
aDiskMAx = 5*R

# funções da métrica

def Sigma(r,theta):
    return r**2 + a**2 * (cos(theta))**2

def dSigmadr(r):
    return 2 * r

def dSigmadtheta(theta):
    return -2 * a**2 * cos(theta) * sin(theta)

def Delta(r):
    return r**2 - R * r + a**2

def dDeltadr(r):
    return 2 * r - R

# plotagem
# cria uma figura 3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# cria uma grade de pontos em uma esfera

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = R * np.outer(np.cos(u), np.sin(v))
y = R * np.outer(np.sin(u), np.sin(v))
z = R * np.outer(np.ones(np.size(u)), np.cos(v))

# plota a esfera

ax.plot_surface(x, y, z, color='black')

# configura a proporção dos eixos para serem iguais

ax.set_box_aspect([1, 1, 1])

# mostra a figura

#plt.show()

# resoluções da tela, em pixels

resolution_height = 2
resolution_width = 2

# dimensões da janela observacional, não sei se é aplicável em python

window_height = 0.00001
window_width = (resolution_width / resolution_height) * window_height
distance_from_window = -1.4e-4

coords_no_aDisk = np.zeros((resolution_height, resolution_width, 3))
coords_aDisk = np.zeros((resolution_height,resolution_width, 3))

# passo do runge-kunta

#stepsize = 0.1
h = 0.1

# sinceramente eu não entendi essa parte do programa

for j in range(1,2):
        for i in range(1,2):
             h = window_height / 2 -(i - 1) * window_height / (resolution_height - 1)
             w = -window_width / 2 + (j - 1) * window_width / (resolution_width - 1)

             r = 70
             theta = math.pi             