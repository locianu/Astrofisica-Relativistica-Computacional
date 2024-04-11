import matplotlib.pyplot as plt
import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D


#dados e equações relevantes pra métrica

G = 1 #constante gravitacional
M = 1 #massa do buraco negro
a = 0.6 #momento angular por massa a = J (Mc)^-1, J é o momento angular, também chamado de "parâmetro de Kerr"
c = 1 #velocidade da luz
R = 2*G*M / c**2 #raio de Schwarzschild
radius_celestial_sphere = 80 #raio da esfera celeste
aDiskMin = 2*R #limite inferior do disco de acreção
aDiskMAx = 5*R #limite superior

#definindo pi para facilitar a vida

pi = math.pi

# funções matemáticas

def sin(theta):
     return np.sin(theta)

def cos(theta):
     return np.cos(theta)

def csc(theta):
     return 1 / np.sin(theta)

def sqrt(theta):
     return np.sqrt(theta)

def tan(theta):
     if theta != pi / 2:
          return np.tan(theta)
     else:
          print("Invalid vector")

def cot(theta):
     if theta != 0:
          return 1 / np.tan(theta)
     else:
          print("Invalid vector")

# funções da métrica

def Sigma(r,theta):
    return r**2 + (a**2) * (cos(theta))**2

def dSigmadr(r):
    return 2 * r

def dSigmadtheta(theta):
    return -2 * (a**2) * cos(theta) * sin(theta)

def Delta(r):
    return (r**2) - R * r + (a**2)

def dDeltadr(r):
    return 2 * r - R

# conversão de coordenadas

def Boyer2Cart(r, theta, phi):
     x = sqrt(r**2 + a**2) * sin(theta) * cos(phi)
     y = sqrt(r**2 + a**2) * sin(theta) * sin(phi)
     z = r * cos(theta)
     return np.column_stack((x, y, z))


# resoluções da tela, em pixels

resolution_height = 4
resolution_width = 4

# dimensões da janela observacional

window_height = 0.00001 #h_P
window_width = (resolution_width / resolution_height) * window_height #w_P
distance_from_window = -1.4e-4 #d_P

# coordenadas para mapeamento
coords_no_aDisk = np.zeros((resolution_height, resolution_width, 3))
coords_aDisk = np.zeros((resolution_height,resolution_width, 3))

# passo do runge-kunta 4

stepsize = 0.1

# lista vazia para adicionar as curvas

all_curves = []

for j in range(resolution_width):
        for i in range(resolution_height): #(i,j) é a localização de um pixel na imagem dada
          #localização dos pixels na janela observacional
             if (j<(resolution_width / 2 - 0.01 * resolution_width)):
               h = (window_height / 2) - (i - 1) * window_height / (resolution_height - 1)
               w = (-window_width / 2) + (j - 1) * window_width / (resolution_width - 1)
             else:
               h = (window_height / 2) - (i - 1) * window_height / (resolution_height - 1)
               w = (-window_width / 2) + (0.01 * window_width) + (j - 1) * window_width / (resolution_width - 1)


             #posição inicial da "câmera" e por consequinte das geodésicas (0, r(0), theta(0), phi(0)) 
             #em coordenadas Boyer-Lindquist
             r = 70
             theta = pi / 2 - pi / 46 # deslocando o buraco negro para visualização com o aDisk (disco de acreção?)
             phi = 0
             t_dot = 1 
             phi_dot = ((csc(theta) * w) / (sqrt((a**2 + r**2) * (distance_from_window**2 + w**2 + h**2))))

             #momento conjugado  
             p_r = ((2 * Sigma(r,theta) * (h * (a**2 + r**2) * cos(theta) + r * sqrt(a**2 + r**2) * sin(theta) * distance_from_window))
                    / (sqrt(distance_from_window**2 + h**2 + w**2) * (a**2 + 2 * r**2 + a**2 * cos(2 * theta)) * Delta(r)))
             p_theta = ((2 * Sigma(r,theta) * (-h * r * sin(theta) + sqrt(a**2 + r**2) * cos(theta) * distance_from_window))
                    / (sqrt(distance_from_window**2 + h**2 + w**2) * (a**2 + 2 * r**2 + a**2 * cos(2 * theta))))
             
             curve = np.array([[r, theta, phi, p_r, p_theta]])

             #quantias conservadas
             E = (1 - R / r) * t_dot + (R * a * phi_dot) / r #energia do vetor de killing com respeito ao tempo t
             L = (-(R * a) / r * t_dot + (r**2 + a**2 + (R * a**2) / r) * phi_dot) #momento angular do vetor de killing com respeito a phi
             
             # geodésicas, sistema de EDOs de primeira ordem
             
             def geodesic():
               def f(r, theta, p_r, p_theta):
                    f1 = (p_r * Delta(r)) / Sigma(r, theta)
                    f2 = (p_theta) / Sigma(r, theta)
                    f3 = ((a * ( -a * L + r * R * E) + L * csc(theta)**2 * Delta(r)) / (Delta(r) * Sigma(r,theta)))
                    f4 = (- (1 / (2 * Delta(r)**2 * Sigma(r,theta)**2)) * (Sigma(r,theta) 
                         * (-E * Delta(r) * (a * R * ( -2 * L + a * E * sin(theta)**2) 
                         + 2 * r * E * Sigma(r,theta)) + (a * (a * L**2 - 2 * L * r * R * E + a * r * R * E**2 
                         * sin(theta)**2) + p_r ** 2 * Delta(r)**2 + (a**2 + r**2) * E**2 * Sigma(r,theta)) 
                         * dDeltadr(r)) + Delta(r) * (a * (L * (a * L - 2 * r * R * E) + a * r * R * E**2 * sin(theta)**2) 
                         - Delta(r) * (p_theta**2 + L**2 * csc(theta)**2 + p_r**2 * Delta(r))) * dSigmadr(r)))
                    f5 = (- (1 / (2 * Delta(r) * Sigma(r,theta)**2)) * (-2 * sin(theta) * (a**2 * r * R * E**2 * cos(theta)
                         + L**2 * cot(theta) * csc(theta)**3 * Delta(r)) * Sigma(r,theta) + (a * (L * (a * L - 2 * r * R * E) 
                         + a * r * R * E**2 * sin(theta)**2) - Delta(r) * (p_theta**2 + L**2 * csc(theta)**2 + p_r**2 * Delta(r))) * dSigmadtheta(theta)))
                    return [f1,f2,f3,f4,f5]
               return f

             f = geodesic()
             
             k = 0
             Nk = 20000

             while((R < curve[k][0]) and (curve[k][0] < radius_celestial_sphere) and (k < Nk)):
                  
                  # valores de coodenadas "limpos"
                  
                  curve[k][1] = curve[k][1] % (2 * pi)
                  curve[k][2] = curve[k][2] % (2 * pi)
                  
                  if curve[k][1] > pi:
                       curve[k][1] = 2 * pi - curve[k][1]
                       curve[k][2] = (pi + curve[k][2]) % (2 * pi)
                  
                  theta = curve[k][1]

                  # runge-kutta
                  step = min([stepsize*Delta(r), stepsize])
                                    
                  k1 = step * f(r, theta, p_r, p_theta)[0]
                  m1 = step * f(r, theta, p_r, p_theta)[1]
                  n1 = step * f(r, theta, p_r, p_theta)[2]
                  s1 = step * f(r, theta, p_r, p_theta)[3]
                  v1 = step * f(r, theta, p_r, p_theta)[4]

                  k2 = step * f(r + k1 / 2, theta + m1 / 2, p_r + n1 / 2, p_theta + s1 / 2)[0]
                  m2 = step * f(r + k1 / 2, theta + m1 / 2, p_r + n1 / 2, p_theta + s1 / 2)[1]
                  n2 = step * f(r + k1 / 2, theta + m1 / 2, p_r + n1 / 2, p_theta + s1 / 2)[2]
                  s2 = step * f(r + k1 / 2, theta + m1 / 2, p_r + n1 / 2, p_theta + s1 / 2)[3]
                  v2 = step * f(r + k1 / 2, theta + m1 / 2, p_r + n1 / 2, p_theta + s1 / 2)[4]

                  k3 = step * f(r + k2 / 2, theta + m2 / 2, p_r + n2 / 2, p_theta + s2 / 2)[0]
                  m3 = step * f(r + k2 / 2, theta + m2 / 2, p_r + n2 / 2, p_theta + s2 / 2)[1]
                  n3 = step * f(r + k2 / 2, theta + m2 / 2, p_r + n2 / 2, p_theta + s2 / 2)[2]
                  s3 = step * f(r + k2 / 2, theta + m2 / 2, p_r + n2 / 2, p_theta + s2 / 2)[3]
                  v3 = step * f(r + k2 / 2, theta + m2 / 2, p_r + n2 / 2, p_theta + s2 / 2)[4]

                  k4 = step * f(r + k3, theta + m3, p_r + n3, p_theta + s3)[0]
                  m4 = step * f(r + k3, theta + m3, p_r + n3, p_theta + s3)[1]
                  n4 = step * f(r + k3, theta + m3, p_r + n3, p_theta + s3)[2]
                  s4 = step * f(r + k3, theta + m3, p_r + n3, p_theta + s3)[3]
                  v4 = step * f(r + k3, theta + m3, p_r + n3, p_theta + s3)[4]

                  r = r + (k1 + (2 * k2) + (2 * k3) + k4) / 6
                  theta = theta + (m1 + (2 * m2) + (2 * m3) + m4) / 6
                  p_r = p_r + (n1 + (2 * n2) + ( 2 * n3) + n4) / 6
                  p_theta = p_theta + (s1 + (2 * s2) + (2 * s3) + s4) / 6
                  phi = phi + (v1 + (2 * v2) + (2 * v3) + v4) / 6
                  
                  k = k + 1

                  x_att = np.array([[r,theta,phi,p_r,p_theta]])
                  curve = np.vstack((curve, x_att))

             # transforma em coordenadas cartesianas  
             cart = Boyer2Cart(curve[:, 0], curve[:, 1], curve[:, 2])

             # adiciona na lista
             all_curves.append(cart)
# plotagem
# cria uma figura 3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# mostra a figura
for curve in all_curves:
     ax.plot3D(curve[:, 0], curve[:, 1], curve[:, 2])

# cria uma grade de pontos em uma esfera

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
u, v = np.meshgrid(u, v)

x = R * np.cos(u) * np.sin(v)
y = R * np.sin(u) * np.sin(v)
z = R * np.cos(v)
# plota a esfera

ax.plot_surface(x, y, z, color='black', alpha = 1)

# configura a proporção dos eixos para serem iguais

plt.axis('equal')
ax.set_box_aspect([1,1,1])


plt.show()

