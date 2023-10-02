import matplotlib.pyplot as plt
import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D


# condições iniciais

G = 1
M = 1 #massa do buraco negro
a = 0.6 #momento angular
R = 2*G*M 
radius_celestial_sphere = 80
aDiskMin = 2*R
aDiskMAx = 5*R

# funções matemáticas

def sin(theta):
     return math.sin(theta)

def cos(theta):
     return math.cos(theta)

def csc(theta):
     return 1 / math.sin(theta)

def sqrt(theta):
     return math.sqrt(theta)

def tan(theta):
     return math.tan(theta)
def cot(theta):
     return 1 / math.tan(theta)

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

plt.show()

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

#sts = 0.1
stepsize = 0.1

#definindo pi para facilitar a vida

pi = math.pi

# sinceramente eu não entendi essa parte do programa

for j in range(resolution_width):
        for i in range(resolution_height):
             h = window_height / 2 - (i - 1) * window_height / (resolution_height - 1)
             w = -window_width / 2 + (j - 1) * window_width / (resolution_width - 1)
             #h = window_height / 2 - i * window_height / (resolution_height - 1)
             #w = -window_width / 2 + j * window_width / (resolution_width - 1)

             r = 70
             theta = pi / 2 - pi / 46 # deslocando o buraco negro para visualização com o aDisk (disco de acreção?)
             phi = 0
             t_dot = 1             

             phi_dot = ((csc(theta) * w) / 
                        (sqrt((a**2 + r**2) * distance_from_window**2 + w**2 + h**2)))
             p_r = (2 * Sigma(r, theta) * (h * (a**2 + r**2) * cos(theta) + r)
                     + r * sqrt(a**2 + r**2) * sin(theta) * distance_from_window 
                     / (sqrt(distance_from_window**2 + h**2 + w**2) * (a**2 +2 
                     * r**2 + a**2 * cos(2 * theta)) * Delta(r)))
             p_theta = (2 * Sigma(r, theta) * 
                        (-h * r * sin(theta) + sin(theta) + sqrt(a**2 + r**2)
                         * cos(theta) * distance_from_window)
                         / (sqrt(distance_from_window**2 + h**2 + w**2) 
                         * (a**2 + 2*r**2 + a**2 * cos(2*theta))))

             E = (1 - R / r) * t_dot + (R * a * phi_dot) / r
             L = (-(R * a) / r * t_dot 
                  + (r**2 + a**2 + (R * a**2) / r) * phi_dot)
             
             # geodésicas, no código original elas são definidas DENTRO do loop de for, dá pra fazer o mesmo aqui? porque tanto E quanto L não são definidos até o momento em que geodesics() é definida
             
             def geodesic(): 
               def f(r, theta, p_r, p_theta):
                    global E
                    global L
                    f1 = (p_r * Delta(r)) / Sigma(r, theta)
                    f2 = (p_theta) / Sigma(r, theta)
                    f3 = (-(1 / (2 * Delta(r)**2 * Sigma(r, theta)**2)) 
                          * (Sigma(r, theta) * (-E * Delta(r) * (a * R * (-2 * L + a*E*sin(theta)**2) 
                          + 2 * r * E * Sigma(r, theta)) + (a * (a * L**2 - 2 * L * r * R * E 
                          + a * r * R * E**2 * sin(theta)**2) + p_r**2 * Delta(r)**2 +(a**2 + r**2) 
                          * E**2 * Sigma(r, theta)) * dDeltadr(r)) + Delta(r)
                          * (a * (L * (a * L - 2 * r * R * E) + a * r * R * E**2 * sin(theta)**2)
                          - Delta(r) * (p_theta**2 + L**2 * csc(theta)**2 + p_r**2 * Delta(r))) * dSigmadr(r)))
                    f4 = (-(1 / (2 * Delta(r) * Sigma(r, theta)**2)) * (-2 * sin(theta) * 
                         (a**2 * r * R * E**2 * cos(theta) + L**2 * cot(theta) * csc(theta)**3 * Delta(r)) 
                          * Sigma(r, theta) + (a * (L * (a * L - 2 * r * R * E) + a * r * R * E**2 
                          * sin(theta)**2) - Delta(r) * (p_theta**2 + L**2 * csc(theta)**2 + p_r**2 
                          * Delta(r))) * dSigmadtheta(theta)))
                    f5 = ((a * (-a * L + r * R * E) + L * csc(theta)**2 * Delta(r)) 
                          / (Delta(r) * Sigma(r, theta)))
                    return [f1, f2, f3, f4, f5]
               return f
             
             f = geodesic()
             
             x_0 = np.array([[r, theta, p_r, p_theta, phi]])
             curve = np.copy(x_0)
             
             k = 0
             Nk = 20000

             while((R<r) and (r<radius_celestial_sphere) and (k<Nk)):
                  
                  # valores de coodenadas "limpos"
                  
                  curve[k][1] = curve[k][1] % (2 * pi)
                  curve[k][4] = curve[k][4] % (2 * pi)
                  
                  if curve[k][1] > pi:
                       curve[k][1] = 2 * pi - curve[k][1]
                       curve[k][4] = (pi + curve[k][4]) % (2 * pi)
                      
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
                  print(f"raio = {r}")

                  x = np.array([r, theta, p_r, p_theta, phi])
                  # isso aqui talvez dê problema no futuro, então testar com
                  #curve[k:] = x
                  curve = np.vstack((curve, x))

             # transformando tudo pra coordenadas euclidianas e plotagem

             #n, m = curve.shape
             #A = a * np.ones((n, 1))

             # por algum motivo, tem isso aqui no código original
             #PHIZ = np.ones(n, 1)

             def Boyer2Cart(r, theta, phi):
                  x = sqrt(r**2 + a**2) * sin(theta) * cos(phi)
                  y = sqrt(r**2 + a**2) * sin(theta) * sin(phi)
                  z = r * cos(theta)
                  return np.column_stack((x, y, z))
             
             cart = Boyer2Cart(curve[:, 0], curve[:, 1], curve[:, 4])
             plt.plot(cart[:, 0], cart[:, 1], cart[:, 2])
             plt.show()

             