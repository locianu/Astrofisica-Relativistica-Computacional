import matplotlib.pyplot as plt
import numpy as np


# resoluções da tela, em pixels, resolution_height x resolution_width = quantidade de curvas calculadas

resolution_height = 2
resolution_width = 2

# dimensões da janela observacional

window_height = 0.00001 #h_P
window_width = (resolution_width / resolution_height) * window_height #w_P
distance_from_window = -1.5e-4 #d_P

# coordenadas para mapeamento
coords_no_aDisk = np.zeros((resolution_height, resolution_width, 3))
coords_aDisk = np.zeros((resolution_height,resolution_width, 3))

# passo do runge-kunta 4

stepsize = 0.1

#dados e equações relevantes pra métrica

G = 1 #constante gravitacional
M = 1 #massa do buraco negro
a = 0.99 #momento angular por massa a = J (Mc)^-1, J é o momento angular, também chamado de "parâmetro de Kerr"
c = 1 #velocidade da luz
R = 2*G*M / c**2 #raio de Schwarzschild
radius_celestial_sphere = 80 #raio da esfera celeste
aDiskMin = 2*R #limite inferior do disco de acreção
aDiskMAx = 5*R #limite superior

#definindo pi para facilitar a vida

pi = np.pi

# cria uma figura 3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# cria uma grade de pontos em uma esfera

u = np.linspace(0, 2 * pi, 100)
v = np.linspace(0, pi, 100)

x = R * np.outer(np.cos(u), np.sin(v))
y = R * np.outer(np.sin(u), np.sin(v))
z = R * np.outer(np.ones(np.size(u)), np.cos(v))
# plota a esfera

ax.plot_surface(x, y, z, rstride=4, cstride=4, color='black')


# funções matemáticas

def sin(theta):
     return np.sin(theta)

def cos(theta):
     return np.cos(theta)

def csc(theta):
     if theta == 0:
          return 0
     else:
          return 1 / np.sin(theta)

def sqrt(theta):
     return np.sqrt(theta)

def tan(theta):
     if theta == pi / 2:
          return 0
     else:
          return np.tan(theta)

def cot(theta):
     if theta == 0:
          return 0
     else:
          return 1 / np.tan(theta)

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
     A = a
     x = sqrt(r**2 + A**2) * sin(theta) * cos(phi)
     y = sqrt(r**2 + A**2) * sin(theta) * sin(phi)
     z = r * cos(theta)
     return np.array([x, y, z])

for j in range(resolution_width):
     for i in range(resolution_height): #(i,j) é a localização de um pixel na imagem dada
     #localização dos pixels na janela observacional
          if (j<(resolution_width / 2 - 0.01 * resolution_width)):
               h = (window_height / 2) - (i - 1) * window_height / (resolution_height - 1)
               w = (-window_width / 2) + (j - 1) * window_width / (resolution_width - 1)
          else:
               h = (window_height / 2) - (i - 1) * window_height / (resolution_height - 1)
               w = (-window_width / 2) + (0.01 * window_width) + (j - 1) * window_width / (resolution_width - 1)


          # condições iniciais
          r = 70
          theta = pi / 2 - pi / 46 # deslocando o buraco negro para visualização com o aDisk (disco de acreção?)
          phi = pi / 2
          t_dot = 1 
          phi_dot = ((csc(theta) * w) / (sqrt((a**2 + r**2) * (distance_from_window**2 + w**2 + h**2))))

          # momento conjugado  
          p_r = ((2 * Sigma(r,theta) * (h * (a**2 + r**2) * cos(theta) + r * sqrt(a**2 + r**2) * sin(theta) * distance_from_window))
                / (sqrt(distance_from_window**2 + h**2 + w**2) * (a**2 + 2 * r**2 + a**2 * cos(2 * theta)) * Delta(r)))
          p_theta = ((2 * Sigma(r,theta) * (-h * r * sin(theta) + sqrt(a**2 + r**2) * cos(theta) * distance_from_window))
                  / (sqrt(distance_from_window**2 + h**2 + w**2) * (a**2 + 2 * r**2 + a**2 * cos(2 * theta))))
          # quantias conservadas
          E = (1 - R / r) * t_dot + (R * a * phi_dot) / r #energia do vetor de killing com respeito ao tempo t
          L = -(R * a) / r * t_dot + (r**2 + a**2 + (R * a**2) / r) * phi_dot #momento angular do vetor de killing com respeito a phi
             
          # geodésicas, sistema de EDOs de primeira ordem, batem com o programa do Rodrigo
             
          def f(r, theta, p_r, p_theta):
               f1 = (p_r * Delta(r)) / Sigma(r, theta) #ok
               f2 = (p_theta) / Sigma(r, theta) #o
               f3 = ((a * ( -a * L + r * R * E) + (L * csc(theta)**2) * Delta(r)) / (Delta(r) * Sigma(r,theta))) #ok
               f4 = (- (1 / ((2 * Delta(r)**2) * Sigma(r,theta)**2)) * (Sigma(r,theta) 
                    * (-E * Delta(r) * (a * R * ( -2 * L + a * E * sin(theta)**2) 
                    + 2 * r * E * Sigma(r,theta)) + (a * (a * L**2 - 2 * L * r * R * E + a * r * R * (E**2) 
                    * sin(theta)**2) + (p_r**2) * Delta(r)**2 + (a**2 + r**2) * (E**2) * Sigma(r,theta)) 
                    * dDeltadr(r)) + Delta(r) * (a * (L * (a * L - 2 * r * R * E) + a * r * R * (E**2) * sin(theta)**2) 
                    - Delta(r) * (p_theta**2 + (L**2) * (csc(theta)**2) + (p_r**2) * Delta(r))) * dSigmadr(r))) 
               f5 = (- (1 / (2 * Delta(r) * Sigma(r,theta)**2)) * (-2 * sin(theta) * ((a**2) * r * R * (E**2) * cos(theta)
                    + (L**2) * cot(theta) * (csc(theta)**3) * Delta(r)) * Sigma(r,theta) + (a * (L * (a * L - 2 * r * R * E) 
                    + a * r * R * (E**2) * sin(theta)**2) - Delta(r) * (p_theta**2 + (L**2) * csc(theta)**2 + (p_r**2) * Delta(r))) * dSigmadtheta(theta))) 
               return np.array([f1,f2,f3,f4,f5])
             
          k = 0
          Nk = 20000
          x0 = np.array([[r, theta, phi, p_r, p_theta]])
          curve = np.zeros((Nk, 5))
          curve[0] = x0
          print(curve[0])

          while R < curve[k, 0] and curve[k, 0] < radius_celestial_sphere and k < Nk-1:
               curve[k, 1] = curve[k, 1] % (2 * pi)
               curve[k, 2] = curve[k, 2] % (2 * pi)
                  
               if curve[k, 1] > pi:
                    curve[k, 1] = 2 * pi - curve[k, 1]
                    curve[k, 2] = (pi + curve[k, 2]) % (2 * pi)
                  
               # runge-kutta
               step = min([stepsize*Delta(curve[k, 0]), stepsize])
               k1 = step * f(r, theta, p_r, p_theta)
               k2 = step * f(r + 0.5 * k1[0], theta + 0.5 * k1[1], p_r + 0.5 * k1[3], p_theta + 0.5 * k1[4])
               k3 = step * f(r + 0.5 * k2[0], theta + 0.5 * k2[1], p_r + 0.5 * k2[3], p_theta + 0.5 * k2[4])
               k4 = step * f(r + k3[0], theta + k3[1], p_r + k3[3], p_theta + k3[4])
               
               r = r + (k1[0] + (2 * k2[0]) + (2 * k3[0]) + k4[0]) / 6
               theta = theta + (k1[1] + (2 * k2[1]) + (2 * k3[1]) + k4[1]) / 6
               phi = phi + (k1[2] + (2 * k2[2]) + (2 * k3[2]) + k4[2]) / 6
               p_r = p_r + (k1[3] + (2 * k2[3]) + (2 * k3[3]) + k4[3]) / 6
               p_theta = p_theta + (k1[4] + (2 * k2[4]) + (2 * k3[4]) + k4[4]) / 6
               
               k = k + 1
               
               x_att = np.array([r, theta, phi, p_r, p_theta])
               #print(f'x_att = {x_att}')
               curve[k] = x_att
               if k == 1:
                    print(curve[1])
          curve = curve[:k]
          cart = np.zeros((3, len(curve)))
          for i in range(len(curve)):
               cart[:, i] = Boyer2Cart(curve[i, 0], curve[i, 1], curve[i, 2])
          ax.plot(cart[0], cart[1], cart[2])
# configura a proporção dos eixos para serem iguais

plt.axis('square')

# plotagem
plt.show()

