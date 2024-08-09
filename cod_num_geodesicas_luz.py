import matplotlib.pyplot as plt
import numpy as np
import skimage # para transformações, mudar o código lá em baixo
from PIL import Image

# resoluções da tela, em pixels, resolution_height x resolution_width = quantidade de curvas calculadas

resolution_height = 32
resolution_width = 32

# dimensões da janela observacional

window_height = 0.00001 #h_P
window_width = (resolution_width / resolution_height) * window_height #w_P
distance_from_window = -1.4e-4 #d_P

# coordenadas para mapeamento
# sem disco coords_no_aDisk = (theta, phi, boolean)
coords_no_aDisk = np.zeros((resolution_height, resolution_width, 3))

# com disco
# com disco coords_no_aDisk = (raio, phi, boolean)
coords_aDisk = np.zeros((resolution_height,resolution_width, 3))

# passo do runge-kunta 4

stepsize = 0.1

# dados e equações relevantes pra métrica

G = 1 #constante gravitacional
M = 1 #massa do buraco negro
a = 0.99 #momento angular por massa a = J (Mc)^-1, J é o momento angular, também chamado de "parâmetro de Kerr"
c = 1 #velocidade da luz
R = 2*G*M / c**2 #raio de Schwarzschild
radius_celestial_sphere = 80 #raio da esfera celeste
aDiskMin = 2*R #limite inferior do disco de acreção
aDiskMax = 5*R #limite superior

# definindo pi para facilitar a vida

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
     for i in range(resolution_height): # (i,j) é a localização de um pixel na imagem dada
     # localização dos pixels na janela observacional
          if (j<(resolution_width / 2 - 0.01 * resolution_width)):
               h = (window_height / 2) - (i - 1) * window_height / (resolution_height - 1)
               w = (-window_width / 2) + (j - 1) * window_width / (resolution_width - 1)
          else:
               h = (window_height / 2) - (i - 1) * window_height / (resolution_height - 1)
               w = (-window_width / 2) + (0.01 * window_width) + (j - 1) * window_width / (resolution_width - 1)


          # condições iniciais
          r = 79.999
          theta = pi / 2 - pi / 40
          phi = pi / 2
          t_dot = 1 
          phi_dot = ((csc(theta) * w) / (sqrt((a**2 + r**2) * (distance_from_window**2 + w**2 + h**2))))

          # momento conjugado  
          p_r = ((2 * Sigma(r,theta) * (h * (a**2 + r**2) * cos(theta) + r * sqrt(a**2 + r**2) * sin(theta) * distance_from_window))
                / (sqrt(distance_from_window**2 + h**2 + w**2) * (a**2 + 2 * r**2 + a**2 * cos(2 * theta)) * Delta(r)))
          p_theta = ((2 * Sigma(r,theta) * (-h * r * sin(theta) + sqrt(a**2 + r**2) * cos(theta) * distance_from_window))
                  / (sqrt(distance_from_window**2 + h**2 + w**2) * (a**2 + 2 * r**2 + a**2 * cos(2 * theta))))
          # quantias conservadas
          E = (1 - R / r) * t_dot + (R * a * phi_dot) / r # energia do vetor de killing com respeito ao tempo t
          L = -(R * a) / r * t_dot + (r**2 + a**2 + (R * a**2) / r) * phi_dot # momento angular do vetor de killing com respeito a phi
             
          # geodésicas, sistema de EDOs de primeira ordem, batem com o programa do Rodrigo
             
          def f(r, theta, p_r, p_theta):
               f1 = (p_r * Delta(r)) / Sigma(r, theta) #
               f2 = (p_theta) / Sigma(r, theta) 
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
          
          # define o tipo de plot
          plot_type = 0

          match plot_type:
               case 1:
                    while R < curve[k, 0] and curve[k, 0] < radius_celestial_sphere and k < Nk:
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
                         curve[k] = x_att          
                    curve = curve[:k]
                    cart = np.zeros((3, len(curve)))
                    for i in range(len(curve)):
                         cart[:, i] = Boyer2Cart(curve[i, 0], curve[i, 1], curve[i, 2])
                    ax.plot(cart[0], cart[1], cart[2])
               
               case 0:
                    passed_through_ADisk = 0

                    while 1.2 * R < curve[0, 0] and curve[0, 0] < radius_celestial_sphere and k < Nk:
                         # runge-kutta
                         temp_curve_step = min([stepsize*Delta(curve[0, 0]), stepsize])
                         k1 = temp_curve_step * f(r, theta, p_r, p_theta)
                         k2 = temp_curve_step * f(r + 0.5 * k1[0], theta + 0.5 * k1[1], p_r + 0.5 * k1[3], p_theta + 0.5 * k1[4])
                         k3 = temp_curve_step * f(r + 0.5 * k2[0], theta + 0.5 * k2[1], p_r + 0.5 * k2[3], p_theta + 0.5 * k2[4])
                         k4 = temp_curve_step * f(r + k3[0], theta + k3[1], p_r + k3[3], p_theta + k3[4])
                         
                         r = r + (k1[0] + (2 * k2[0]) + (2 * k3[0]) + k4[0]) / 6
                         theta = theta + (k1[1] + (2 * k2[1]) + (2 * k3[1]) + k4[1]) / 6
                         phi = phi + (k1[2] + (2 * k2[2]) + (2 * k3[2]) + k4[2]) / 6
                         p_r = p_r + (k1[3] + (2 * k2[3]) + (2 * k3[3]) + k4[3]) / 6
                         p_theta = p_theta + (k1[4] + (2 * k2[4]) + (2 * k3[4]) + k4[4]) / 6
                         
                         k = k + 1
                         
                         x_att = np.array([r, theta, phi, p_r, p_theta])               
                         curve[k] = x_att

                         if passed_through_ADisk == 0:
                              # checa se o step temporário passa pelo disco de acreção
                              if (curve[k, 1] - pi / 2) * (curve[k, 1] - pi / 2) < 0:
                                   if aDiskMin < curve[k, 0] and curve[k, 0] < aDiskMax:
                                        coords_aDisk[i,j,:] = [curve[k, 0], curve[k, 2], 1]
                                        passed_through_ADisk = 1
                         
                         # arrumando os valores das coordenadas
                         curve[k, 1] = curve[k, 1] % (2 * pi)
                         curve[k, 2] = curve[k, 2] % (2 * pi)
                         if curve[k, 1] > pi:
                              curve[k, 1] = 2 * pi - curve[k, 1]
                              curve[k, 2] = (pi + curve[k, 2]) % (2 * pi)
                         
                         k = k + 1
                    
                         if (1.2 * R) < curve[0, 0]:
                              coords_no_aDisk[i, j, :] = np.array([curve[0, 1], curve[0, 2], 0])
                              print(coords_no_aDisk)
                         else:
                              coords_no_aDisk[i, j, :] = np.array([curve[0, 1], curve[0, 2], 1])
                              print(coords_no_aDisk)
# configura a proporção dos eixos para serem iguais
print("acabei de rodar!")

# salvando as coordenadas num arquivo 
np.savez('imgMap_no_BN.npz',
          coords_no_aDisk = coords_no_aDisk,
          coords_aDisk = coords_aDisk,
          aDiskMin = aDiskMin,
          aDiskMax = aDiskMax)

# gerando imagens de diferentes cenas celestiais

image_type_w_disk = 1
boolean_aDisk = 1

# loadando os dados salvos
dados = np.load('imgMap_no_BN.npz')
coords_no_aDisk = dados['coords_no_aDisk']
coords_aDisk = dados['coords_aDisk']
aDiskMin = dados['aDiskMin']
aDiskMax = dados['aDiskMax']

# dimensões da imagem gerada
resolution_height, resolution_width, _ = coords_no_aDisk.shape
print(coords_no_aDisk.shape)

# lendo a imagem sem BN

accretion_disk_scene = Image.open('adisk.jpg')
accretion_disk_scene_array = np.array(accretion_disk_scene)

# dando as dimensões da cena com o disco de acreção

accretion_disk_scene_heigth_res, accretion_disk_scene_width_res, _ = accretion_disk_scene_array.shape

print(accretion_disk_scene_array.shape)
print(accretion_disk_scene_heigth_res)
print(accretion_disk_scene_width_res)

# mesma coisa pra imagem da simulação
celestial_scene = Image.open('kevin_gill_tycho_2_star_map.jpg')
celestial_scene_array = np.array(celestial_scene)
celestial_scene_height_res, celestial_scene_width_res, _ = celestial_scene_array.shape

print(celestial_scene_array.shape)
print(celestial_scene_height_res)
print(celestial_scene_width_res)

# iniciando a imagem gerada como uma matriz de zeros, totalmente preta
IMG = np.zeros((resolution_height, resolution_width, 3), dtype = np.uint8)

print(IMG.shape) 
print(IMG.dtype)

batch_size = 100
for j in range(resolution_width):
     for i in range(resolution_height):
          if image_type_w_disk == 1:
               if boolean_aDisk == 1:
                    if coords_aDisk[i, j, 2] == 1:

                         # limpar os valores de coordenadas
                         coords_aDisk[i, j, 1] = coords_aDisk[i, j, 1] % [2 * pi]
                         if coords_aDisk[i, j, 1] > pi:
                              coords_aDisk[i, j, 1] = 2 * pi - coords_aDisk[i, j, 1]
                         
                         if coords_no_aDisk[i, j, 2] == 1:
                              IMG[i, j, :] = accretion_disk_scene_array[round((coords_aDisk[i,j,0] - aDiskMin) * ((accretion_disk_scene_heigth_res - 1) / (aDiskMax - aDiskMin)) + 1),
                                                                        round(coords_aDisk[i, j, 1] * ((accretion_disk_scene_width_res - 1) / (2 * pi)) + 1), :]
                         
                         else:
                              IMG[i, j, :] = max(accretion_disk_scene_array[round((coords_aDisk[i, j, 0] - aDiskMin) * ((accretion_disk_scene_heigth_res - 1) / (aDiskMax - aDiskMin)) + 1),
                                                                            round(coords_aDisk[i, j, 1] * ((accretion_disk_scene_width_res - 1) / (2 * pi)) + 1), :],
                                                 celestial_scene_array[round(coords_no_aDisk[i, j, 0] * ((celestial_scene_height_res - 1) / pi) + 1),
                                                                        round(coords_no_aDisk[i, j, 1] * ((celestial_scene_width_res - 1) / (2 * pi)) + 1), :])
                         
                         boolean_aDisk = 0
          
          if boolean_aDisk == 1:
               if coords_no_aDisk[i, j, 2] == 1:
                    IMG[i, j, 0] = 0
                    IMG[i, j, 1] = 0
                    IMG[i, j, 2] = 0
               
               else:
                    IMG[i, j, :] = celestial_scene_array[round(coords_no_aDisk[i, j, 0] * ((celestial_scene_height_res - 1) / pi) + 1),
                                                         round(coords_no_aDisk[i, j, 1] * ((celestial_scene_width_res - 1) / (2 * pi)) + 1), :]
          
          boolean_aDisk = 1

print(f'dimensão da imagem: {IMG.shape}')
# manipulando a imagem

IMG = np.concatenate([
     IMG[:, :min(int(resolution_width / 2 - 3), IMG.shape[1] - 1), :],
     IMG[:, min(int(resolution_width / 2 - 3), IMG.shape[1] - 1):min(int(resolution_width / 2 - 3), IMG.shape[1] - 1) + 1, :],
     np.maximum(IMG[:, min(int(resolution_width / 2 - 3), IMG.shape[1] - 1):min(int(resolution_width / 2 - 3), IMG.shape[1] - 1) + 1, :],
                IMG[:, min(int(resolution_width / 2 + 5), IMG.shape[1] - 1):min(int(resolution_width / 2 + 5), IMG.shape[1] - 1) + 1, :]),
     IMG[:, min(int(resolution_width / 2 + 5), IMG.shape[1] - 1):min(int(resolution_width / 2 + 5), IMG.shape[1] - 1) + 1, :],
     IMG[:, min(int(resolution_width / 2 + 5), IMG.shape[1] - 1):, :]
], axis=1)

img_pil = Image.fromarray(IMG)

img_pil.save('imagem.jpg', format='JPEG')