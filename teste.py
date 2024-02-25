import matplotlib.pyplot as plt
import numpy as np

R = 2
# plotagem
# cria uma figura 3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

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

ax.set_box_aspect([1,1,1])

plt.show()