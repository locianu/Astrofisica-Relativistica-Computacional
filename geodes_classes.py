import numpy as np
#dados e equações relevantes pra métrica

G = 1 #constante gravitacional
M = 1 #massa do buraco negro
a = 0.6 #momento angular por massa a = J (Mc)^-1, J é o momento angular, também chamado de "parâmetro de Kerr"
c = 1 #velocidade da luz
R = 2*G*M / c**2 #raio de Schwarzschild
radius_celestial_sphere = 80 #raio da esfera celeste
aDiskMin = 2*R #limite inferior do disco de acreção
aDiskMAx = 5*R #limite superior

class metrica:
    def __init__(self, radius, theta, phi):
        self.radius = radius
        self.theta = theta
        self.phi = phi
    
    # funções matemáticas

    def sin(self):
        return np.sin(self.theta)

    def cos(self):
        return np.cos(self.theta)

    def csc(self):
        return 1 / np.sin(self.theta)

    def sqrt(self):
        return np.sqrt(self.theta)

    def tan(self):
        return np.tan(self.theta)
    
    def cot(self):
        return 1 / np.tan(self.theta)

    # funções da métrica

    def Sigma(self):
        return self.radius**2 + a**2 * (cos(self.theta))**2

    def dSigmadr(self):
        return 2 * self.radius

    def dSigmadtheta(self):
        return -2 * a**2 * cos(self.theta) * sin(self.theta)

    def Delta(self):
        return self.radius**2 - R * self.radius + a**2

    def dDeltadr(self):
        return 2 * self.radius - R

    # conversão de coordenadas

    def Boyer2Cart(self):
        x = sqrt(self.radius**2 + a**2) * sin(self.theta) * cos(self.phi)
        y = sqrt(self.radius**2 + a**2) * sin(self.theta) * sin(self.phi)
        z = r * cos(self.theta)
        return np.column_stack((x, y, z))

    