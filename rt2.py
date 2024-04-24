#Ray tracing
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import csv
ptype=1
rh=2
rw=2
#Dados metrica
G=1.
M=1.
#Parametro de Kerr
a=0.97
R=2*G*M
rcs=80.
aDiskMin=2*R
aDiskMax=5*R
#PLOTAR ESFERA REPRESENTANDO BURACO NEGRO
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = R * np.outer(np.cos(u), np.sin(v))
y = R * np.outer(np.sin(u), np.sin(v))
z = R * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, rstride=4, cstride=4, color='black')
#Formulas
def Sigma(r,theta):
    Sigma=r**2+(a**2)*math.cos(theta)**2
    return Sigma
def drSigma(r):
    drSigma=2*r
    return drSigma
def dthetaSigma(theta):
    dthetaSigma=-2*(a**2)*math.cos(theta)*math.sin(theta)
    return dthetaSigma
def Delta(r):
    Delta=(r**2)-R*r+(a**2)
    return Delta
def drDelta(r):
    drDelta=2*r-R
    return drDelta
#COORDENADAS JANELA OBSERVACAO
wh=0.00001
ww=(rw/rh)*wh
dfw=-2*0.00007
#COORDENADAS PARA MAPEAMENTO
cDisk=np.zeros((rh,rw,3))
cnDisk=np.zeros((rh,rw,3))
#passo de Runge-Kutta
p=0.1
j=0
i=0
while j<rw:
 j=j+1
 i=0
 while i<rh:
  i=i+1
  #COLOCAR IF NA SINGULARIDADE PI/2
  sh=0.01
  if j<(rw/2-sh*rw):
   h=wh/2-(i-1)*wh/(rh-1)
   w=-ww/2+(j-1)*ww/(rw-1)
  else:
   h=wh/2-(i-1)*wh/(rh-1)
   w=-ww/2+(sh*ww)+(j-1)*ww/(rw-1)
  #Condicoes iniciais
  r=70.
  theta=((math.pi)/2.)-((math.pi)/46.)
  phi=0.
  t_dot=1.
  phi_dot=(w/math.sin(theta))/math.sqrt((a**2+r**2)*(dfw**2+w**2+h**2))
  p_r=2*Sigma(r,theta)*(h*(a**2+r**2)*math.cos(theta)+r*math.sqrt(a**2+r**2)\
*math.sin(theta)*dfw)/(math.sqrt(dfw**2+h**2+w**2)*(a**2+2*r**2+(a**2)\
*math.cos(2*theta))*Delta(r))
  p_theta=2*Sigma(r,theta)*(-h*r*math.sin(theta)+math.sqrt(a**2+r**2)\
*math.cos(theta)*dfw)/(math.sqrt(dfw**2+h**2+w**2)*(a**2+2*r**2+(a**2)\
*math.cos(2*theta)))
  #Quantidades conservadas
  E=(1-R/r)*t_dot+(R*a*phi_dot)/r
  L=-(R*a)/r*t_dot+(r**2+a**2+(R*a**2)/r)*phi_dot
  #Equacoes da Geodesica
  #x=[r,theta,phi,pr,ptheta]
  #dx=[dr,dtheta,dphi,dpr,dptheta]
  def f(x1,x2,x4,x5): 
      f=np.array([(x4*Delta(x1))/Sigma(x1,x2),
      x5/Sigma(x1,x2),
      (a*(-a*L+x1*R*E)+(L/math.sin(x2)**2)*Delta(x1))/(Delta(x1)*Sigma(x1,x2)),
      -(1/((2*Delta(x1)**2)*Sigma(x1,x2)**2))*(Sigma(x1,x2)*(-E*Delta(x1)*(a*R*(-2*L+a*E*math.sin(x2)**2)+2*x1*E*Sigma(x1,x2))\
+(a*(a*L**2-2*L*x1*R*E+a*x1*R*(E**2)*math.sin(x2)**2)+(x4**2)*Delta(x1)**2+(a**2+x1**2)*(E**2)*Sigma(x1,x2))*drDelta(x1))+Delta(x1)*(a*(L*(a*L-2*x1*R*E)+a*x1*R*(E**2)*math.sin(x2)**2)\
-Delta(x1)*(x5**2+(L**2)*(1/math.sin(x2))**2+(x4**2)*Delta(x1)))*drSigma(x1)),
-(1/(2*Delta(x1)*Sigma(x1,x2)**2))*(-2*math.sin(x2)*((a**2)*x1*R*(E**2)*math.cos(x2)+(L**2)*(1/math.tan(x2))\
*((1/math.sin(x2))**3)*Delta(x1))*Sigma(x1,x2)+(a*(L*(a*L-2*x1*R*E)+a*x1*R*(E**2)\
*math.sin(x2)**2)-Delta(x1)*(x5**2+(L**2)*(1/math.sin(x2))**2+(x4**2)*Delta(x1)))*dthetaSigma(x2))])
      return f
  def rk4step(x1,x2,x3,x4,x5,step):
      f1=step*f(x1,x2,x4,x5)
      f2=step*f(x1+(1/2.)*f1[0],x2+(1/2.)*f1[1],x4+(1/2.)*f1[3],x5+(1/2.)*f1[4])
      f3=step*f(x1+(1/2.)*f2[0],x2+(1/2.)*f2[1],x4+(1/2.)*f2[3],x5+(1/2.)*f2[4])
      f4=step*f(x1+f3[0],x2+f3[1],x4+f3[3],x5+f3[4])
      z=np.array([x1+(1/6.)*(f1[0]+2*f2[0]+2*f3[0]+f4[0]),
                  x2+(1/6.)*(f1[1]+2*f2[1]+2*f3[1]+f4[1]),
                  x3+(1/6.)*(f1[2]+2*f2[2]+2*f3[2]+f4[2]),
                  x4+(1/6.)*(f1[3]+2*f2[3]+2*f3[3]+f4[3]),
                  x5+(1/6.)*(f1[4]+2*f2[4]+2*f3[4]+f4[4])])
      return z
  #Resolvendo por RK4
  #parametros
  #numero maximo de passos
  ktot=20000
  x0=np.array([[r],[theta],[phi],[p_r],[p_theta]])
  x0l=np.array([r,theta,phi,p_r,p_theta])
  curve=np.arange(ktot*5,dtype=float).reshape(ktot,5)
  curve[0]=x0l
  print(curve[0])
  if ptype==1:
   k=0
   while R<curve[k][0] and curve[k][0]<rcs and k<ktot-1:
    curve[k][1]=curve[k][1]%(2*math.pi)
    curve[k][2]=curve[k][2]%(2*math.pi)
    if curve[k][1]>(math.pi):
     curve[k][1]=(2*math.pi)-curve[k][1]
     curve[k][2]=(math.pi+curve[k][2])%(2*math.pi)
    k=k+1
    #ESCREVER FUNCAO RK4STEP
    rk4=rk4step(curve[k-1][0],curve[k-1][1],curve[k-1][2],curve[k-1][3],curve[k-1][4],min([p*Delta(curve[k-1][0]),p]))
    curve[k]=rk4
    if k==1:
     print(curve[1])
   print(k)
   kpep=k
   #Tranformar para coordenadas euclidianas e plot
   def boyer2cart(r,theta,phi):
    A=a
    boyer2cart=np.array([math.sqrt(r**2+A**2)*math.sin(theta)*math.cos(phi),math.sqrt(r**2+A**2)*math.sin(theta)*math.sin(phi),r*math.cos(theta)])
    return boyer2cart
   cart=np.arange((kpep+1)*3,dtype=float).reshape(3,(kpep+1))
   k=0
   while k<kpep+1:
    cart[:,k]=boyer2cart(curve[k,0],curve[k,1],curve[k,2])
    k=k+1
   print(k)
   # ~ ax.plot_wireframe(cart[0], cart[1], cart[2])
   ax.plot(cart[0], cart[1], cart[2])

"""  else:
   k=0
   ptd=0
   #Curvas computadas na esfera celestial onde ela existe
   while 1.2*R<curve[k][0] and curve[k][0]<rcs and k<ktot-1:
    tcs=rk4step(curve[k-1][0],curve[k-1][1],curve[k-1][2],curve[k-1][3],curve[k-1][4],min([p*Delta(curve[k-1][0]),p]))
    if ptd==0:
     #checar se o passo temporario passa pelo disco de acrescimo
     if (tcs[1]-math.pi/2)*(curve[1]-math.pi/2)<0:
      if aDiskMin<tcs[0] and tcs[0]<aDiskMax:
       cDisk[i,j,:]=np.array([tcs[0],tcs[2],1])
       ptd=1
    #atualizar curva com tcs
    curve[k,:]=tcs
    #limpar coordenadas
    curve[k][1]=curve[k][1]%(2*math.pi)
    curve[k][2]=curve[k][2]%(2*math.pi)
    if curve[k][1]>(math.pi):
     curve[k][1]=(2*math.pi)-curve[k][1]
     curve[k][2]=(math.pi+curve[k][2])%(2*math.pi)
    k=k+1
   #Salvar coordenadas de acordo com os casos
   if 1.2*R<curve[k,0]:
    cnDisk[i,j,:]=np.array([curve[k,1],curve[k,2],0])
   else:
    cnDisk[i,j,:]=np.array([curve[k,1],curve[k,2],1])"""
plt.axis('square')
plt.show()

#Salvar coordenadas em arquivo
with open('raytracer.csv','w') as csvFile:
    writer=csv.writer(csvFile)
    writer.writerows(cnDisk)
csvFile.close()
with open('raytracer1.csv','w') as csvFile:
    writer=csv.writer(csvFile)
    writer.writerows(cDisk)
csvFile.close()
with open('raytracer2.csv','w') as csvFile:
    writer=csv.writer(csvFile)
    writer.writerows(aDiskMin)
csvFile.close()
with open('raytracer3.csv','w') as csvFile:
    writer=csv.writer(csvFile)
    writer.writerows(aDiskMax)
csvFile.close()
#Gerar imagens finais das cenas de diferentes esferas celestiais

