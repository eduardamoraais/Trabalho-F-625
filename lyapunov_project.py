from numpy import array, arange, max, mean, log, sum, empty, logspace, linspace
from numpy.linalg import norm
from math import cos, sin, pi, sqrt
from pylab import plot, ylabel, xlabel, show, imshow, colorbar, loglog, savefig, figure, subplots, title
from tqdm import tqdm

g = 9.8
mass_ratio = [2*10**i for i in linspace(-2, 2, 100)] #razões entre massa a serem usadas
length_ratio = [2*10**i for i in linspace(-2, 2, 100)] #razões entre comprimentos dos fios a serem usadas
m1, l1 = 1.0, 1.0 #massa e comprimento do fio do pêndulo interno

#parâmetros para método de RK4

N = 6000
a, b = 0, 30
tau = b
h = (b - a)/N
T = arange(a, b, h)
T_lamb = arange(25, b, h)

d0 = 1e-12 #distância entre condições iniciais

R1 = array([pi/2, 0.0, pi/2, 0.0], float) #theta1, omega1, theta2, omega2
R2 = array([pi/2 - d0/sqrt(2), 0.0, pi/2 - d0/sqrt(2), 0.0], float)

n = len(mass_ratio)
LYAPUNOV = empty([n, n], float) #inicialização da matriz dos expoentes a serem plotados no gráfico de densidade

#equações do movimento para método RK4

def f(r, t, m1, m2, l1, l2):
  t1 = r[0]
  w1 = r[1]
  t2 = r[2]
  w2 = r[3]
  ft1 = w1
  fw1 = -(g*(2*m1 + m2)*sin(t1) + m2*g*sin(t1 - 2*t2) + 2*sin(t1 - t2)*m2*(w2**2*l2 + w1**2*l1*cos(t1 - t2)))/(l1*(2*m1 + m2 - m2*cos(2*(t1 - t2))))
  ft2 = w2
  fw2 = (2*sin(t1 - t2)*(w1**2*l1*(m1 + m2) + g*(m1 + m2)*cos(t1) + w2**2*l2*m2*cos(t1 - t2)))/(l2*(2*m1 + m2 - m2*cos(2*(t1 - t2))))
  return array([ft1, fw1, ft2, fw2], float)

#método de RK4

def RK4(R, m2, l2):
  k1 = h*f(R, t, m1, m2, l1, l2)
  k2 = h*f(R + 0.5*k1, t + 0.5*h, m1, m2, l1, l2)
  k3 = h*f(R + 0.5*k2, t + 0.5*h, m1, m2, l1, l2)
  k4 = h*f(R + k3, t + h, m1, m2, l1, l2)
  R += (k1 + 2*k2 + 2*k3 + k4)/6
  return R

#atualização do segundo pêndulo a ser feita a cada iteração, a fim de manter válida a aproximação dos expoentes

def R2_update(R1, R2, d):
  t11, t12 = R1[0], R2[0]
  w11, w12 = R1[1], R2[1]
  t21, t22 = R1[2], R2[2]
  w21, w22 = R1[3], R2[3]

  dt1 = (d0/d)*(t12 - t11)
  dw1 = (d0/d)*(w12 - w11)
  dt2 = (d0/d)*(t22 - t21)
  dw2 = (d0/d)*(w22 - w21)

  R2 = array([t11 + dt1, w11 + dw1, t21 + dt2, w21 + dw2], float)

  return R2

#loop para calcular expoentes para cada ângulo inicial de interesse

name = int(1)
for angle in tqdm([pi/40, pi/6, pi/4, pi/3, pi/2], desc = 'Angle loop'):
  R1 = array([angle, 0.0, angle, 0.0], float) #theta1, omega1, theta2, omega2
  R2 = array([angle - d0/sqrt(2), 0.0, angle - d0/sqrt(2), 0.0], float)
  for i in tqdm(range(n), desc = 'Mass loop', leave = False):
    for j in tqdm(range(n), desc = 'Length loop', leave = False):

      mr = mass_ratio[i]
      lr = length_ratio[j]
      m2, l2 = m1*mr, l1*lr #massa e comprimento do fio do pêndulo externo
      LAMB = []

      for t in T: #um passo de RK4 é realizado a cada iteração dentro do tempo de simulação
        R1 = RK4(R1, m2, l2)
        R2 = RK4(R2, m2, l2)

        d = norm(R1 - R2)
        LAMB.append(log(d/d0)) #cálculo do expoente
        
        R2 = R2_update(R1, R2, d) #atualização do segundo pêndulo duplo

      LAMBf = []

      for k in range(int(6000*(25/30)), len(T)): #média sobre os expoentes na região em que já se estailizaram
        LAMBf.append(sum(LAMB[:k + 1])/((k + 1)*h)) #armazenamento das médias

      lyapunov = mean(LAMBf) #média sobre todos os expoentes no intervalo considerado

      LYAPUNOV[i, j] = lyapunov

      #reinicialização das posições iniciais para testar novos parâmetros

      R1 = array([angle, 0.0, angle, 0.0], float) #theta1, omega1, theta2, omega2
      R2 = array([angle - d0/sqrt(2), 0.0, angle - d0/sqrt(2), 0.0], float)

  lmax = max(LYAPUNOV) #apenas caso se queira definir um vmax na escala do gráfico de densidade, proporcional ao máximo valor dos expoentes

  #plot do gráfico de densidade dos expoentes

  fig, ax = subplots()
  ax = imshow(LYAPUNOV, origin = 'lower', extent = [-2, 2, -2, 2])
  xlabel('$\log_{10}(L_2/L_1)$')
  ylabel('$\log_{10}(M_2/M_1)$')
  title(f'$Angle = {round(angle/pi, 3)} \pi$')
  cb = fig.colorbar(ax)
  cb.set_label('$\lambda$')
  savefig('lyapunov' + str(name))
  name += int(1)