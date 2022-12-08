from vpython import cylinder, vector, sphere, rate
from numpy import array, arange
from math import cos, sin, pi
from pylab import plot, ylabel, xlabel, show

g = 9.8
m1, m2, l1, l2 = 1.0, 10.0, 1.0, 1.0 #massas e comprimentos dos fios

R = array([pi/2, 0.0, pi/2, 0.0], float) #theta1, omega1, theta2, omega2

#parâmetros para método RK4

N = 100000
a, b = 0, 100
h = (b - a)/N
T = arange(a, b, h)

T1, W1, T2, W2, Elist = [], [], [], [], []

#equações de movimento para RK4

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

for t in T:
  T1.append(R[0])
  W1.append(R[1])
  T2.append(R[2])
  W2.append(R[3])
  k1 = h*f(R, t, m1, m2, l1, l2)
  k2 = h*f(R + 0.5*k1, t + 0.5*h, m1, m2, l1, l2)
  k3 = h*f(R + 0.5*k2, t + 0.5*h, m1, m2, l1, l2)
  k4 = h*f(R + k3, t + h, m1, m2, l1, l2)
  R += (k1 + 2*k2 + 2*k3 + k4)/6

#inicialização dos elementos da animação (esferas e fios)

s1 = sphere(pos = vector(l1, 0, 0), radius = 0.09*(m1)**(1/3))
s2 = sphere(pos = vector(l1 + l2, 0, 0), radius = 0.09*(m2)**(1/3))
c1 = cylinder(pos = vector(0, 0, 0), axis = vector(l1, 0, 0), radius = 0.008)
c2 = cylinder(pos = vector(l1, 0, 0), axis = vector(l1 + l2, 0, 0), radius = 0.008)

#animação do movimento

fps = 50
di = int((1/fps)/(100/len(T)))

for i, t in enumerate(T[::di]):
    rate(fps)
    x1 = l1*sin(T1[i*di])
    y1 = -l1*cos(T1[i*di])
    x2 = x1 + l2*sin(T2[i*di])
    y2 = y1 - l2*cos(T2[i*di])

    s1.pos = vector(x1, y1, 0)
    s2.pos = vector(x2, y2, 0)
    c1.pos, c1.axis = vector(0, 0, 0), vector(x1, y1, 0)
    c2.pos, c2.axis = vector(x1, y1, 0), vector(x2 - x1, y2 - y1, 0)