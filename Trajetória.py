from numpy import array, arange
from math import cos, sin, pi,sqrt
from pylab import plot, ylabel, xlabel, show, legend


g = 9.8

#equações de movimento para método de RK4

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

def grafico(R, ml):

	#parâmetros para método de RK4

	m1 = ml[0]
	m2= ml[1]
	l1 = ml[2]
	l2 = ml[3]
	N = 60000
	a, b = 0, 60
	h = (b - a)/N
	T = arange(a, b, h)

	T1, W1, T2, W2 = [], [], [], []

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
	return(T1,T2)


def movimento(R, ml):

	d0 = 1e-12 #distância entre condições iniciais (pode ser definida com um valor maior para visualizar mais caos em menos tempo de simulação do RK4)
	eps = array([d0/sqrt(2), 0.0, d0/sqrt(2),0.0])
	Reps = R + eps #segundo pêndulo duplo
	l1 = ml[2]
	l2 = ml[3]

	T11, T12 = grafico(R,ml)
	T21, T22 = grafico(Reps, ml)
	X1, Y1, X2, Y2 = [],[],[],[]


	for i in range(len(T11)): #posições dos pêndulos
		X1.append(l1*sin(T11[i]) + l2*sin(T12[i]))
		Y1.append(-l1*cos(T11[i]) - l2*cos(T12[i]))
		X2.append(l1*sin(T21[i]) + l2*sin(T22[i]))
		Y2.append(-l1*cos(T21[i]) - l2*cos(T22[i])) 

	#plot das posições

	plot(X1,Y1, 'r', label = 'R')
	plot(X2, Y2, label = 'R + $d_0$')
	xlabel('x')
	ylabel('y')
	legend()
	show()

R1 = array([pi/2, 0.0, pi/2, 0.0], float) #theta1, omega1, theta2, omega2
ml = array([1.0, 10.0, 1.0, 1.0], float)   #m1, m2, l1, l2

movimento(R1, ml)