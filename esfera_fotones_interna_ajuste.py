import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipkinc#, ellipkm1
import math

def P_real(R, a, b):

    b2 = b * b
    # P^3 - b^2 * N + b^2

    # p0 * P^3 + p1 * P^2 + p2 * P + p3
    J = b2 * (R - 1.0) / (R - a)
    n = 4
    p = np.zeros(n)
    p[0] = 1.0
    p[1] = 0.0
    p[2] = - J
    p[3] = a * J

    x = np.roots(p)
    x_value = np.isreal(x) #return a true if x is not complex

    y = []
    for i in range(0, n - 1):
        if x_value[i]:
            A = np.real(x[i])
            if A >= 1.5 * a: # P >= 3 / 2 * a
                y = np.append(y, [A], 0)
	
    return y

def N_complejo(R, a, b):

    b2 = b * b
    # P^3 - b^2 * N + b^2

    # p0 * P^3 + p1 * P^2 + p2 * P + p3
    J = b2
    n = 4
    p = np.zeros(n)
    p[0] = 1.0
    p[1] = 0.0
    p[2] = - J
    p[3] = - J

    x = np.roots(p)

    x_value = np.isreal(x) #return a true if x is not complex

    y = []
    for i in range(0, n - 1):
        if x_value[i]:
            A = np.real(x[i])
            if A <= 3: # N <= 3 * a
                y = np.append(y, [A], 0)
	
    return y

def I_menos(R, a, b):

	P = P_real(R, a, b)
	#print(P)
	P = float(P)
	Q = np.sqrt((P - a) * (P + 3 * a))
	k = np.sqrt((Q - P + 3 * a) / (2 * Q))
	mu = 2.0 * np.sqrt(P * a / Q)
	M = 2 * a * (P / R - 1)
	B = Q - P + 3 * a
	S = np.sqrt(1 + (M / B))
	XR = np. arcsin(S)
	I = mu * (ellipk(k * k) - ellipkinc(XR, k * k)) / np.sqrt(a)

	return I

def I_mas(R, b):

	N = N_complejo(R, 1, b)
	#print(P)
	N = float(N)
	a1 = N + 3
	a2 = np.sqrt(2 * N + 3)
	k = np.sqrt(0.5 + 0.25 * a1 / a2)
	mu = np.sqrt(N / a2)
	TgR = np.sqrt((1 + N / R) / a2)
	XR = 2.0 * np. arctan(TgR)
	TgR = np.sqrt(1 / a2)
	Xinf = 2.0 * np. arctan(TgR)
	I = mu * (ellipkinc(XR, k * k) - ellipkinc(Xinf, k * k))

	return I



def carga_datos(R12, R):
	#str_R12 = str('%.2f' % R12)
	str_R = str('%.1f' % R)
	s = f'_R12_{R12}_R_{R}'
	ubicacion = '/Users/pedroschmied/Desktop/Defensa/R_' + str_R
	archivo = ubicacion + '/plano_thin-shell' + s + '.txt'
	
	M = np.loadtxt(archivo)
	b = M[:, 0]
	end = M[:, 1]
	phi = M[:, 2]
	r_transferencia = M[:, 3:len(M[0,:])]
	
	return b, end, phi, r_transferencia

def deflexion_analitica(R12, R, lim_xi, lim_xf):
	bc = 1.5 * np.sqrt(3)
	bc_interno = bc * R12 * np.sqrt((R - R12) / (R - 1))
	#############################################
	b, end, phi, r_transferencia = carga_datos(R12, R)
	laps = phi / (2.0 * np.pi)
	axes.grid()	
	axes.plot(b, laps, 'r', lw = 1.5)#, label = r'$\phi$')
	#############################################
	#############curvas "analíticas"#############
	
	#analitico:
	#x1 = np.arange( bc - 0.2, bc - 0.00001, 0.001)
	x2 = np.arange(bc_interno + 0.0001, bc - 0.0001, 0.005)

	#j01 = np.zeros(len(x1))
	j02 = np.zeros(len(x2))
	a = R12
	
	# para bc^-
	for i in range(0, len(x2)):
		I_menos_b = 0#I_menos(R, R12, x2[i])
		I_mas_b = I_mas(R, x2[i])

		j02[i] = 2 * (I_menos_b + I_mas_b) / (2.0 * np.pi)
	'''
	for i in range(0, len(x1)):
		I_menos_b = I_complejo(R, a, x1[i], 1 / a, 1 / R)
		I_mas_b = I_complejo(R, 1, x1[i], 1 / R, 0)
		j01[i] = (I_menos_b + I_mas_b - 0. * np.pi) / (2.0 * np.pi)
	'''
	
	#axes.plot(x1, j01, 'g', lw = 1.8)
	axes.plot(x2, j02, 'g', lw = 1.8)
	
	#plt.ylim(0, 1)
	axes.set_xlabel('b', fontsize = 12, weight = 'bold')
	axes.set_ylabel('n', fontsize = 12, weight = 'bold')#Δϕ


	axes.set_xlim(lim_xi, lim_xf)
	plt.ylim(0.3, 3.5)
	footer_text = f'R ~ {round(R, 3)}'
	plt.text(.2, 0.65, footer_text, horizontalalignment='center', verticalalignment='center', transform = axes.transAxes,size=9, weight='semibold', color = 'k')
	footer_text = f'R_ ~ {round(R12, 3)}'
	plt.text(.2, 0.58, footer_text, horizontalalignment='center', verticalalignment='center', transform = axes.transAxes,size=9, weight='semibold', color = 'k')

def deflexion(R12, R, lim_xi, lim_xf):
	bc = 1.5 * np.sqrt(3)
	bc_interno = bc * R12 * np.sqrt((R - R12) / (R - 1))
	##############################################
	b, end, phi, r_transferencia = carga_datos(R12, R)
	laps = phi / (2.0 * np.pi)
	axes.grid()	
	axes.plot(b, laps, 'r', lw = 1.5)
	#############################################
	#############curvas "analíticas"#############
	B = bc_interno
	x1 = np.arange(B - 0.2, B - 0.00001, 0.001)
	x2 = np.arange(B + 0.0001, B + 0.3, 0.001)
	j1 = np.zeros(len(x1))
	j2 = np.zeros(len(x2))

	cociente_mayor = R / R12
	C1 = np.log(18**2 * np.sqrt(3) * R12 * (7 - 4 * np.sqrt(3)) * np.sqrt((R - R12) / (R - 1.0)))
	C1 = C1 - np.log((2 * cociente_mayor + 1.5 + np.sqrt(3 * cociente_mayor * (3 + cociente_mayor))) / (cociente_mayor - 1.5))
	cociente = R12 / R
	C2 = np.log(36**2 * np.sqrt(3) * R12 *  np.sqrt((R - R12) / (R - 1.0))) + 2.0 * np.log(1 - 1.5 * cociente) - 4.0 * np.log(np.sqrt(3) + np.sqrt(1 + 3 * cociente))

	I_mas_b = I_mas(R, bc_interno)
	print(I_mas_b)
	for i in range(0, len(x2)):
		j2[i] = (C2 - np.log(x2[i] - bc_interno) + 2 * I_mas_b) / (2.0 * np.pi)

	for i in range(0, len(x1)):
		j1[i] = (C1 - np.log(bc_interno - x1[i]) + I_mas_b) / (2.0 * np.pi)

	axes.plot(x1, j1, '--k', lw = 1.8)
	axes.plot(x2, j2, '--k', lw = 1.8)#, label = r'$\phi$')
	
	#plt.ylim(0, 1)
	axes.set_xlabel('b', fontsize = 12, weight = 'bold')
	axes.set_ylabel('n', fontsize = 12, weight = 'bold')#Δϕ
	#plt.title(f'Vueltas', fontsize = 14, weight = 'bold')
	#axes[row, col].set_xlim(0, 4)
	#plt.legend()
	axes.set_yticks(np.arange(0.25, 3.5, 0.5))
	axes.set_xlim(lim_xi, lim_xf)
	plt.ylim(0.3, 3)
	footer_text = f'R ~ {round(R, 3)}'
	plt.text(.2, 0.65, footer_text, horizontalalignment='center', verticalalignment='center', transform = axes.transAxes,size=9, weight='semibold', color = 'k')
	footer_text = f'R_ ~ {round(R12, 3)}'
	plt.text(.2, 0.58, footer_text, horizontalalignment='center', verticalalignment='center', transform = axes.transAxes,size=9, weight='semibold', color = 'k')


R = 1.3
# 27 imagenes
# 2 partes: 1 de 3 x 5 y otra de 3 x 4
parte = 2

nrow = 1
ncol = 1
fig = plt.figure(figsize = (5, 4))#3 * ncol, 2.2 * nrow))

dR12 = 0.02
R12 = 0.6
bc = np.sqrt(3) * 1.5
bc_interno = bc * R12 * np.sqrt((R - R12) / (R - 1))

lim_xi = bc_interno - 0.3
lim_xf = bc_interno + 0.3


for i in range(0, nrow * ncol):
    if R12 == 1:
        R12 = int(1)
    axes = fig.add_subplot(nrow, ncol, i + 1)
    deflexion(R12, R,lim_xi, lim_xf)
    #deflexion_analitica(R12, R, lim_xi, lim_xf)
    R12 = '%0.3f'%(R12 + dR12)
    R12 = float(R12)


plt.subplots_adjust(left = 0.15, right = 0.97, bottom = 0.15, top = 0.97)
plt.show()
