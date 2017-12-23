# -*- coding: utf-8 -*-
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
import pdes


def u0(x):
    y =  x - 2 * np.floor((x+1)/2)
    return np.where(abs(y)<=1/2, 1-2*abs(y), 0)

def u(x,t):
    return u0(x-2*t)


c = 2
Ix = np.array([-1, 1])
It = np.array([0, 10])
M = 40

#LAS LÍNEAS 22-42 SON LA ANIMACIÓN
# N = 400
#
# W = pdes.lax_friedrich(u0, c, Ix, It, M, N)
#
# grid_space = np.linspace(Ix[0], Ix[1], M+1)
# delta_time = (It[1] - It[0])/N;
# min_y = min( u(grid_space, 0) ) - 0.1
# max_y = max( u(grid_space, 0) ) + 0.1
#
# for i in range(N+1):
#     plt.plot(grid_space, u(grid_space, i*delta_time), 'k',
#         grid_space, W[:,i], 'bo');
#     plt.axis([Ix[0], Ix[1], min_y, max_y]);
#     plt.show(block = False)
#     plt.pause(0.00000001)
#     plt.close()
#
# plt.plot(grid_space, u(grid_space, N*delta_time), 'k',
#           grid_space, W[:,N], 'bo');
# plt.axis([Ix[0], Ix[1], min_y, max_y]);
# plt.show(block = True)

#Partes dos y tres
print()
x = np.linspace(Ix[0], Ix[1], M+1)

t =  np.zeros(5)
f = np.zeros(5)
h_max = np.zeros(5)
x_max = np.zeros(5)
volumen = np.zeros(5)

N = [500, 400]
for n in N:
    print('Para sigma = {}'.format(400/n))
    W = pdes.lax_friedrich(u0, c, Ix, It, M, n)
    for i in range(5):
        t[i] = int(2*(i+1))
        f = W[:,int(n/5)*(i+1)]
        h_max[i] = max(f)
        x_max[i] = Ix[0] + np.argmax(f)*1/20
        volumen[i] = simps(f,x)

    print('t = 2n |     Altura máxima   | x de altura máxima |     "Volumen"')
    print('--------------------------------------------------------------------')
    for i in range(4):
        print(t[i],'   |','  {0:.13f}'.format(h_max[i]),'  |        ','{}'.format(x_max[i]),'      |','{0:.13f}'.format(volumen[i]))
    print(t[-1],'  |','  {0:.13f}'.format(h_max[-1]),'  |        ','{0:0.3}'.format(x_max[-1]),'      |','{0:.13f}'.format(volumen[-1]))
    print()

#Parte 5
N = 200
i = 0
W = pdes.lax_friedrich(u0, c, Ix, It, M, N)

while True:
    if max(abs(W[:,i])) > 5:
        break
    i = i + 1
print()

print('El primer momento en que se excede 5 es t = {}'.format(It[0]+1/20*i))
print()
