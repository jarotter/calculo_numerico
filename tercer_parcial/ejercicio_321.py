# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 17:16:50 2017

@author: SaulAlvarez
"""

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import pdes

def b(x):
    return np.cos(2.0*np.pi*x)

def l(t):
    return np.exp(-4.0*t)

def r(t):
    return 0

def u(x,t):
    return np.exp(-4*t)*np.cos(2.0*np.pi*x)

# Forward method
M=20
N=100

X=np.linspace(0,1,M+1);
T=np.linspace(0,1,N+1);
U=np.array([[u(x,t) for t in T] for x in X])
W=pdes.m_heat_exp(b,l,r,1.0/(np.pi**2),(0,1),(0,1),M,N)
print('Error en (3/10,1) con método forward:',abs(U[3*M//10][N]-W[3*M//10][N]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X, T = np.meshgrid(X, T)
ax.plot_wireframe(X,T,np.transpose(U))
plt.show()
