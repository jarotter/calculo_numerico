#-*- coding: utf-8 -*-
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

#Este va en pdes
def m_heat_CN(b,l,r,D,Ix,It,M,N):

    h=1.0*(Ix[1]-Ix[0])/M; #Delta x
    k=1.0*(It[1]-It[0])/N; #Delta t

    sigma=1.0*D*k/(h*h);

    m=M-1; n=N;

    #Matriz tridiagonal A
    A = np.diag([2 + 2*sigma for i in range(m)])
    A = A + np.diag([-sigma for i in range(m-1)],1)
    A = A + np.diag([-sigma for i in range(m-1)],-1)

    #Aproximaciones para la condición de frontera de Neuman
    A[m-1,m-1] = 2 + (2.0/3)*sigma;
    A[m-1,m-2] = -(2.0/3)*sigma

    #Matriz tridiagonal B
    B= np.diag([2 - 2*sigma for i in range(m)])
    B = B + np.diag([sigma for i in range(m-1)],1)
    B = B + np.diag([sigma for i in range(m-1)],-1)

    #Aproximaciones para la condición de frontera de Neuman
    B[m-1,m-1] = 2 -(2.0/3)*sigma;
    B[m-1,m-2] = (2.0/3)*sigma

    left_side= [l(It[0]+k*i) for i in range(n+1)];

    temp_w = np.zeros((n+1,m));
    temp_w[0,:]= np.transpose( [ b(Ix[0] + h*(i+1)) for i in range(m)])

    v = np.transpose(np.zeros(m))

    for j in range (n):
        v[0]=sigma *(left_side[j]+left_side[j+1]);
        v[m-1] = 2 * h *sigma* (r(j)+r(j-1))/3.0

        aux = np.array(B.dot(np.transpose(temp_w[j,:])))
        aux = aux + v

        temp_w[j+1,:]= np.linalg.solve(A,aux)

    W = np.zeros((n+1,M+1))
    W[:,0] = np.transpose(left_side)

    for i in range(1,m+1):
        W[:,i] = temp_w[:,i-1];

    for i in range(0,n+1):
        W[i,m+1]=(4.0/3)*W[i,m] - (1.0/3)*W[i,m-1]

    return W


def b(x):
    return np.cos(2.0*np.pi*x)

def l(t):
    return np.exp(-4.0*t)

def r(t):
    return 0

#Para este problema en particular
def error(b,l,r,D,Ix,It,M,N,W):
    arr = abs(np.linspace(Ix[0], Ix[1], M) - 0.3);
    _min = min(arr)
    j = [i for i in range(len(arr)) if _min == arr[i]][0]

    valor_exacto = b(3.0/10)*l(1.0)
    aproximacion = W[N,j]
    return valor_exacto, aproximacion


D=1.0/((np.pi)*(np.pi));
Ix = [0,1]
It = [0,1]

###
print('\n Inciso a)\n')
# CN con  ∆x=1/20 y ∆t=1/10 para realizar la gráfica
W = m_heat_CN(b,l,r,D,Ix,It,20,10)

# Realizando la gráfica
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X=np.linspace(0,1,20+1);
T=np.linspace(0,1,10+1);
X, T = np.meshgrid(X, T)
ax.plot_wireframe(X, T, W)
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u(x,y)')
plt.show()

###
print('\n Inciso b)\n')
W = m_heat_CN(b,l,r,D,Ix,It,10,10)
exacto, aprox = error(b,l,r,D,Ix,It,10,10,W)
print('h=1/10, k=1/10')
print('valor exacto: ', exacto)
print('aproximación: ', aprox)
e1 = abs(exacto - aprox)
print('error absoluto: ', e1 )

W= m_heat_CN(b,l,r,D,Ix,It,20,20)
exacto, aprox = error(b,l,r,D,Ix,It,20,20,W)
print('\nh=1/20, k=1/20')
print('valor exacto: ', exacto)
print('aproximación: ', aprox)
e2 = abs(exacto - aprox)
print('error absoluto: ', e2 )
print ('eoc = ', abs(np.log(e2/e1)/np.log(2)))

W = m_heat_CN(b,l,r,D,Ix,It,40,40)
exacto, aprox = error(b,l,r,D,Ix,It,40,40,W)
print('\nh=1/40, k=1/40')
print('valor exacto: ', exacto)
print('aproximación: ', aprox)
e3 = abs(exacto - aprox)
print('error absoluto: ', e3 )
print ('eoc = ', abs(np.log(e3/e2)/np.log(2)))

W = m_heat_CN(b,l,r,D,Ix,It,80,80)
exacto, aprox = error(b,l,r,D,Ix,It,80,80,W)
print('\nh=1/80, k=1/80')
print('valor exacto: ', exacto)
print('aproximación: ', aprox)
e4 = abs(exacto - aprox)
print('error absoluto: ', e4 )
print ('eoc = ', abs(np.log(e4/e3)/np.log(2)))


###
print('\n Inciso c)\n')
W = m_heat_CN(b,l,r,D,Ix,It,20,25)
exacto, aprox = error(b,l,r,D,Ix,It,20,25,W)
print('h=1/20, k=1/25')
print('valor exacto: ', exacto)
print('aproximación: ', aprox)
print('error absoluto: ', abs(exacto - aprox))

W = m_heat_CN(b,l,r,D,Ix,It,20,50)
exacto, aprox = error(b,l,r,D,Ix,It,20,50,W)
print('\nh=1/20, k=1/50')
print('valor exacto: ', exacto)
print('aproximación: ', aprox)
print('error absoluto: ', abs(exacto - aprox))

W = m_heat_CN(b,l,r,D,Ix,It,20,100)
exacto, aprox = error(b,l,r,D,Ix,It,20,100,W)
print('\nh=1/20, k=1/100')
print('valor exacto: ', exacto)
print('aproximación: ', aprox)
print('error absoluto: ', abs(exacto - aprox))

W = m_heat_CN(b,l,r,D,Ix,It,20,200)
exacto, aprox = error(b,l,r,D,Ix,It,20,200,W)
print('\n h=1/20, k=1/200')
print('valor exacto: ', exacto)
print('aproximación: ', aprox)
print('error absoluto: ', abs(exacto - aprox))
