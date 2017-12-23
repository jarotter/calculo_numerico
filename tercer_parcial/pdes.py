# -*- coding: utf-8 -*-
import numpy as np


def lax_friedrich(f, c, Ix, It, M, N):
    """Solves the transport equation.

    Lax-Friedrich method for solving linear, hyperbolic partial differential
    equations with periodic boundary conditions.


    Parameters
    ----------
    f: function
        Initial condition
    c: double
        Transport speed
    Ix: numpy vector
        Interval of interest in space
    It: numpy vector
        Interval of interest in time
    M: int
        Number of subintervals in space
    N:
        Number of subintervals in time

    Returns
    -------
    numpy vector
       W, the solution matrix of the PDEs whose i,jth entry is the
       numerical approximation at the ith point in space and jth point in time.
    """
    h = (Ix[1]-Ix[0])/M;
    k = (It[1]-It[0])/N;
    x = np.linspace(Ix[0], Ix[1], M+1);

    W = np.zeros([M+1, N+1]);
    W[:, 0] = f(x);

    A = np.diag((1/2 + (c*k)/(2*h))*np.ones(M-2), -1) + \
        np.diag((1/2 - (c*k)/(2*h))*np.ones(M-2), 1);

    u = np.zeros(M-1);
    v = np.zeros(M-1);

#Aquí cambié a N+1
    for i in range(1, N+1):
        W[0,i] = (1/2)*(W[1,i-1]+W[M-1,i-1]) - \
                 ((c*k)/(2*h))*(W[1,i-1]-W[M-1,i-1]);
        W[M,i] = W[0,i];
        u[0] = W[0, i-1];
        v[M-2] = W[M, i-1];
        W[1:M, i] = np.dot(A,W[1:M, i-1]) + (1/2 + (c*k)/(2*h))*u + (1/2 - (c*k)/(2*h))*v;

    return W

def leapfrog(f, c, Ix, It, M, N):
    """Solves the transport equation.

    Leapfrog method for solving linear, hyperbolic partial differential
    equations with periodic boundary conditions.


    Parameters
    ----------
    f: function
        Initial condition
    c: double
        Transport speed
    Ix: numpy vector
        Interval of interest in space
    It: numpy vector
        Interval of interest in time
    M: int
        Number of subintervals in space
    N:
        Number of subintervals in time

    Returns
    -------
    numpy vector
       W, the solution matrix of the PDEs whose i,jth entry is the
       numerical approximation at the ith point in space and jth point in time.
    """

    h = (Ix[1] - Ix[0])/M;
    k = (It[1] - It[0])/N;
    x = np.linspace(Ix[0], Ix[1], M+1);

    W = np.zeros([M+1, N+1]);
    W[:,0] = f(x);

    A = np.diag(((c*k)/(h))*np.ones(M-2), -1) + np.diag(((-c*k)/(h))*np.ones(M-2), 1);
    u = np.zeros(M-1);
    v = np.zeros(M-1);

    Itprima = np.array([It[0], It[0]+2*k])

    W[:,1] = lax_friedrich(f, c, Ix, Itprima, M, 2)[:,1];

    for i in range(2, N+1):
        W[0, i] = W[0, i-2] + ((c*k)/h)*(W[M-1, i-1]-W[1, i-1]);
        W[M, i] = W[M, i-2] + ((c*k)/h)*(W[M-1, i-1] - W[1, i-1]);
        u[0] = W[0, i-1];
        v[M-2] = W[M, i-1];
        W[1:M, i] = W[1:M, i-2] + np.dot(A, W[1:M, i-1]) + ((c*k)/h)*u - ((c*k)/h)*v;

    return W

def m_heat_exp(b,l,r,D,Ix,It,M,N):
    # Solves u_t=Du_xx in Ix \times It
    #
    # b is boundary condition b(x)=u(x,t0) (function)
    # l is boundary condition l(t)=u(a,t) (function)
    # r is boundary condition r(t)=u_x(b,t) (function)
    # D is the heat equation constant (float)
    # Ix=(a,b) (numpy array)
    # It=(t_0,t_0+T) (numpy array)
    # M is the number of space intervals (int)
    # N is the number of time intervals (int)
    #
    # Returns W, a numpy array where W[i][j] represents an approximation for
    # the solution after i space steps and j time steps.

    dt=float(It[1]-It[0])/N
    dx=float(Ix[1]-Ix[0])/M
    W=[[0 for j in range(N+1)] for i in range(M+1)]
    Xs=np.linspace(Ix[0],Ix[1],M+1)
    Ts=np.linspace(It[0],It[1],N+1)

    for x in range(M+1):
        W[x][0]=b(Xs[x])

    for t in range(N+1):
        W[0][t]=l(Ts[t])

    for t in range(1,N+1):
        for x in range(1,M):
            W[x][t]=D*dt*float(W[x-1][t-1]-2*W[x][t-1]+W[x+1][t-1])/(dx**2)+W[x][t-1]

        W[M][t]=D*dt*float(r(Ts[t-1]))/dx-D*dt*float(W[M][t-1]-W[M-2][t-1])/(2*dx**2)+W[M][t-1]

    return W
