import numpy as np
import math

# 6. Nonlinear Finite Difference Method for BVP
# Uses Multivariate Newton’s Method to solve a nonlinear equation
# Inputs: interval (interval), boundary values (boundary_values),
#         number of steps (n)

def exact_solution(t):
    n = len(t)
    y = np.zeros(n)

    for i in range(n):
        y[i] = 1 - math.pi/20 * math.tan(math.pi*t[i]/4)
    return y

def f(w,h,ya,yb,n):
    y=np.zeros(n); # initialize solution array

    y[0]= (ya - 2*w[0] + w[1]) - h*5*(w[1]-ya)*(1-w[0]);
    y[n-1]= (w[n-2] - 2*w[n-1] + yb) -h*5*(yb-w[n-2])*(1-w[n-1]);

    for j in range(1, n-1):
        y[j] = (w[j-1] - 2*w[j] + w[j+1]) - h*5*(w[j+1]-w[j-1])*(1-w[j]);

    return y

def jacobian(w,h,ya,yb,n):
    M= np.zeros((n,n));

    M[0][0]= -2 +5*h*(w[1]-ya);
    M[n-1][n-1]= -2 +5*h*(yb-w[n-3]);

    for i in range(1,n-1):
        M[i][i]= -2 +5*h*(w[i+1]-w[i-1]);

    for i in range(0,n-1):
        M[i][i+1]=1-5*h*(1-w[i]);
        M[i+1][i]=1+5*h*(1-w[i]);

    return M


def non_linear_finite_difference_method(interval,boundary_values,n):
    a=interval[0]; b=interval[1];
    ya=boundary_values[0]; yb=boundary_values[1];

    h=(b-a)/(n+1); # Step size (uniform)
    w=np.zeros(n);

    # loop of Newton step
    for i in range(20):
        e = np.linalg.solve(jacobian(w,h,ya,yb,n),f(w,h,ya,yb,n))
        w = w-e;
        if (np.linalg.norm(e,np.inf)< 1e-9):
            break

    #concatenating the boundary values to the aproximation vector
    ans = np.zeros(n+2);
    ans[0] = ya;
    for i in range(1,n+1):
        ans[i] = w[i-1]
    ans[n+1] = yb;

    return ans

interval = [0,1];
boundary_values = [1,1-math.pi/20];

print('Errores globales:')

# h = 1/20
n = 19
h=1/n
# Aproximation with finite differences
w = non_linear_finite_difference_method(interval,boundary_values,n)
# Exact solution
t = [ i*h for i in range(n+2)]
y = exact_solution(t)
# Global error
print('h=1/20: ', max(abs(w-y)))
print(abs(w[-1] - y[-1]))

# h=1/40
n = 39
h=1/n
w = non_linear_finite_difference_method(interval,boundary_values,n)
t = [ i*h for i in range(n+2)]
y = exact_solution(t)
print('h=1/40: ',max(abs(w-y)))
print(abs(w[-1] - y[-1]))

# h=1/80
n = 79
h=1/n
w = non_linear_finite_difference_method(interval,boundary_values,n)
t = [ i*h for i in range(n+2)]
y = exact_solution(t)
print('h=1/80: ',max(abs(w-y)))
print(abs(w[-1] - y[-1]))
