import odes
import numpy as np
import math

def exact_solution(t) :
    y = 2 * (1-1*math.exp(4*t)) /(1+1*math.exp(4*t))
    return y

def y_prime(t,y):
    return y*y - 4

t0=0
y0=0
h = [1/50,1/100,1/200]
n = [100,200,400]

errors_euler = np.zeros(3)
errors_Trapezoid = np.zeros(3)
errors_midpoint = np.zeros(3)
errors_rk4 = np.zeros(3)

for i in range(3):

    euler_approximation = odes.explicitEuler(y_prime,t0=t0,y0=y0,stepSize=h[i],stepNum=n[i])
    trapezoid_approximation = odes.explicitTrapezoid(y_prime,t0=t0,y0=y0,stepSize=h[i],stepNum=n[i])
    midpoint_approximation = odes.midpointMethod(y_prime,t0=t0,y0=y0,stepSize=h[i],stepNum=n[i])
    rk4_approximation = odes.RK4(y_prime,t0=t0,y0=y0,stepSize=h[i],stepNum=n[i])

    y = exact_solution(2)
    errors_euler[i] = (abs(euler_approximation[n[i]]-y))
    errors_Trapezoid[i] = (abs(trapezoid_approximation[n[i]]-y))
    errors_midpoint[i] = (abs(midpoint_approximation[n[i]]-y))
    errors_rk4[i] = (abs(rk4_approximation[n[i]]-y))


print('Global error using the explicit euler method:\n')
for i in range(2):
    print('   h =', h[i],': ', errors_euler[i], ', eoc: ',
          np.log(errors_euler[i+1]/errors_euler[i])/np.log(h[i+1]/h[i]))
print('   h =', h[2],': ', errors_euler[2], ', eoc: ','---------------')


print('\nGlobal error using the explicit trapezoid method:\n')
for i in range(2):
    print('   h =', h[i],': ', errors_Trapezoid[i], ', eoc: ',
          np.log(errors_Trapezoid[i+1]/errors_Trapezoid[i])/np.log(h[i+1]/h[i]))
print('   h =', h[2],': ', errors_Trapezoid[2], ', eoc: ','---------------')


print('\nGlobal error using the explicit midpoint method:\n')
for i in range(2):
    print('   h =', h[i],': ', errors_midpoint[i], ', eoc: ',
          np.log(errors_midpoint[i+1]/errors_midpoint[i])/np.log(h[i+1]/h[i]) )
print('   h =', h[2],': ', errors_midpoint[2], ', eoc: ','---------------')


print('\nGlobal error using the Runge Kutta 4 method:\n')
for i in range(2):
    print('   h =', h[i],': ', errors_rk4[i], ', eoc: ',
          np.log(errors_rk4[i+1]/errors_rk4[i])/np.log(h[i+1]/h[i]))
print('   h =', h[2],': ', errors_rk4[2], ', eoc: ','---------------')
