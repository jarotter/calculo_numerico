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

errorsEuler = np.zeros(3)
errorsTrapezoid = np.zeros(3)
errorsMidpoint = np.zeros(3)
errorsRK4 = np.zeros(3)

for i in range(3):

    eulerAproximation = odes.explicitEuler(y_prime,t0=t0,y0=y0,stepSize=h[i],stepNum=n[i])
    trapezoidAproximation = odes.explicitTrapezoid(y_prime,t0=t0,y0=y0,stepSize=h[i],stepNum=n[i])
    midpointAproximation = odes.midpointMethod(y_prime,t0=t0,y0=y0,stepSize=h[i],stepNum=n[i])
    RK4Aproximation = odes.RK4(y_prime,t0=t0,y0=y0,stepSize=h[i],stepNum=n[i])

    y = exact_solution(2)
    errorsEuler[i] = (abs(eulerAproximation[n[i]]-y))
    errorsTrapezoid[i] = (abs(trapezoidAproximation[n[i]]-y))
    errorsMidpoint[i] = (abs(midpointAproximation[n[i]]-y))
    errorsRK4[i] = (abs(RK4Aproximation[n[i]]-y))


print('Global error using the explicit euler method:\n')
for i in range(2):
    print('   h =', h[i],': ', errorsEuler[i], ', eoc: ',
          np.log(errorsEuler[i+1]/errorsEuler[i])/np.log(h[i+1]/h[i]))
print('   h =', h[2],': ', errorsEuler[2], ', eoc: ','---------------')


print('\nGlobal error using the explicit trapezoid method:\n')
for i in range(2):
    print('   h =', h[i],': ', errorsTrapezoid[i], ', eoc: ',
          np.log(errorsTrapezoid[i+1]/errorsTrapezoid[i])/np.log(h[i+1]/h[i]))
print('   h =', h[2],': ', errorsTrapezoid[2], ', eoc: ','---------------')


print('\nGlobal error using the explicit midpoint method:\n')
for i in range(2):
    print('   h =', h[i],': ', errorsMidpoint[i], ', eoc: ',
          np.log(errorsMidpoint[i+1]/errorsMidpoint[i])/np.log(h[i+1]/h[i]) )
print('   h =', h[2],': ', errorsMidpoint[2], ', eoc: ','---------------')


print('\nGlobal error using the Runge Kutta 4 method:\n')
for i in range(2):
    print('   h =', h[i],': ', errorsRK4[i], ', eoc: ',
          np.log(errorsRK4[i+1]/errorsRK4[i])/np.log(h[i+1]/h[i]))
print('   h =', h[2],': ', errorsRK4[2], ', eoc: ','---------------')
