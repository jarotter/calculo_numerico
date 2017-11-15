import odes
import numpy as np
from numpy import linalg

def correct_y(t):
    return np.array([np.exp(-t)*np.sin(t),np.exp(-t)*np.cos(t)])

def prime_y(t,y1,y2):
    return np.array([-y1+y2,-y1-y2])

# Ejercicio 3.3.1
w10=odes.midpointMethod(prime_y,t0=0,y0=np.array([0,1]),stepSize=0.1,stepNum=10)
w100=odes.midpointMethod(prime_y,t0=0,y0=np.array([0,1]),stepSize=0.01,stepNum=100)

print('Error global para método del punto medio con h=1/10:',linalg.norm(w10[-1]-correct_y(1)))
print('Error global para método del punto medio con h=1/100:',linalg.norm(w100[-1]-correct_y(1)))
print('(Errores bajo la norma 2)')

# Ejercicio 3.3.2
# Sí es consistente la reducción del error, porque la función que define la derivada es analítica
