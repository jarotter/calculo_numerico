import odes
import numpy as np

def prime_y(t,y):
    return 2*(t+1)*y

def correct_y(t0,y0,t):
    return y0*np.exp(t*t-t0*t0+2*(t-t0))

# Ejercicio 3.2.1
wT1=odes.explicitTrapezoid(prime_y,t0=0,y0=1,stepSize=0.1,stepNum=10)
errorT1=abs(wT1[-1]-correct_y(0,1,1))
print('Error global de trapecio explícito en t=1 con h=1/10:',errorT1,'\n')

# Ejercicio 3.2.2
wT2=odes.explicitTrapezoid(prime_y,t0=0,y0=1,stepSize=0.05,stepNum=20)
errorT2=abs(wT2[-1]-correct_y(0,1,1))
print('Error global de trapecio explícito en t=1 con h=1/20:',errorT2)

wRK4=odes.RK4(prime_y,t0=0,y0=1,stepSize=0.1,stepNum=10)
errorRK=abs(wRK4[-1]-correct_y(0,1,1))
print('Error global de RK4 en t=1 con h=1/10:',errorRK,'\n')

# Ejercicio 3.2.3
print('                    | Trapecio(1/10) | Trapecio(1/20) |    RK4(1/10)')
print('Error global en t=1 | ','{0:.10f}'.format(errorT1),' | ','{0:.10f}'.format(errorT2),' | ','{0:.10f}'.format(errorRK))
print('Número de estados   |       20       |       40       |       40')

# Ejercicio 3.2.4
# En este caso, entre reducir los intervalos y cambiar de método a RK4,
# la forma más conveniente de reducir el error es cambiar a un método
# de orden mayor, ya que para el método del trapecio, el error global es
# O(h^2), así que al dividir el intervalo entre 2, reducimos el error en
# un factor de 4. En cambio, al cambiar a RK4, el error es O(h^4), así que
# manteniendo la misma h=1/10, reducimos el error global en un factor de
# (1/100)/(1/10000)=100
