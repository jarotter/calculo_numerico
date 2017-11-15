import odes
import numpy as np

def y_prime(t,y):
    return 5*y-3*y*y

t0=0
y0=1/2

# Ejercicio 3.4.1
wImp=odes.implicitEuler(y_prime,t0=t0,y0=y0,stepSize=1/2,stepNum=40)
print('y(20) aproximado con Euler implícito y h=1/2:',wImp[-1],'\n')

# Por lo tanto, la solución de este problema de valores iniciales tiende a
# 5/3

#Ejercicio 3.4.2
wExp6=odes.explicitEuler(y_prime,t0=t0,y0=y0,stepSize=1/6,stepNum=120)
wExp5=odes.explicitEuler(y_prime,t0=t0,y0=y0,stepSize=1/5,stepNum=100)
wExp4=odes.explicitEuler(y_prime,t0=t0,y0=y0,stepSize=1/4,stepNum=80)
wExp2=odes.explicitEuler(y_prime,t0=t0,y0=y0,stepSize=1/2,stepNum=40)

print('y(20) aproximado con Euler explícito y h=1/6:',wExp6[-1])
print('y(20) aproximado con Euler explícito y h=1/5:',wExp5[-1])
print('y(20) aproximado con Euler explícito y h=1/4:',wExp4[-1])
print('y(20) aproximado con Euler explícito y h=1/2:',wExp2[-1])

print()
#Ejercicio 3.4.3
wRK23_10=odes.RK23(y_prime,t0=t0,y0=y0,t=20,tol=0.1)
wRK23_100=odes.RK23(y_prime,t0=t0,y0=y0,t=20,tol=0.01)

print('y(20) aproximado con RK23 y tol=1/10:',wRK23_10[1][-1])
print('y(20) aproximado con RK23 y tol=1/100:',wRK23_100[1][-1])

# Es mejor el método implícito que el método explícito adaptativo ya que
# la ecuación presentada es muy "stiff", lo que hace que los métodos adaptativos
# tengan que trabajar con números muy pequeños, lo que acrecenta los errores en
# las FLOPS. Además, la ecuación resultante en el método implícito es cuadrática,
# lo que la hace muy fácil de aproximar con el método de Newton.
