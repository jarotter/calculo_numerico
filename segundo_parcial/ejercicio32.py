import odes
import numpy as np

def y_prime(t,y):
    return 2*(t+1)*y

def exact_solution(t0,y0,t):
    return y0*np.exp(t*t-t0*t0+2*(t-t0))

# Ejercicio 3.2.1
wT1=odes.explicit_trapezoid(y_prime,t0=0,y0=1,step_size=0.1,step_num=10)
error_t1=abs(wT1[-1]-exact_solution(0,1,1))
print('Error global de trapecio explícito en t=1 con h=1/10:',error_t1,'\n')

# Ejercicio 3.2.2
wT2=odes.explicit_trapezoid(y_prime,t0=0,y0=1,step_size=0.05,step_num=20)
errorT2=abs(wT2[-1]-exact_solution(0,1,1))
print('Error global de trapecio explícito en t=1 con h=1/20:',errorT2)

wrk4=odes.rk4(y_prime,t0=0,y0=1,step_size=0.1,step_num=10)
error_rk=abs(wrk4[-1]-exact_solution(0,1,1))
print('Error global de rk4 en t=1 con h=1/10:',error_rk,'\n')

# Ejercicio 3.2.3
print('                    | Trapecio(1/10) | Trapecio(1/20) |    rk4(1/10)')
print('Error global en t=1 | ','{0:.10f}'.format(error_t1),' | ','{0:.10f}'.format(errorT2),' | ','{0:.10f}'.format(error_rk))
print('Número de estados   |       20       |       40       |       40')

# Ejercicio 3.2.4
# En este caso, entre reducir los intervalos y cambiar de método a rk4,
# la forma más conveniente de reducir el error es cambiar a un método
# de orden mayor, ya que para el método del trapecio, el error global es
# O(h^2), así que al dividir el intervalo entre 2, reducimos el error en
# un factor de 4. En cambio, al cambiar a rk4, el error es O(h^4), así que
# manteniendo la misma h=1/10, reducimos el error global en un factor de
# (1/100)/(1/10000)=100
