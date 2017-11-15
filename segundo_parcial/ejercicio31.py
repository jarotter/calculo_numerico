import odes
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def prime_y(t,y):
    return 2*(t+1)*y

# Exercise 3.1.1
t0=0
y0=1
stepSize=0.1
stepNum=10
w=odes.explicitEuler(prime_y,t0=t0,y0=y0,stepSize=0.1,stepNum=10)
print('Resultados de Euler explícito:')
print(' t  |     w_t')
for i in range(11):
    print('{0:.1f}'.format(t0+i*stepSize),'|','{0:10f}'.format(w[i]))

print()

# Exercise 3.1.2
# Part 1
def correct_y(t0,y0,t):
    return y0*np.exp(t*t-t0*t0+2*(t-t0))
print('Error global de la aproximación:',abs(correct_y(t0,y0,1)-w[-1]),'\n')

# Part 2
maxError=0
for i in range(1,stepNum+1):
    maxError=max(maxError,abs(correct_y(t0+(i-1)*stepSize,w[i-1],t0+i*stepSize)-w[i]))

print('Error máximo local:',maxError,'\n')

# Part 3
h=[0.1*2**(-i) for i in range(6)]
errors=[]
for i in range(6):
    errors.append(abs(correct_y(t0,y0,1)-odes.explicitEuler(prime_y,t0=t0,y0=y0,stepSize=h[i],stepNum=10*(2**i))[-1]))

#for i in range(6):
    #errors[i]=np.log(errors[i])

# print('Log de los errores globales como función de k:')
# print('k | log de error')
# for i in range(6):
#     print(i,'|','{0:10f}'.format(errors[i]))
#
# print()

plt.loglog(h, errors, marker = '.')
plt.xlabel('h')
plt.ylabel('global error at 1')
plt.title('Global error at ' r'$t=1$' + ' as function of ' + r'$h = 0.1×2^{-k}$')
plt.savefig('plot31.eps', format = 'eps')
plt.show()

# Part 4
errors=[]
for i in range(6):
    maxError=0
    w=odes.explicitEuler(prime_y,t0=t0,y0=y0,stepSize=h[i],stepNum=10*2**(i))
    for j in range(1,10*(2**i)):
        maxError=max(maxError,abs(correct_y(t0+(j-1)*h[i],w[j-1],t0+j*h[i])-w[j]))

    errors.append(maxError)

print('k |   paso   | máximo de los errores locales del método de Euler | eoc')
for i in range(5):
    print(i,'|','{0:.6f}'.format(h[i]),'|                 ','{0:.13f}'.format(errors[i]),'                 |','{0:.10f}'.format(np.log(errors[i+1]/errors[i])/np.log(h[i+1]/h[i])))

print(5,'|','{0:.6f}'.format(h[5]),'|                 ','{0:.13f}'.format(errors[5]),'                 |','------------')
