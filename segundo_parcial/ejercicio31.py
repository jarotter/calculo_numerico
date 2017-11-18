import odes
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def y_prime(t,y):
    return 2*(t+1)*y

# Exercise 3.1.1
t0=0
y0=1
step_size=0.1
step_num=10
w=odes.explicit_euler(y_prime,t0=t0,y0=y0,step_size=0.1,step_num=10)
print('Resultados de Euler explícito:')
print(' t  |     w_t')
for i in range(11):
    print('{0:.1f}'.format(t0+i*step_size),'|','{0:10f}'.format(w[i]))

print()

# Exercise 3.1.2
# Part 1
def exact_solution(t0,y0,t):
    return y0*np.exp(t*t-t0*t0+2*(t-t0))
print('Error global de la aproximación:',abs(exact_solution(t0,y0,1)-w[-1]),'\n')

# Part 2
max_error=0
for i in range(1,step_num+1):
    max_error=max(max_error,abs(exact_solution(t0+(i-1)*step_size,w[i-1],t0+i*step_size)-w[i]))

print('Error máximo local:',max_error,'\n')

# Part 3
h=[0.1*2**(-i) for i in range(6)]
errors=[]
for i in range(6):
    errors.append(abs(exact_solution(t0,y0,1)-odes.explicit_euler(y_prime,t0=t0,y0=y0,step_size=h[i],step_num=10*(2**i))[-1]))

plt.loglog(h, errors, marker = '.')
plt.xlabel('h')
plt.ylabel('global error at 1')
plt.title('Global error at ' r'$t=1$' + ' as a function of ' + r'$h = 0.1×2^{-k}$')
plt.grid(True)
plt.savefig('plot31.eps', format = 'eps')

# Part 4
local_errors=[]
for i in range(6):
    max_error=0
    w=odes.explicit_euler(y_prime,t0=t0,y0=y0,step_size=h[i],step_num=10*2**(i))
    for j in range(1,10*(2**i)+1):
        max_error=max(max_error,abs(exact_solution(t0+(j-1)*h[i],w[j-1],t0+j*h[i])-w[j]))

    local_errors.append(max_error)

print('k |   paso   | máximo de los errores locales del método de Euler | eoc')
for i in range(5):
    print(i,'|','{0:.6f}'.format(h[i]),'|                 ','{0:.13f}'.format(local_errors[i]),'                 |','{0:.10f}'.format(np.log(errors[i+1]/errors[i])/np.log(h[i+1]/h[i])))

print(5,'|','{0:.6f}'.format(h[5]),'|                 ','{0:.13f}'.format(local_errors[5]),'                 |','------------')
