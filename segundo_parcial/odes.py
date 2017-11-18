import numpy as np

def derivative(func,x0,h=1e-6):
    # Approximates the derivative of func at x0

    return (func(x0+h)-func(x0-h))/(2*h)

def newton_solve(func,val,tol=1e-15):
    # Approximates a root of function func with initial approximation val
    # using Newton's method

    prevval=val
    for i in range(100):
        while abs(derivative(func,val))<1e-10: val=val+0.1
        val=val-func(val)/derivative(func,val)


        if abs(prevval-val)<tol: return val
        prevval=val

###############################################################################
#*****************************************************************************#
###############################################################################


def explicit_euler_1d(func,t0,y0,step_size=0.1,step_num=50):
    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        y_values[i]=y_values[i-1]+step_size*func(t_values[i-1],y_values[i-1])

    return y_values

def explicit_trapezoid_1D(func,t0,y0,step_size=0.1,step_num=50):
    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        nextYAprox=y_values[i-1]+step_size*func(t_values[i-1],y_values[i-1])
        y_values[i]=y_values[i-1]+step_size*(func(t_values[i-1],y_values[i-1])+func(t_values[i],nextYAprox))/2

    return y_values

def midpoint_method_1D(func,t0,y0,step_size=0.1,step_num=50):
    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        m=y_values[i-1]+step_size*func(t_values[i-1],y_values[i-1])/2
        y_values[i]=y_values[i-1]+step_size*func(t_values[i-1]+step_size/2,m)

    return y_values

def rk4_1D(func,t0,y0,step_size=0.1,step_num=50):
    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        k1=func(t_values[i-1],y_values[i-1])
        k2=func(t_values[i-1]+step_size/2,y_values[i-1]+step_size*k1/2)
        k3=func(t_values[i-1]+step_size/2,y_values[i-1]+step_size*k2/2)
        k4=func(t_values[i],y_values[i-1]+step_size*k3)

        y_values[i]=y_values[i-1]+step_size*(k1+2*k2+2*k3+k4)/6

    return y_values

def rk23_1D(func,t0,y0,t,tol=1e-1):
    t_values=[t0]
    y_values=[y0]

    while t_values[-1]+1e-6<t:
        h=0.1

        for i in range(5):
            k1=func(t_values[-1],y_values[-1])
            k2=func(t_values[-1]+h/2,y_values[-1]+h*k1/2)
            k3=func(t_values[-1]+3*h/4,y_values[-1]+3*h*k2/4)
            k4=func(t_values[-1]+h,y_values[-1]+h*(2*k1+3*k2+4*k3)/9)

            ord2=y_values[-1]+h*(7*k1+6*k2+8*k3+3*k4)/24
            ord3=y_values[-1]+h*(2*k1+3*k2+4*k3)/9

            if abs(ord2-ord3)==0:
                h=1
                break
            elif abs(ord2-ord3)>tol*0.95 and abs(ord2-ord3)<1.05*tol:
                break
            else:
                h=min(1,t-t_values[-1],0.1*(tol/abs(ord2-ord3))**(1/2))
                h=max(0.01,h)
                break

            if h<=0.011 or h>=1.001: break

        k1=func(t_values[-1],y_values[-1])
        k2=func(t_values[-1]+h/2,y_values[-1]+h*k1/2)
        k3=func(t_values[-1]+3*h/4,y_values[-1]+3*h*k2/4)
        k4=func(t_values[-1]+h,y_values[-1]+h*(2*k1+3*k2+4*k3)/9)

        t_values.append(h+t_values[-1])
        y_values.append(y_values[-1]+h*(7*k1+6*k2+8*k3+3*k4)/24)

    return t_values,y_values


###############################################################################
#*****************************************************************************#
###############################################################################

def explicit_euler_md(func,t0,y0,step_size=0.1,step_num=50):
    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        y_values[i]=y_values[i-1]+step_size*func(t_values[i-1],*y_values[i-1])

    return y_values

def explicit_trapezoid_md(func,t0,y0,step_size=0.1,step_num=50):
    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        nextYAprox=y_values[i-1]+step_size*func(t_values[i-1],*y_values[i-1])
        y_values[i]=y_values[i-1]+step_size*(func(t_values[i-1],*y_values[i-1])+func(t_values[i],*nextYAprox))/2

    return y_values

def midpoint_method_md(func,t0,y0,step_size=0.1,step_num=50):
    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        m=y_values[i-1]+step_size*func(t_values[i-1],*y_values[i-1])/2
        y_values[i]=y_values[i-1]+step_size*func(t_values[i-1]+step_size/2,*m)

    return y_values

def rk4_md(func,t0,y0,step_size=0.1,step_num=50):
    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        k1=func(t_values[i-1],*y_values[i-1])
        k2=func(t_values[i-1]+step_size/2,*(y_values[i-1]+step_size*k1/2))
        k3=func(t_values[i-1]+step_size/2,*(y_values[i-1]+step_size*k2/2))
        k4=func(t_values[i],*(y_values[i-1]+step_size*k3))

        y_values[i]=y_values[i-1]+step_size*(k1+2*k2+2*k3+k4)/6

    return y_values

###############################################################################
#*****************************************************************************#
###############################################################################

def explicit_euler(func,t0,y0,step_size=0.1,step_num=50):
    # Approximates the value of y at time t0+step_size*k, for k=0,...,step_num,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using Euler's explicit method

    if (type(y0)==int) or (type(y0)==float): return explicit_euler_1d(func,t0,y0,step_size,step_num)
    else: return explicit_euler_md(func,t0,y0,step_size,step_num)

def explicit_trapezoid(func,t0,y0,step_size=0.1,step_num=50):
    # Approximates the value of y at time t0+step_size*k, for k=0,...,step_num,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using trapezoid explicit method

    if (type(y0)==int) or (type(y0)==float): return explicit_trapezoid_1D(func,t0,y0,step_size,step_num)
    else: return explicit_trapezoid_md(func,t0,y0,step_size,step_num)

def implicit_euler(func,t0,y0,step_size=0.1,step_num=50):
    # Approximates the value of y at time t0+step_size*k, for k=0,...,step_num,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using Euler's explicit method

    t_values=[t0+step_size*i for i in range(step_num+1)]
    y_values=[0 for i in range(step_num+1)]

    y_values[0]=y0
    for i in range(1,step_num+1):
        def find_sig(y):
            return y-y_values[i-1]-step_size*func(t_values[i],y)

        y_values[i]=newton_solve(find_sig,y_values[i-1]+step_size*func(t_values[i-1],y_values[i-1]))

    return y_values

def midpoint_method(func,t0,y0,step_size=0.1,step_num=50):
    # Approximates the value of y at time t0+step_size*k, for k=0,...,step_num,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using midpoint method

    if (type(y0)==int) or (type(y0)==float): return midpoint_method_1D(func,t0,y0,step_size,step_num)
    else: return midpoint_method_md(func,t0,y0,step_size,step_num)

def rk23(func,t0,y0,t,tol):
    # Approximates a trajectory y considering that between any two step the approximation
    # is of the same order as tol and that the derivative of y equals func(t,y) and y0=y(t0)
    # using Bogacki-Shampine method

    return rk23_1D(func,t0,y0,t,tol)

def rk4(func,t0,y0,step_size=0.1,step_num=50):
    # Approximates the value of y at time t0+step_size*k, for k=0,...,step_num,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using Runge Kutta-4

    if (type(y0)==int) or (type(y0)==float): return rk4_1D(func,t0,y0,step_size,step_num)
    else: return rk4_md(func,t0,y0,step_size,step_num)
