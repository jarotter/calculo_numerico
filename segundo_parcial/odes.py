import numpy as np

def derivate(func,x0,h=1e-6):
    # Approximates the derivative of func at x0

    return (func(x0+h)-func(x0-h))/(2*h)

def newtonSolve(func,val,tol=1e-15):
    # Approximates a root of function func with initial approximation val
    # using Newton's method

    prevval=val
    for i in range(100):
        while abs(derivate(func,val))<1e-10: val=val+0.1
        val=val-func(val)/derivate(func,val)


        if abs(prevval-val)<tol: return val
        prevval=val

###############################################################################
#*****************************************************************************#
###############################################################################


def explicitEuler_1D(func,t0,y0,stepSize=0.1,stepNum=50):
    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        yValues[i]=yValues[i-1]+stepSize*func(tValues[i-1],yValues[i-1])

    return yValues

def explicitTrapeze_1D(func,t0,y0,stepSize=0.1,stepNum=50):
    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        nextYAprox=yValues[i-1]+stepSize*func(tValues[i-1],yValues[i-1])
        yValues[i]=yValues[i-1]+stepSize*(func(tValues[i-1],yValues[i-1])+func(tValues[i],nextYAprox))/2

    return yValues

def midpointMethod_1D(func,t0,y0,stepSize=0.1,stepNum=50):
    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        m=yValues[i-1]+stepSize*func(tValues[i-1],yValues[i-1])/2
        yValues[i]=yValues[i-1]+stepSize*func(tValues[i-1]+stepSize/2,m)

    return yValues

def RK4_1D(func,t0,y0,stepSize=0.1,stepNum=50):
    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        k1=func(tValues[i-1],yValues[i-1])
        k2=func(tValues[i-1]+stepSize/2,yValues[i-1]+stepSize*k1/2)
        k3=func(tValues[i-1]+stepSize/2,yValues[i-1]+stepSize*k2/2)
        k4=func(tValues[i],yValues[i-1]+stepSize*k3)

        yValues[i]=yValues[i-1]+stepSize*(k1+2*k2+2*k3+k4)/6

    return yValues

def RK23_1D(func,t0,y0,t,tol=1e-1):
    tValues=[t0]
    yValues=[y0]

    while tValues[-1]+1e-6<t:
        h=0.1

        for i in range(5):
            k1=func(tValues[-1],yValues[-1])
            k2=func(tValues[-1]+h/2,yValues[-1]+h*k1/2)
            k3=func(tValues[-1]+3*h/4,yValues[-1]+3*h*k2/4)
            k4=func(tValues[-1]+h,yValues[-1]+h*(2*k1+3*k2+4*k3)/9)

            ord2=yValues[-1]+h*(7*k1+6*k2+8*k3+3*k4)/24
            ord3=yValues[-1]+h*(2*k1+3*k2+4*k3)/9

            if abs(ord2-ord3)==0:
                h=1
                break
            elif abs(ord2-ord3)>tol*0.95 and abs(ord2-ord3)<1.05*tol:
                break
            else:
                h=min(1,t-tValues[-1],0.1*(tol/abs(ord2-ord3))**(1/2))
                h=max(0.01,h)
                break

            if h<=0.011 or h>=1.001: break

        k1=func(tValues[-1],yValues[-1])
        k2=func(tValues[-1]+h/2,yValues[-1]+h*k1/2)
        k3=func(tValues[-1]+3*h/4,yValues[-1]+3*h*k2/4)
        k4=func(tValues[-1]+h,yValues[-1]+h*(2*k1+3*k2+4*k3)/9)

        tValues.append(h+tValues[-1])
        yValues.append(yValues[-1]+h*(7*k1+6*k2+8*k3+3*k4)/24)

    return tValues,yValues


###############################################################################
#*****************************************************************************#
###############################################################################

def explicitEuler_MD(func,t0,y0,stepSize=0.1,stepNum=50):
    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        yValues[i]=yValues[i-1]+stepSize*func(tValues[i-1],*yValues[i-1])

    return yValues

def explicitTrapeze_MD(func,t0,y0,stepSize=0.1,stepNum=50):
    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        nextYAprox=yValues[i-1]+stepSize*func(tValues[i-1],*yValues[i-1])
        yValues[i]=yValues[i-1]+stepSize*(func(tValues[i-1],*yValues[i-1])+func(tValues[i],*nextYAprox))/2

    return yValues

def midpointMethod_MD(func,t0,y0,stepSize=0.1,stepNum=50):
    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        m=yValues[i-1]+stepSize*func(tValues[i-1],*yValues[i-1])/2
        yValues[i]=yValues[i-1]+stepSize*func(tValues[i-1]+stepSize/2,*m)

    return yValues

def RK4_MD(func,t0,y0,stepSize=0.1,stepNum=50):
    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        k1=func(tValues[i-1],*yValues[i-1])
        k2=func(tValues[i-1]+stepSize/2,*(yValues[i-1]+stepSize*k1/2))
        k3=func(tValues[i-1]+stepSize/2,*(yValues[i-1]+stepSize*k2/2))
        k4=func(tValues[i],*(yValues[i-1]+stepSize*k3))

        yValues[i]=yValues[i-1]+stepSize*(k1+2*k2+2*k3+k4)/6

    return yValues

###############################################################################
#*****************************************************************************#
###############################################################################

def explicitEuler(func,t0,y0,stepSize=0.1,stepNum=50):
    # Approximates the value of y at time t0+stepSize*k, for k=0,...,stepNum,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using Euler's explicit method

    if (type(y0)==int) or (type(y0)==float): return explicitEuler_1D(func,t0,y0,stepSize,stepNum)
    else: return explicitEuler_MD(func,t0,y0,stepSize,stepNum)

def explicitTrapeze(func,t0,y0,stepSize=0.1,stepNum=50):
    # Approximates the value of y at time t0+stepSize*k, for k=0,...,stepNum,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using trapeze explicit method

    if (type(y0)==int) or (type(y0)==float): return explicitTrapeze_1D(func,t0,y0,stepSize,stepNum)
    else: return explicitTrapeze_MD(func,t0,y0,stepSize,stepNum)

def implicitEuler(func,t0,y0,stepSize=0.1,stepNum=50):
    # Approximates the value of y at time t0+stepSize*k, for k=0,...,stepNum,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using Euler's explicit method

    tValues=[t0+stepSize*i for i in range(stepNum+1)]
    yValues=[0 for i in range(stepNum+1)]

    yValues[0]=y0
    for i in range(1,stepNum+1):
        def findSig(y):
            return y-yValues[i-1]-stepSize*func(tValues[i],y)

        yValues[i]=newtonSolve(findSig,yValues[i-1]+stepSize*func(tValues[i-1],yValues[i-1]))

    return yValues

def midpointMethod(func,t0,y0,stepSize=0.1,stepNum=50):
    # Approximates the value of y at time t0+stepSize*k, for k=0,...,stepNum,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using midpoint method

    if (type(y0)==int) or (type(y0)==float): return midpointMethod_1D(func,t0,y0,stepSize,stepNum)
    else: return midpointMethod_MD(func,t0,y0,stepSize,stepNum)

def RK23(func,t0,y0,t,tol):
    # Approximates a trajectory y considering that between any two step the approximation
    # is of the same order as tol and that the derivative of y equals func(t,y) and y0=y(t0)
    # using Bogacki-Shampine method

    return RK23_1D(func,t0,y0,t,tol)

def RK4(func,t0,y0,stepSize=0.1,stepNum=50):
    # Approximates the value of y at time t0+stepSize*k, for k=0,...,stepNum,
    # considering that the derivative of y equals func(t,y) and y0=y(t0)
    # using Runge Kutta-4

    if (type(y0)==int) or (type(y0)==float): return RK4_1D(func,t0,y0,stepSize,stepNum)
    else: return RK4_MD(func,t0,y0,stepSize,stepNum)
