import numpy as np
from numpy.linalg import inv 
from numpy.linalg import norm 
import time

def driver():
    x0 = np.array([.1,.1,-.1])
    
    Nmax = 100
    tol = 10e-10
    
    t = time.time()
    for j in range(20):
      [xstar,ier,its] = Newton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Lazy Newton: the error message reads:',ier)
    print('Lazy Newton: took this many seconds:',elapsed/20)
    print('Lazy Newton: number of iterations is:',its)


def evalF(x): 

    F = np.zeros(3)
    
    F[0] = 4**x[0]**2+x[1]**2-4
    F[1] = x[0]+x[1]-np.sin(x[0]-x[1])
    
    return F
    
def aproxevalJ(x,h): 

    j1=(4*(x[0]+h)**2+x[1]**2-4-(4**x[0]**2+x[1]**2-4))/h
    j2=(4*x[0]**2+(x[1]+h)**2-4-(4**x[0]**2+x[1]**2-4))/h
    j3=(x[0]+h+x[1]-np.sin(x[0]+h-x[1])-(x[0]+x[1]-np.sin(x[0]-x[1])))/h
    j4=(x[0]+x[1]+h-np.sin(x[0]-x[1]-h)-(x[0]+x[1]-np.sin(x[0]-x[1])))/h
    J = np.array([[j1,j2], 
        [j3, j4]]) 
    
    return J

def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
       J = aproxevalJ(x0,10e-3)
       Jinv = inv(J)
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]

driver()