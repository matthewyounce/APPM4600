import numpy as np
import time
from numpy.linalg import inv 
from numpy.linalg import norm

def driver():
    
    x0 = np.array([1,1])
    
    Nmax = 100
    tol = 1e-10
    
    t = time.time()
    for j in range(50):
      [xstar,ier,its,normlst] =  Newton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Newton: the error message reads:',ier) 
    print('Newton: took this many seconds:',elapsed/50)
    print('Netwon: number of iterations is:',its)
    compute_order(normlst[:-1],norm(xstar))
    
    
    for j in range(20):
      [xstar,ier,its,normlst] =  LazyNewton(x0,tol,Nmax)
    elapsed = time.time()-t
    print(xstar)
    print('Lazy Newton: the error message reads:',ier)
    print('Lazy Newton: took this many seconds:',elapsed/20)
    print('Lazy Newton: number of iterations is:',its)
    compute_order(normlst[:-1],norm(xstar))
    
def evalF(x): 

    F = np.zeros(2)
    
    F[0]=3*x[0]**2-x[1]**2
    F[1]=3*x[0]*x[1]**2-x[0]**3-1
    return F
    
def evalJ(x):
    
    J=np.array([[6*x[0],-2*x[1]],
                [3*x[1]**2-3*x[0]**2,6*x[0]*x[1]]])
    
    return J


def LazyNewton(x0,tol,Nmax):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    normlst=[norm(x0)]
    Jinv = np.array([[1/6,1/18],
                     [0,1/6]])
    for its in range(Nmax):

       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       normlst.append(norm(x1))
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its, normlst]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]

def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    normlst=[norm(x0)]
    for its in range(Nmax):
       J = evalJ(x0)
       Jinv = inv(J)
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       normlst.append(norm(x1))
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its,normlst]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]

def compute_order(x,xstar):
    diff1 = np.abs(x[1::]-xstar)
    diff2 = np.abs(x[0:-1]-xstar)
    
    fit=np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)
    _lambda=np.exp(fit[1])
    alpha=fit[0]
    print(f'lambda is {_lambda}')
    print(f'alpha is {alpha}')
    return fit

driver()