import numpy as np
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

def driver():

    #part 1
    #x0 = np.array([1, 1])
    
    #part 2
    #x0 = np.array([1, -1])
    
    #part 3
    x0 = np.array([0, 0])
    
    Nmax = 200
    tol = 1e-8
    
    #t = time.time()
    #for j in range(50):
      #[xstar,ier,its,seq] =  Newton(x0,tol,Nmax)
    #elapsed = time.time()-t
    #print(xstar)
    #print('Newton: the error message reads:',ier) 
    #print('Newton: took this many seconds:',elapsed/50)
    #print('Netwon: number of iterations is:',its)
    #compute_order(seq[:-1],norm(xstar))
    
     
    #t = time.time()
    #for j in range(20):
      #[xstar,ier,its,seq] =  LazyNewton(x0,tol,Nmax)
    #elapsed = time.time()-t
    #print(xstar)
    #print('Lazy Newton: the error message reads:',ier)
    #print('Lazy Newton: took this many seconds:',elapsed/20)
    #print('Lazy Newton: number of iterations is:',its)
    #compute_order(seq[:-1],norm(xstar))
    
     
    t = time.time()
    for j in range(20):
      [xstar,ier,its,seq] = Broyden(x0, tol,Nmax)     
    elapsed = time.time()-t
    print(xstar)
    print('Broyden: the error message reads:',ier)
    print('Broyden: took this many seconds:',elapsed/20)
    print('Broyden: number of iterations is:',its)
    compute_order(seq[:-1],norm(xstar))
     
def evalF(x): 

    F = np.zeros(2)
    
    F[0] = x[0]**2+x[1]**2-4
    F[1] = np.e**x[0]+x[1]-1
    return F
    
def evalJ(x): 

    
    J = np.array([[2*x[0],2*x[1]],
                  [np.e**x[0],1]])
    return J


def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    seq=[norm(x0)]
    for its in range(Nmax):
       J = evalJ(x0)
       Jinv = inv(J)
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       seq.append(norm(x1))
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its,seq]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its,seq]
           
def LazyNewton(x0,tol,Nmax):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    seq=[norm(x0)]
    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):

       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       seq.append(norm(x1))
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its,seq]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its,seq]   
    
def Broyden(x0,tol,Nmax):
    '''tol = desired accuracy
    Nmax = max number of iterations'''

    '''Sherman-Morrison 
   (A+xy^T)^{-1} = A^{-1}-1/p*(A^{-1}xy^TA^{-1})
    where p = 1+y^TA^{-1}Ax'''

    '''In Newton
    x_k+1 = xk -(G(x_k))^{-1}*F(x_k)'''


    '''In Broyden 
    x = [F(xk)-F(xk-1)-\hat{G}_k-1(xk-xk-1)
    y = x_k-x_k-1/||x_k-x_k-1||^2'''

    ''' implemented as in equation (10.16) on page 650 of text'''
    
    '''initialize with 1 newton step'''
    seq=[norm(x0)]
    A0 = evalJ(x0)

    v = evalF(x0)
    A = np.linalg.inv(A0)

    s = -A.dot(v)
    xk = x0+s
    seq.append(norm(xk))
    for  its in range(Nmax):
       '''(save v from previous step)'''
       w = v
       ''' create new v'''
       v = evalF(xk)
       '''y_k = F(xk)-F(xk-1)'''
       y = v-w;                   
       '''-A_{k-1}^{-1}y_k'''
       z = -A.dot(y)
       ''' p = s_k^tA_{k-1}^{-1}y_k'''
       p = -np.dot(s,z)                 
       u = np.dot(s,A) 
       ''' A = A_k^{-1} via Morrison formula'''
       tmp = s+z
       tmp2 = np.outer(tmp,u)
       A = A+1./p*tmp2
       ''' -A_k^{-1}F(x_k)'''
       s = -A.dot(v)
       xk = xk+s
       seq.append(norm(xk))
       if (norm(s)<tol):
          alpha = xk
          ier = 0
          return[alpha,ier,its,seq]
    alpha = xk
    ier = 1
    return[alpha,ier,its,seq]


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