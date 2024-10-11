import numpy as np
from numpy.linalg import norm

def driver():
    x0=np.array([1,1,1])
    
    Nmax = 100
    tol = 10e-8
    
    [xstar,ier,its,lst] =  it_meth(x0,tol,Nmax)
    print(xstar)
    print('it_meth: the error message reads:',ier)
    print('it_meth: number of iterations is:',its)
    compute_order(lst[:-1],norm(xstar))
    

def evald(x): 

    F = (x[0]**2+4*x[1]**2+4*x[2]**2-16)/(4*x[0]**2+64*x[1]**2+64*x[2]**2)
    return F

def it_meth(x0,tol,Nmax):

    normlst=[norm(x0)]
    for its in range(Nmax):

       d = evald(x0)
       x1 =np.zeros(3)
       
       x1[0]=x0[0]-d*2*x0[0]
       x1[1]=x0[1]-d*8*x0[1]
       x1[2]=x0[2]-d*8*x0[2]
       
       normlst.append(norm(x1))
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its,normlst]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its, normlst]

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
    