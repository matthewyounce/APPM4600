import numpy as np

def driver():

# test functions 
     f1 = lambda x: (10/(x+4))**(1/2)
# fixed point is alpha1 = 1.36523

     Nmax = 100
     tol = 1e-10

# test f1 '''
     x0 = 1.5
     [xstar,ier,n] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f1(xstar):',f1(xstar))
     print('Error message reads:',ier)
     print('squence is', n)
     
     
     y=Aikens(n)
     print(y)
     
     xstar=n[-1]
     fit=compute_order(y,xstar)



# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    n=np.zeros((Nmax+1,1))
    n[0,0]=x0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       n[count,0]=x1
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          n=n[:count+1,0]
          return [xstar,ier, n]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier,n]
    
def Aikens(x):
    x0=x[:-2]
    x1=x[1:-1]
    x2=x[2:]
    y=x0-((x1-x0)**2)/(x2-2*x1+x0)
    return(y)
    
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



    