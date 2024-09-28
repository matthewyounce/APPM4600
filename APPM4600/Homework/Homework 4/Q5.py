# import libraries
import numpy as np
import matplotlib.pyplot as plt
        
def driver():
  f = lambda x: x**6-x-1
  fp = lambda x: 6*x**5-1
  p0 = 2

  Nmax = 100
  tol = 1e-8

  (pn,pstarn,info,it) = newton(f,fp,p0,tol, Nmax)
  print('the approximate root is', pstarn)
  print('the error message reads:', '%d' % info)
  print('Number of iterations:', '%d' % it)
  #fetting the values for my table
  print('Newton row of table =',pn-pstarn)
  n=pn-pstarn
  
  #fit=compute_order(pn,pstarn)
  
  (ps,pstars,info) = secant(2,1,f,tol, Nmax)
  print('the approximate root is', '%16.16e' % pstars)
  print('the error message reads:', '%d' % info)
  #getting the values for my table
  print('Secant row of table',ps-pstars)
  s=ps-pstars
  xn=abs(n[:-1])
  yn=abs(n[1:])
  xs=abs(s[:-1])
  ys=abs(s[1:])
  #fit=compute_order(ps[:-1],pstars)
  
  xline=np.linspace(.00001,1,1000)
  ylinen=xline**2
  ylines=xline**1.5
  
  #ploting newton method
  plt.plot(xline,ylinen,'red')
  plt.plot(xn,yn)
  plt.xscale('log')
  plt.yscale('log')
  plt.title('Newton')
  plt.ylabel('|x{k+1}-alpha|')
  plt.xlabel('|x{k}-alpha|')
  plt.show()
  
  #plotting secend method
  plt.plot(xline,ylines,'red')
  plt.plot(xs,ys)
  plt.xscale('log')
  plt.yscale('log')
  plt.title('Secant')
  plt.ylabel('|x{k+1}-alpha|')
  plt.xlabel('|x{k}-alpha|')
  plt.show()
  
  
  


def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1);
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p[:it+1],pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]
        

def secant(x0,x1,f,tol,Nmax):
    p=np.zeros(Nmax+2)
    p[0]=x0
    p[1]=x1
    if abs(f(x1)-f(x0))==0:
        return [p,x1,1]
    for j in range(Nmax):
        x2=x1-f(x1)*(x1-x0)/(f(x1)-f(x0))
        p[j+2]=x2
        if abs(x2-x1)<tol:
            return [p[:j+3],x2,0]
        x0=x1
        x1=x2
        if abs(f(x1)-f(x0))==0:
            return [p,x1,1]
    return [p,x2,1]

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
