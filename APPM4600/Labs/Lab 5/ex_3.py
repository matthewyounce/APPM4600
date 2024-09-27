# import libraries
import numpy as np

def driver():

# use routines    
    f = lambda x: np.e**(x**2+7*x-30)-1
    df=lambda x: np.e**(x**2+7*x-30)*(2*x+7)
    ddf=lambda x:np.e**(x**2+7*x-30)*2+np.e**(x**2+7*x-30)*(2*x+7)**2
    a = 2
    b = 4.5

    tol = 1e-8
    Nmax=100

    [astar,ier,count] = bisection(f,a,b,tol,df,ddf,Nmax)
    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))
    print('# of iterations =',count)




# define routines
def bisection(f,a,b,tol,df,ddf,Nmax):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 



    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    
    while (abs(1-(df(d)**2-f(d)*ddf(d))/df(d)**2)> 1):
      fd = f(d)
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    
    if ier==0:
          p = np.zeros(Nmax+1);
          p[0] = astar
          for it in range(Nmax):
              
              p1 = astar-f(astar)/df(astar)
              p[it+1] = p1
              count+=1
              if (abs(p1-astar) < tol):
                  pstar = p1
                  return [pstar,ier,count]
              astar = p1
          pstar = p1
          return [pstar,ier,count]

    return[astar,1,0]
      
driver()