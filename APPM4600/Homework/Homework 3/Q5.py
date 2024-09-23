#import libraries
import numpy as np
import matplotlib.pyplot as plt
    
def driver():
    
    
     #part A
     
     #creating my vectors
     x=np.linspace(-3,5,1000)
     y=x-4*np.sin(2*x)-3
     zero=np.zeros(1000)
     
     #plotting my vectors
     plt.plot(x,y)
     plt.plot(x,zero, 'red')
     plt.xlabel('x')
     plt.ylabel('y')
     plt.title('y=x+4sin(2x)-3')
     plt.show()

# test functions 
     f1 = lambda x: 5*x/4-np.sin(2*x)-3/4


     Nmax = 100
     tol = 1e-13

# test f1 '''
     x0 = -1
     [xstar1,ier,count] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar1)
     print('f1(xstar1):',f1(xstar1))
     print('Error message reads:',ier)
     
     x0 = -.5
     [xstar2,ier,count] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar2)
     print('f1(xstar2):',f1(xstar2))
     print('Error message reads:',ier)
     
     x0 = 1.5
     [xstar3,ier,count] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar3)
     print('f1(xstar3):',f1(xstar3))
     print('Error message reads:',ier)
     
     x0 = 3
     [xstar4,ier,count] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar4)
     print('f1(xstar4):',f1(xstar4))
     print('Error message reads:',ier)
     
     x0 = 4.5
     [xstar5,ier,count] = fixedpt(f1,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar5)
     print('f1(xstar5):',f1(xstar5))
     print('Error message reads:',ier)

# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier,count]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier,count]
    

driver()