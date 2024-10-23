import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from numpy.linalg import norm

def driver():


    f = lambda x: 1/(1+(10*x)**2)

    ''' interval'''
    a = -1
    b = 1
   
    for N in range(199,200):
       ''' create equispaced interpolation nodes'''
       # Question 2
       #xint = np.linspace(a,b,N+1)
       
       #Question 3
       xint=np.zeros(N+1)
       for i in range(N+1):
           xint[i]=np.cos((2*i-1)*np.pi/(2*N))
       xint=xint[1:]
    
       ''' create interpolation data'''
       yint = f(xint)
    
       ''' create points for evaluating the Lagrange interpolating polynomial'''
       Neval = 1000
       xeval = np.linspace(a,b,Neval+1)
       yeval_l= np.zeros(Neval+1)
  
       '''Initialize and populate the first columns of the 
       divided difference matrix. We will pass the x vector'''
       y = np.zeros( (N, N) )
     
       for j in range(N):
           y[j][0]  = yint[j]

       ''' evaluate lagrange poly '''
       for kk in range(Neval+1):
           yeval_l[kk] = eval_lagrange(xeval[kk],xint,yint,N)
          

    


       ''' create vector with exact values'''
       fex = f(xeval)
    
       err_l =  norm(fex-yeval_l) 
       print(f'langrange err for N={N} = ', err_l)
       

       plt.figure()    
       plt.plot(xeval,yeval_l,'o',color='red')
       plt.plot(xeval,fex) 
       plt.legend()

       plt.figure() 
       err_l = abs(yeval_l-fex)
       plt.semilogy(xeval,err_l)
       plt.show()

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N)
    
    for count in range(N):
       for jj in range(N):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)
  
driver()