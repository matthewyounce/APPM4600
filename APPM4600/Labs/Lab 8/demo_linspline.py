import matplotlib.pyplot as plt
import numpy as np


def driver():
    
    f = lambda x: 1/(1+(10*x)**2)
    a = -1
    b = 1
    
    ''' create points you want to evaluate at'''
    Neval = 100
    xeval =  np.linspace(a,b,Neval)
    
    for i in range(2,11):
        
        ''' number of intervals'''
        Nint = i
    
        '''evaluate the linear spline'''
        yeval = eval_lin_spline(xeval,Neval,a,b,f,Nint)
    
        ''' evaluate f at the evaluation points'''
        fex = np.zeros(Neval)
        for j in range(Neval):
            fex[j] = f(xeval[j]) 
      
        plt.figure()
        plt.plot(xeval,fex)
        plt.plot(xeval,yeval)
        plt.show()
        err = abs(yeval-fex)
        plt.figure()
        plt.semilogy(xeval,err)
        plt.show()
    
    
def  eval_lin_spline(xeval,Neval,a,b,f,Nint):

    '''create the intervals for piecewise approximations'''
    #part 1
    #xint = np.linspace(a,b,Nint+1)
    
    #part 2
    xint=np.zeros(Nint+1)
    for i in range(Nint+1):
        xint[i]=np.cos((2*(i+1)-1)*np.pi/(2*(Nint+1)))
    xint=xint[::-1]
        
    print(xint)
   
    '''create vector to store the evaluation of the linear splines'''
    yeval = np.zeros(Neval) 
    
    for jint in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        '''let n denote the length of ind'''
        
        '''temporarily store your info for creating a line in the interval of 
         interest'''
        a1= xint[jint]
        fa1 = f(a1)
        b1 = xint[jint+1]
        fb1 = f(b1)
        
        for kk in range(len(xeval)):
           '''use your line evaluator to evaluate the lines at each of the points 
           in the interval'''
           '''yeval(ind(kk)) = call your line evaluator at xeval(ind(kk)) with 
           the points (a1,fa1) and (b1,fb1)''' 
           print(xeval[kk],a1,b1)
           if xeval[kk]>=a1 and xeval[kk]<b1:
               yeval[kk]=create_line(a1,fa1,b1,fb1,xeval[kk])
               print(create_line(a1,fa1,b1,fb1,xeval[kk]))
               
                  
    return yeval


def create_line(a,fa,b,fb,alpha):
    return ((fb-fa)/(b-a))*(alpha-a)+fa
           
           
driver()               
