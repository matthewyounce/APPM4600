import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad


def driver():
    
    f = lambda x: 1/(1+x**2)
    a = -5
    b = 5
    
    # exact integral
    I_ex = 2.7468
    
#    N =100
#    ntest = np.arrange(0,N,step=2)
    
#    errorT = np.zeros(len(ntest))
#    errorS = np.zeros(len(ntest))
    
#    for j in range(0,len(ntest)):
#        n = ntest[j]

# for simpson's n must be even.        
# n+1 = number of pts.
    n = 1291
    I_trap = CompTrap(a,b,n,f)
    print('I_trap= ', I_trap)
    
    err = abs(I_ex-I_trap)   
    
    print('absolute error = ', err)    
    n=110
    I_simp = CompSimp(a,b,n,f)

    print('I_simp= ', I_simp)
    
    err = abs(I_ex-I_simp)   
    
    print('absolute error = ', err)
    
    [I,error,d]=quad(f,-5,5,epsabs=10**-4,full_output=1)
    print('scipy 10^-4 = ',I)
    print('absolute error =',abs(I_ex-I))
    print('func calls =',d['neval'])
    [I,error,d]=quad(f,-5,5,epsabs=10 **-6,full_output=1)
    print('scipy 10^-6 = ',I)
    print('absolute error =',abs(I_ex-I))
    print('func calls =',d['neval'])
        

     
def CompTrap(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    
    I_trap = h*f(xnode[0])*1/2
    
    for j in range(1,n):
         I_trap = I_trap+h*f(xnode[j])
    I_trap= I_trap + 1/2*h*f(xnode[n])
    
    return I_trap     

def CompSimp(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    I_simp = f(xnode[0])

    for j in range(1,n):
         # even part
         if j%2==0:
             I_simp = I_simp+2*f(xnode[j])
         # odd part
         else:
             I_simp = I_simp +4*f(xnode[j])
    I_simp= I_simp + f(xnode[n])
    
    I_simp = h/3*I_simp
    
    return I_simp     


    
    
driver()    

#checking to find the largest value
f = lambda x:(2*(1+x**2)**2-8*x**2*(1+x**2))/(1+x**2)**4
x=np.linspace(-5,5,1000000)
y=abs(f(x))
print(max(y))