import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
import time 



def driver():

     ''' create  matrix for testing different ways of solving a square 
     linear system'''

     '''' N = size of system'''
    
     for N in [100, 500, 1000, 2000, 4000, 5000]:
 
         ''' Right hand side'''
         b = np.random.rand(N,1)
         A = np.random.rand(N,N)
         t=time.time()
         x = scila.solve(A,b)
         norm_t=time.time()-t
         t=time.time()
         (lu,piv)=scila.lu_factor(A)
         lu_time=time.time()-t
         t=time.time()
         v= scila.lu_solve((lu,piv),b)
         lu_solve_t=time.time()-t
         
         print('N=',N,' regular time=',norm_t,'lu decomp time=',lu_time,'lu solve time=',lu_solve_t)
     
         test = np.matmul(A,x)
         lu_test=np.matmul(A,v)
         r = la.norm(test-b)
         vr=la.norm(lu_test-b)
     
         #print(r)
         #print(vr)

     ''' Create an ill-conditioned rectangular matrix '''
     N = 10
     M = 5
     A = create_rect(N,M)     
     b = np.random.rand(N,1)


     
def create_rect(N,M):
     ''' this subroutine creates an ill-conditioned rectangular matrix'''
     a = np.linspace(1,10,M)
     d = 10**(-a)
     
     D2 = np.zeros((N,M))
     for j in range(0,M):
        D2[j,j] = d[j]
     
     '''' create matrices needed to manufacture the low rank matrix'''
     A = np.random.rand(N,N)
     Q1, R = la.qr(A)
     test = np.matmul(Q1,R)
     A =    np.random.rand(M,M)
     Q2,R = la.qr(A)
     test = np.matmul(Q2,R)
     
     B = np.matmul(Q1,D2)
     B = np.matmul(B,Q2)
     return B     
          
  
if __name__ == '__main__':
      # run the drivers only if this is called from the command line
      driver()       
