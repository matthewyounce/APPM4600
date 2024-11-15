import numpy as np

def lgwt(N,a,b):
  """ 
   This script is for computing definite integrals using Legendre-Gauss 
   Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
   [a,b] with truncation order N
  
   Suppose you have a continuous function f(x) which is defined on [a,b]
   which you can evaluate at any x in [a,b]. Simply evaluate it at all of
   the values contained in the x vector to obtain a vector f. Then compute
   the definite integral using np.sum(f*w)
  
   Written by Greg von Winckel - 02/25/2004
   translated to Python - 10/30/2022
  """
  N = N-1
  N1 = N+1
  N2 = N+2
  eps = np.finfo(float).eps  
  xu = np.linspace(-1,1,N1)
  
  # Initial guess
  y = np.cos((2*np.arange(0,N1)+1)*np.pi/(2*N+2))+(0.27/N1)*np.sin(np.pi*xu*N/N2)

  # Legendre-Gauss Vandermonde Matrix
  L = np.zeros((N1,N2))
  
  # Compute the zeros of the N+1 Legendre Polynomial
  # using the recursion relation and the Newton-Raphson method
  
  y0 = 2.
  one = np.ones((N1,))
  zero = np.zeros((N1,))

  # Iterate until new points are uniformly within epsilon of old points
  while np.max(np.abs(y-y0)) > eps:
      
    L[:,0] = one
    
    L[:,1] = y
    for k in range(2,N1+1): 
      L[:,k] = ((2*k-1)*y*L[:,k-1]-(k-1)*L[:,k-2])/k
    
    lp = N2*(L[:,N1-1]-y*L[:,N2-1])/(1-y**2)   
    
    y0 = y
    y = y0-L[:,N2-1]/lp
    
    
  
  # Linear map from[-1,1] to [a,b]
  x=(a*(1-y)+b*(1+y))/2
  
  # Compute the weights
  w=(b-a)/((1-y**2)*lp**2)*(N2/N1)**2
  return x,w

def eval_composite_trap(n,a,b,f):
  h = (b-a)/n
  w=[1/2]
  xnode = a+np.arange(0,n+1)*h
  
  I_trap = h*f(xnode[0])*1/2
  
  for j in range(1,n):
       I_trap = I_trap+h*f(xnode[j])
  I_trap= I_trap + 1/2*h*f(xnode[n])
  w.append(1)
  w.append(w)
  
  return I_trap,xnode,w

def eval_composite_simpsons(n,a,b,f):
  h = (b-a)/n
  xnode = a+np.arange(0,n+1)*h
  I_simp = f(xnode[0])
  w=[1/3]

  nhalf = n/2
  for j in range(1,int(nhalf)+1):
       # even part 
       I_simp = I_simp+2*f(xnode[2*j])
       w.append(4/3)
       # odd part
       I_simp = I_simp +4*f(xnode[2*j-1])
       w.append(2/3)
  I_simp= I_simp + f(xnode[n])
  
  I_simp = h/3*I_simp
  w.append(1/3)
  
  return I_simp,xnode,w   

def eval_gauss_quad(M,a,b,f):
  """
  Non-adaptive numerical integrator for \int_a^b f(x)w(x)dx
  Input:
    M - number of quadrature nodes
    a,b - interval [a,b]
    f - function to integrate
  
  Output:
    I_hat - approx integral
    x - quadrature nodes
    w - quadrature weights

  Currently uses Gauss-Legendre rule
  """
  x,w = lgwt(M,a,b)
  I_hat = np.sum(f(x)*w)
  return I_hat,x,w

def adaptive_quad(a,b,f,tol,M,method):
  """
  Adaptive numerical integrator for \int_a^b f(x)dx
  
  Input:
  a,b - interval [a,b]
  f - function to integrate
  tol - absolute accuracy goal
  M - number of quadrature nodes per bisected interval
  method - function handle for integrating on subinterval
         - eg) eval_gauss_quad, eval_composite_simpsons etc.
  
  Output: I - the approximate integral
          X - final adapted grid nodes
          nsplit - number of interval splits
  """
  # 1/2^50 ~ 1e-15
  maxit = 50
  left_p = np.zeros((maxit,))
  right_p = np.zeros((maxit,))
  s = np.zeros((maxit,1))
  left_p[0] = a; right_p[0] = b;
  # initial approx and grid
  s[0],x,_ = method(M,a,b,f);
  # save grid
  X = []
  X.append(x)
  j = 1;
  I = 0;
  nsplit = 1;
  while j < maxit:
    # get midpoint to split interval into left and right
    c = 0.5*(left_p[j-1]+right_p[j-1]);
    # compute integral on left and right spilt intervals
    s1,x,_ = method(M,left_p[j-1],c,f); X.append(x)
    s2,x,_ = method(M,c,right_p[j-1],f); X.append(x)
    if np.max(np.abs(s1+s2-s[j-1])) > tol:
      left_p[j] = left_p[j-1]
      right_p[j] = 0.5*(left_p[j-1]+right_p[j-1])
      s[j] = s1
      left_p[j-1] = 0.5*(left_p[j-1]+right_p[j-1])
      s[j-1] = s2
      j = j+1
      nsplit = nsplit+1
    else:
      I = I+s1+s2
      j = j-1
      if j == 0:
        j = maxit
  return I,np.unique(X),nsplit

def driver():
    f=lambda x:np.sin(1/x)
    a=.1
    b=2
    M=5
    tol=10**(-3)
    [I,X,nsplit]=adaptive_quad(a,b,f,tol,M,eval_composite_trap)
    print('Trap',I,nsplit)
    [I,X,nsplit]=adaptive_quad(a,b,f,tol,M,eval_composite_simpsons)
    print('simpsons',I,nsplit)
    [I,X,nsplit]=adaptive_quad(a,b,f,tol,M,eval_gauss_quad)
    print('guass',I,nsplit)
    
    
driver()    
