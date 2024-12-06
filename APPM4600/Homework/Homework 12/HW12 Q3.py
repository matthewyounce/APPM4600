import numpy as np 
  
# Define the tolerance for the eigenvalue 
# and eigenvector approximations 
# (i.e. the maximum allowed difference between 
# the approximations and the actual values) 
#tol = 1e-6
tol=1e-16
  
# Define the maximum number of iterations 
max_iter = 10000
  
def power(A,x,tol,max_iter):
    lam_prev=1e99
    for i in range(max_iter): 
    # Compute the updated approximation for the eigenvector 
        x = A @ x / np.linalg.norm(A @ x) 
  
    # Compute the updated approximation for the largest eigenvalue 
        lam = (x.T @ A @ x) / (x.T @ x) 
  
    # Check if the approximations have converged 
        if np.abs(lam - lam_prev) < tol: 
            return lam,x,i+1,0
  
    # Store the current approximation for the largest eigenvalue 
        lam_prev = lam 
    return lam,x,max_iter,1

def Hilbert(n):
    A=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            A[i,j]=1/(i+j+1)
    return A

def create_x(n):
    return np.random.uniform(-5,5,n)

#Part A
#for n in range(4,21,4):
#    A=Hilbert(n)
#    v=create_x(n)
#    (lam,v,itr,err)=power(A,v,tol,max_iter)
#    print(lam,v,itr,err,n)    


#Part B
np.random.seed(1)
A=Hilbert(16)
v=create_x(16)
(lam,v,itr,err)=power(np.linalg.inv(A-1e-12*np.identity(16)),v,tol,max_iter)
#print(1e-12+1/lam,itr,err) 
#val,vec = np.linalg.eig(A)
#print(abs(1e-12+1/lam-min(abs(val))))
#print(abs(1e-12+1/lam-min(abs(val)))/(1e-12+1/lam))

#Part C
E=np.random.uniform(-1,1, (16, 16))
val,vec = np.linalg.eig(A+E)
U,S,V= np.linalg.svd(E)
print(abs(1e-12+1/lam-max(abs(val)))<S[0]/S[-1])


#Part D 
#having intial guess orthognoal to eigenvector of largest eigenvalue
#A=np.array([[5,0],[0,1]])
#x=np.array([0,1])
#(lam,v,itr,err)=power(A,x,tol,max_iter)
#print(lam,v,itr,err)  

  