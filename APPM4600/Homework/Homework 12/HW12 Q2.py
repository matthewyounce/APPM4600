import numpy as np

v=np.array([10-(116)**(1/2),4])
v=v/2
B=np.identity(2)-2*np.outer(v,v)/np.inner(v,v)
H=np.zeros([3,3])
H[0,0]=1
H[1:,1:]=B

A=np.array([[12,10,4],[10,8,-5],[4,-5,3]])
Ht=np.transpose(H)

print(H@A@Ht)


