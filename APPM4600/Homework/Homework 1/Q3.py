# This script is for finding z in taylor's remainder formula 

# initializing libraries
import numpy as np

# intializing indepdent varible for part a
x1=np.linspace(0,.5,1000)

# intializing third derivative
y1=abs((1-17*x1+x1**3)*np.sin(x1)+(3-9*x1**2)*np.cos(x1))

#finding z and f(z)
fz1=y1.max()
z1=x1[y1.argmax()]

print('z = ',z1, ' f(z) = ',fz1) 

#repeating the process for part d where z can now range from 0 to 1
x2=np.linspace(0,1,1000)

# intializing third derivative
y2=abs((1-17*x2+x2**3)*np.sin(x2)+(3-9*x2**2)*np.cos(x2))

#finding z and f(z)
fz2=y2.max()
z2=x2[y2.argmax()]

print('z = ',z2, ' f(z) = ',fz2) 