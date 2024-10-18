import numpy as np
import matplotlib.pyplot as plt
def create_line(a,fa,b,fb,alpha):
    return ((fb-fa)/(b-a))*(alpha-a)+fa

a=0
b=1
fa=0
fb=2
y=np.zeros(11)
for i,j in enumerate(np.linspace(0,1,11)):
    y[i]=create_line(a,fa,b,fb,j)
print(y)

plt.plot(np.linspace(0,1,11),y)
    
    

    