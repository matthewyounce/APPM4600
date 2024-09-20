#Importing libraries
import numpy as np
import matplotlib.pyplot as plt

#Part A

#Creating my vectors
t=np.linspace(0,np.pi,31)
y=np.cos(t)

#Finding my sum
S=0
for i in range(len(t)):
    S+=t[i]*y[i]

#printing my awnser
print('The sum is: S', S)

#Part B

#Creating theta
t=np.linspace(0,2*np.pi,1000)

#Creating my functions
x=1.2*(1+.1*np.sin(15*t))*np.cos(t)
y=1.2*(1+.1*np.sin(15*t))*np.sin(t)

#ploting my functions
plt.plot(x,y)
plt.title('Fixed Values')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('square')
plt.show()

#Creating my functions with a foor loop
for i in range(1,11):
    x=i*(1+.05*np.sin((2+i)*t+np.random.uniform(0,2)))*np.cos(t)
    y=i*(1+.05*np.sin((2+i)*t+np.random.uniform(0,2)))*np.sin(t)
    plt.plot(x,y)

plt.title('For loop')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('square')
plt.show()