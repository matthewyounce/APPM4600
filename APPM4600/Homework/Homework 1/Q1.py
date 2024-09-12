#initialize libraries 
import numpy as np
import matplotlib.pyplot as plt

# initializing my independent varaible
x=np.linspace(1.92,2.08,161)

# initializing my response varaibles for part i and ii
y1=x**9-18*x**8+144*x**7-672*x**6+2016*x**5-4032*x**4+5376*x**3-4608*x**2+2304*x-512
y2=(x-2)**9

#plotting y1
plt.plot(x,y1)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Expanded Polynomial')
plt.show()

#plotting y2
plt.plot(x,y2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Factored Polynomial')
plt.show()

#Comparing the difference between y1 and y2
y3=y1-y2

#Plotting the difference between y1 and y2
plt.plot(x,y3)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Difference in Polynomial Forms')
plt.show()
