# Intializing my libraries
import numpy as np
import matplotlib.pyplot as plt

# Intializing my given x values
x1=np.pi
x2=10**6

# Creating an list of my delta values
delta=[10**i for i in range(-16,1)]

# Creating my lists of differnces for each x
y1=[-2*np.sin((2*x1+d)/2)*np.sin(d/2) for d in delta]
y2=[-2*np.sin((2*x2+d)/2)*np.sin(d/2) for d in delta]

# Plotting each of the differences for x1
plt.semilogx(delta,y1)
plt.xlabel('Delta')
plt.ylabel('Difference')
plt.title('Difference for x=pi')
plt.show()

# Plotting each of the differences for x2
plt.semilogx(delta,y2)
plt.xlabel('Delta')
plt.ylabel('Difference')
plt.title('Difference for x=10^6')
plt.show()


# Creating my algorthim for part c for both x values
y3=[-d*np.sin(x1)-d**2/2*np.cos(x1) for d in delta]
y4=[-d*np.sin(x2)-d**2/2*np.cos(x2) for d in delta]

# Plotting the comparison for x1
plt.semilogx(delta,y3)
plt.xlabel('Delta')
plt.ylabel('Difference')
plt.title('Taylor Expansion x=pi')
plt.show()

# Plotting the comparison for x2
plt.semilogx(delta,y4)
plt.xlabel('Delta')
plt.ylabel('Difference')
plt.title('Taylor Expansion for x=10^6')
plt.show()


