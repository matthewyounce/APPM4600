#importing libraries
import math

#part b
#creating my function
f=lambda x: math.e**x-1


#testing the accuracy of the function
print('f(x) =',f(9.999999995*10**(-10)))

#veryifying
x=9.999999995*10**(-10)
print('T(x) =',x+x**2/2)