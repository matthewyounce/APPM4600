import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(0,5,1000)
y=np.sin(x)

approx1=(x-7/60*x**3)/(1+1/20*x**2)
approx2=x/(1+x**2/6+7*x**4/360)
approx3=(x-7*x**3/60)/(1+x**2/20)
mac=x-x**3/6+x**5/120

error1=abs(y-approx1)
error2=abs(y-approx2)
error3=abs(y-approx3)
macerror=abs(y-mac)
plt.plot(x,error1,label='[3,3] pade',linewidth=4)
plt.plot(x,error2,label='[2,4] pade')
plt.plot(x,error3,'black', label='[4,2] pade')
plt.plot(x,macerror,label='maclaurin')
plt.title('Error')
plt.legend()
plt.show()