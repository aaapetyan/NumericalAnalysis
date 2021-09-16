import json
import math
import numpy as np
import matplotlib.pyplot as plt

def simple_iter(x,y,z,h):
      tmp = 4.0*(1.0-h/2)/(4.0 - 2.0*h - 4.0*h*AA - 3.0*AA*h**2.0)*(z + 5.0*h*y/4.0/(1.0-h/2) + 2.0*h*x/(2.0-h)/np.sqrt(x**2.0+1.0))
      zn = 0.0; zz = tmp
      while abs(zz-zn) > eps:
            zn = zz
            zz = tmp + 4.0*(1.0-h/2)/(4.0 - 2.0*h - 4.0*h*AA - 3.0*AA*h**2.0) * h * zz / np.sqrt(zz**2.0+1.0)
      return zz

eps = 0.000001
a = 0.0; b = 1.0
AA = raw_input('A = '); AA = int(AA)
n1 = raw_input('n1 = '); n1 = int(n1)
n = raw_input('n = '); n = int(n)
y = 1.0; z = -1.0; x = a

f = open('results', 'w+')

xx=[]; yy=[]; zz=[]
xx.append(x); yy.append(y); zz.append(z)

h = 10.0/(abs(AA)*n1)

if h > 0:
      for i in range(n1):
            x = x + h
            z = simple_iter(x,y,z,h)
            y = 1 / (1 - h / 2.0) * ( y + h * AA * z + h * x / np.sqrt (x**2.0 + 1.0) )
            xx.append(x); yy.append(y); zz.append(z)
            print>>f, xx[i], ' ', yy[i], ' ', zz[i]

h = (1.0 - 10.0/abs(AA) ) / n

if h > 0:
      for i in range(n1-1,n1+n):
            x = x + h
            z = simple_iter(x,y,z,h)
            y = 1 / (1 - h / 2.0) * ( y + h * AA * z + h * x / np.sqrt (x**2.0 + 1.0) )
            xx.append(x); yy.append(y); zz.append(z)
            print>>f, xx[i+1], ' ', yy[i+1], ' ', zz[i+1]

plt.figure(1)
plt.subplot(211)
plt.plot(xx,yy)
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0,1)
plt.subplot(212)
plt.plot(xx,zz)
plt.xlabel('x')
plt.ylabel('z')
plt.xlim(0,1)
plt.show()
