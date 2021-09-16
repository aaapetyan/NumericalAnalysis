import math
import numpy as np
from numpy import *
import matplotlib.pyplot as plt

#eps = raw_input('eps = '); eps = float(eps)
eps = 0.000000001
n = raw_input('n = '); n = int(n)
f = ones(n); AA = zeros((n,n))
alp = 0.01
res = open('results', 'w+')

def poli_lezh(x,n):
      p0=1; p1=x
      for i in range(2,n+1):
            pi = 1.0/float(i) * ( (2.0*i - 1.0)*x*p1 - (i - 1.0)*p0 )
            p0 = p1
            p1 = pi
      pp = n / (1 - x**2) * ( p0 - x*p1)
      return pi, pp

def roots_lezh(n):
      t = []
      for i in range(1,n+1):
            tk = -2
            t0 = np.cos (np.pi * (4.0*i - 1.0) / (4.0 * n + 2) ) 
            while abs(tk - t0) > eps:
                  tk = t0 - poli_lezh(t0,n)[0] / poli_lezh(t0,n)[1]
                  t0 = tk
            t.append(tk)
      return t

def weights_lezh(n):
      a = []
      r = roots_lezh(n)
      for i in range(n):
            a.append ( 2.0 / (1 - r[i]**2) / poli_lezh(r[i],n)[1]**2 )
      return a

roots = roots_lezh(n)
weights = weights_lezh(n)


for i in range(n):
      roots[i] = 0.5*roots[i] + 0.5
      weights[i] = weights[i]/2.0
print roots
print weights

for i in range(n):
      for k in range(n):
            AA[i,k] = weights[k]*(1.0-alp*np.exp(-(roots[k]+roots[i])))/(1.0+alp*np.exp(-(roots[k]+roots[i])))
      AA[i,i] = AA[i,i] + 1.0

B = linalg.inv(AA)

u = dot (B,f)

for i in range(n):
      print>>res, roots[i], ' ', u[i]

plt.figure(1)
plt.subplot(111)
plt.plot(roots, u)
plt.xlim([0,1])
plt.xlabel('x')
plt.ylabel('u')
plt.show()
