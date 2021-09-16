import numpy as np
from numpy import zeros

n = raw_input('n = '); n = int(n)
M = raw_input('M = '); M = int(M)
f = open('results.dat','w+')
a = 0.0; b = 1.0; T = 1.0
h = (b - a) / n; tau = T / M; sigma = tau / h**2

def fi(x):
      return 1.0 / (1.0 + x**2)**2

u_minus = zeros(n+1); u = zeros(n+1)

for i in range(n+1):
      u_minus[i] = fi(a+i*h)

for k in range(1,M+1):
      alpha = [0]; beta = [0]
      alpha.append( 2.0*sigma / (2.0*sigma + 3.0*tau + 3.0) )
      beta.append( 3.0*u_minus[1] / (2.0*sigma + 3.0*tau + 3.0) )
      for i in range(2,n):
            alpha.append( sigma / ( 2.0*sigma + tau + 1.0 - sigma * alpha[i-1] ) )
            beta.append ( (sigma*beta[i-1] + u_minus[i]) / ( 2.0*sigma + tau + 1.0 - sigma * alpha[i-1] ) )
      u[n] = ((4.0-alpha[n-2])*beta[n-1]-beta[n-2])/(3.0+alpha[n-1]*(alpha[n-2]-4.0))
      for i in range(n-1,0,-1):
            u[i] = alpha[i]*u[i+1]+beta[i]
      u[0] = (4.0*u[1]-u[2]) / 3.0
      
      for i in range(n+1):
            u_minus[i] = u[i]

      if ((k+1)%(M/10)==0):
            print>>f, k+1,'    ',u
