import math
import numpy as np
from numpy import *

n = raw_input('n = '); n = int(n)
M = raw_input('M = '); M = int(M)
f = open('results1.dat','w+')

a = 0.0; b = 1.0; T = 1.0
h = (b - a) / n; tau = T / M; sigma = tau / h**2

if (sigma > 0.5) is True:
      print 'sigma =',sigma,' -- disturbed stability condition warning'

def fi(x):
      return 1.0 / (1.0 + x**2)**2

uk = zeros(n+1); ukk = zeros(n+1)

for i in range(n+1):
      uk[i] = fi(a+i*h)

for k in range(M):
      for i in range(1,n):
            ukk[i] = sigma * uk[i-1] + (1 - 2*sigma-tau) * uk[i] + sigma * uk[i+1]
      ukk[0] = (4.0*ukk[1] - ukk[2]) / 3.0
      ukk[n] = (4.0*ukk[n-1] - ukk[n-2]) / 3.0

      for i in range(n+1):
            uk[i] = ukk[i]

      if ((k+1)%(M/10)==0):
            print>>f, k+1,'    ',uk
