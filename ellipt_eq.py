import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

n = int(raw_input('n= '))
omega = float(raw_input('omega='))
eps = 0.000001
h = 1.0 / n
u = np.zeros((n+1,n+1)); u_next = np.zeros((n+1,n+1))
A = np.zeros((n+1,n+1)); B = np.zeros((n+1,n+1))
C = np.zeros((n+1,n+1)); D = np.zeros((n+1,n+1))
E = np.zeros((n+1,n+1))

for i in range(n+1):
      for k in range(n+1):
            x = i*h; y = k*h
            A[i][k] = (1.0+x**2)/h**2 - x/h; D[i][k] = (1.0+x**2)/h**2 + x/h
            B[i][k] = math.log(4.0+y) / h**2 - 1.0/((4.0+y)*2.0*h)
            C[i][k] = math.log(4.0+y) / h**2 + 1.0/((4.0+y)*2.0*h)
            E[i][k] = 2.0/h**2*(1.0+x**2+math.log(4.0+y))

for i in range(n+1):
      A[0][i] = 0; D[0][i] = 0
      B[0][i] = math.log(4.0+(i*h)) / h**2 - 1.0/((4.0+(i*h))*2.0*h)
      C[0][i] = math.log(4.0+(i*h)) / h**2 + 1.0/((4.0+(i*h))*2.0*h)
      E[0][i] = 2.0/h**2*math.log(4.0+(i*h))

A[n] = A[0]; B[n] = B[0]; C[n] = C[0]; D[n] = D[0]; E[n] = E[0]

for i in range(n+1):
      u_next[i][0] = 1.0+4.0*(i*h)*(1-i*h)
      u_next[i][n] = 1.0+4.0*(i*h)*(1-i*h)
            
def meth_iter(curr,next):
      n_iter = 0
      u = np.array(curr); u_next = np.array(next)
      while (np.max(abs(u_next - u)) > eps):
            u = [[u_next[i][k] for k in range(n+1)] for i in range(n+1)]
            for i in range(n+1):
                  u_next[i][0] = 1.0+4.0*(i*h)*(1.0-i*h)
            for i in range(1,n):
                  for k in range(1,n):
                        u_next[i][k] = 1.0/E[i][k] * (A[i][k]*u[i-1][k]+B[i][k]*u[i][k-1]+C[i][k]*u[i][k+1]+D[i][k]*u[i+1][k])
            for i in range(n+1):
                  u_next[i][n] = 1.0+4.0*(i*h)*(1.0-i*h)
            n_iter = n_iter + 1
      print 'number of iterations = ',n_iter
      return u_next

def meth_zeidel(curr,next):
      n_iter = 0
      u = np.array(curr); u_next = np.array(next)
      while (np.max(abs(u_next - u)) > eps):
            u = [[u_next[i][k] for k in range(n+1)] for i in range(n+1)]
            for i in range(n+1):
                  u_next[i][0] = 1.0+4.0*(i*h)*(1-i*h)
            for i in range(1,n):
                  for k in range(1,n):
                        u_next[i][k] = 1.0/E[i][k] * (A[i][k]*u_next[i-1][k]+B[i][k]*u_next[i][k-1]+C[i][k]*u[i][k+1]+D[i][k]*u[i+1][k])
            for i in range(n+1):
                  u_next[i][n] = 1.0+4.0*(i*h)*(1-i*h)
            n_iter = n_iter + 1
      print 'number of iterations = ',n_iter
      return u_next

def meth_relax(curr,next,omega,return_iter=False):
      n_iter = 0
      u = np.array(curr); u_next = np.array(next)
      while (np.max(abs(u_next - u)) > eps):
            u = [[u_next[i][k] for k in range(n+1)] for i in range(n+1)]
            for i in range(n+1):
                  u_next[i][0] = 1.0+4.0*(i*h)*(1-i*h)
            for i in range(1,n):
                  for k in range(1,n):
                        u_rel = 1.0/E[i][k] * (A[i][k]*u_next[i-1][k]+B[i][k]*u_next[i][k-1]+C[i][k]*u[i][k+1]+D[i][k]*u[i+1][k])
                        u_next[i][k] = u[i][k] + omega * (u_rel - u[i][k])
            for i in range(n+1):
                  u_next[i][n] = 1.0+4.0*(i*h)*(1-i*h)
            n_iter = n_iter + 1
      if return_iter:
            return n_iter
      else:
            print 'number of iterations = ',n_iter
            return u_next


one = meth_iter(u, u_next)
print one
print
two = meth_zeidel(u,u_next)
print two
print

three = meth_relax(u, u_next, omega)
print 'omega=',omega
print three

def plot_u():
      fig = plt.figure(figsize=plt.figaspect(0.3))
      #iter
      ax = fig.add_subplot(1, 3, 1, projection='3d')
      X = np.arange(0,n+1,1)
      Y = np.arange(0,n+1,1)
      X, Y = np.meshgrid(X, Y)
      Z = one

      surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
      #ax.set_zlim(-1.01, 1.01)

      ax.zaxis.set_major_locator(LinearLocator(10))
      ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

      #fig.colorbar(surf, shrink=0.5, aspect=5)
      
      #zeidel
      ax = fig.add_subplot(1, 3, 2, projection='3d')
      X = np.arange(0,n+1,1)
      Y = np.arange(0,n+1,1)
      X, Y = np.meshgrid(X, Y)
      Z = two

      surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
      #ax.set_zlim(-1.01, 1.01)

      ax.zaxis.set_major_locator(LinearLocator(10))
      ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

      #fig.colorbar(surf, shrink=0.5, aspect=5)

      #relax
      ax = fig.add_subplot(1, 3, 3, projection='3d')
      X = np.arange(0,n+1,1)
      Y = np.arange(0,n+1,1)
      X, Y = np.meshgrid(X, Y)
      Z = three

      surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
      #ax.set_zlim(-1.01, 1.01)

      ax.zaxis.set_major_locator(LinearLocator(10))
      ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

      #fig.colorbar(surf, shrink=0.5, aspect=5)

      plt.show()

plot_u()
