## APPM 4600 Homework 2, Prob 6
# Cooper Wark
import numpy as np
import matplotlib.pyplot as plt

## part a
f = lambda x: x-4*np.sin(2*x)-3
eps = 0.5*10**-10
# Covers all roots
a = -2
b = 5

numpts = 100
x = np.linspace(a,b,numpts)

# plt.plot(x,f(x))
# plt.axhline(0, color="black", linewidth = 1)
# plt.axvline(0, color="black", linewidth = 1)
# plt.xlabel("x")
# plt.ylabel("y")
# plt.title("Plot of x-4sin(2x)-3 showing all roots")
# plt.axis("on")
# plt.show()

## Part b
# From code provided in class
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''
    x_arr = np.zeros((Nmax,1))

    count = 0
    while (count <Nmax):
       x_arr[count] = x0
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier, x_arr, count]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier, x_arr, count]

# other stuff for part b
Nmax = 200
tol = 10**-12 # make sure at least 10 accurate digits

g = lambda x: -np.sin(2*x)+5*x/4-3/4
x0 = 0
x0_arr = [-1,-0.89, 0, 2, 3, 4]

for i in range(len(x0_arr)):
    [xstar,ier, x_arr, count] = fixedpt(g,x0_arr[i],tol,Nmax)
    print("For x0 = ", x0_arr[i])
    print("xstar = ", xstar)
    print("ier = ", ier)
    print("count = ", count)
    print("")
