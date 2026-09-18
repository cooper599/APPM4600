## APPM 4600 Hw 3 prob 3
# Cooper Wark

import numpy as np

def driver():
    fa = lambda x: x*(1+(7-x**5)/x**2)**3
    fb = lambda x: x-(x**5-7)/x**2
    fc = lambda x: x - (x**5-7)/(5*x**4)
    fd = lambda x: x-(x**5-7)/12

    x = 7**(1/5)

    print("x = ", x)
    print("fa(x) = ", fa(x))
    print("fb(x) = ", fb(x))
    print("fc(x) = ", fc(x))
    print("fd(x) = ", fd(x))

    x0 = 1
    tol = 1e-10 #10**-10 same
    Nmax = 1000
    # [astar, iera] = fixedpt(fa,x0,tol,Nmax) # Overflow error
    # [bstar, ierb] = fixedpt(fb,x0,tol,Nmax)
    [cstar, ierc] = fixedpt(fc,x0,tol,Nmax)
    [dstar, ierd] = fixedpt(fd,x0,tol,Nmax)

    # print("Root approx of fa = ", astar, "Error Message: ", iera)
    # print("Root approx of fb = ", bstar, "Error Message: ", ierb)
    print("Root approx of fc = ", cstar, "Error Message: ", ierc)
    print("Root approx of fd = ", dstar, "Error Message: ", ierd)

def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]

driver()