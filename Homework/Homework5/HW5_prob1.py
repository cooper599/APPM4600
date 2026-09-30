## APPM 4600 HW 5, Prob 1
# Cooper Wark

import numpy as np

def driver():
    x0 = 1
    y0 = 1
    Nmax = 100
    tol = 1e-12
    [xstar,ystar,nit,error] = numIt(x0,y0,tol,Nmax)
    print("x approx: ", xstar)
    print("y approx: ", ystar)
    print("Num iterations: ", nit)
    print("Error message: ", error)

def numIt(x,y,tol,Nmax):
    nit = 0
    for i in range(Nmax):
        fx = f(x,y)
        gx = g(x,y)
        xn = x - (1/6)*fx - (1/18)*gx
        yn = y - (1/6)*gx
        if (abs(xn-x) < tol) and (abs(yn-y) < tol):
            error = 0
            return [xn,yn,nit,error]
        nit = nit + 1
        x = xn
        y = yn
    error = 1
    return [xn, yn, nit, error]

def f(x,y):
    return 3*x**2 - y**2

def g(x,y):
    return 3*x*y**2 - x**3 - 1

driver()