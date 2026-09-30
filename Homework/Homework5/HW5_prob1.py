## APPM 4600 HW 5, Prob 1
# Cooper Wark

import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm

def driver():
    print(np.sqrt(3)/2)
    x0it = np.array([1.0,1.0])
    x0newt = np.array([1.0,1.0])
    Nmax = 100
    tol = 1e-12
    [xstar,nit,error] = numIt(x0it,tol,Nmax)
    print("x approx: ", xstar[0])
    print("y approx: ", xstar[1])
    print("Num iterations: ", nit)
    print("Error message: ", error)
    print()
    [xstar,ier,its] = Newton(x0newt,tol,Nmax)
    print("Newton x approx: ", xstar[0])
    print("Newton y approx: ", xstar[1])
    print("Num iterations: ", its)
    print("Error message: ", ier)


def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    for its in range(Nmax):
       J = evalJ(x0)
       Jinv = inv(J)
       F = evalF(x0)
       
       x1 = x0 - Jinv.dot(F)
       
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier, its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]

def evalF(x): 
# vector function that you want to find the roots of
    F = np.zeros(2)
    F[0] = f(x[0],x[1])
    F[1] = g(x[0],x[1])
    return F

def evalJ(x): 
# Jacobian of the vector function you want to find the roots of
    J = np.array([[6*x[0],-2*x[1]], 
        [3*x[1]**2-3*x[0]**2,6*x[0]*x[1]]])
    return J


def numIt(x,tol,Nmax):
    xn = np.zeros(2)
    nit = 0
    for i in range(Nmax):
        fx = f(x[0],x[1])
        gx = g(x[0],x[1])
        xn[0] = x[0] - (1/6)*fx - (1/18)*gx
        xn[1] = x[1] - (1/6)*gx
        if norm(xn-x) < tol:
            error = 0
            return [xn,nit,error]
        nit = nit + 1
        x[0] = xn[0]
        x[1] = xn[1]
    error = 1
    return [xn,nit, error]

def f(x,y):
    return 3*x**2 - y**2

def g(x,y):
    return 3*x*y**2 - x**3 - 1

driver()