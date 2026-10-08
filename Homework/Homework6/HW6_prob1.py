## APPM4600 Homework 6 Prob 1
# Cooper Wark

import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm

def driver():
    tol = 1e-4
    Nmax = 100

    x01 = [-0.5,0.25]
    x02 = [-1.5, 1.25]
    x03 = [-3.0,3.0]

    print("Fixed Its: ")
    [xstar1,ier1,its1] = fixedIt(x01,tol,Nmax)
    [xstar2,ier2,its2] = fixedIt(x02,tol,Nmax)
    [xstar3,ier3,its3] = fixedIt(x03,tol,Nmax)
    printRes(xstar1,ier1,its1)
    printRes(xstar2,ier2,its2)
    printRes(xstar3,ier3,its3)

    print("Newton: ")
    [xstar1,ier1,its1] = Newton(x01,tol,Nmax)
    [xstar2,ier2,its2] = Newton(x02,tol,Nmax)
    [xstar3,ier3,its3] = Newton(x03,tol,Nmax)
    printRes(xstar1,ier1,its1)
    printRes(xstar2,ier2,its2)
    printRes(xstar3,ier3,its3)

def printRes(xstar,ier,its):
    print("x approx: ", xstar[0])
    print("y approx: ", xstar[1])
    print("Error: ", ier)
    print("Num Its: ", its)
    print()

def fixedIt(x0,tol,Nmax):
    A = np.array([[0.016,-0.17],
         [0.52,-0.26]])
    for its in range(Nmax):
        F = evalF(x0)
        x1 = x0 - A.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier = 0
            return [xstar,ier,its]
        x0 = x1
    xstar = x1
    ier = 1
    return [xstar,ier,its]

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
    F[0] = 3*x[0]**2 + 4*x[1]**2 - 1
    F[1] = x[1]**3 - 8*x[0]**3 -1
    return F

def evalJ(x): 
# Jacobian of the vector function you want to find the roots of
    J = np.array([[6*x[0],8*x[1]], 
        [-24*x[0]**2,3*x[1]**2]])
    return J

driver()