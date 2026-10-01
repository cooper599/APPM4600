## APPM 4600 Lab 6
# Cooper Wark

import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm
import time

def driver():
    prelab = False
    lab = True
    if prelab:
        tol = 1e-6
        Nmax = 1000
        x01 = [2.0,0.5]
        x02 = [3.0,5.0]
        # (1.5,1.5) -> (1,1)
        x03 = [1.5,2.5]

        # Newton calls
        [xstarN1,ierN1,itsN1] = Newton(x01,tol,Nmax)
        [xstarN2,ierN2,itsN2] = Newton(x02,tol,Nmax)
        # Lazy Newton
        [xstarLN1,ierLN1,itsLN1] = LazyNewton(x01,tol,Nmax)
        [xstarLN2,ierLN2,itsLN2] = LazyNewton(x02,tol,Nmax)

        [xstarN3,ierN3,itsN3] = Newton(x03,tol,Nmax)
        [xstarLN3,ierLN3,itsLN3] = LazyNewton(x03,tol,Nmax)

        print("x0_1 (2.0,0.5) Outputs")
        print("Newton:")
        print("xstar: ", xstarN1)
        print("Error message: ", ierN1)
        print("Num Iterations: ", itsN1)
        print("Lazy Newton:")
        print("xstar: ", xstarLN1)
        print("Error message: ", ierLN1)
        print("Num Iterations: ", itsLN1)
        print()

        print("x0_2 (3,5) Outputs")
        print("Newton:")
        print("xstar: ", xstarN2)
        print("Error message: ", ierN2)
        print("Num Iterations: ", itsN2)
        print("Lazy Newton:")
        print("xstar: ", xstarLN2)
        print("Error message: ", ierLN2)
        print("Num Iterations: ", itsLN2)
        print()

        print("x0_3 (1.5,2.5) Outputs")
        print("Newton:")
        print("xstar: ", xstarN3)
        print("Error message: ", ierN3)
        print("Num Iterations: ", itsN3)
        print("Lazy Newton:")
        print("xstar: ", xstarLN3)
        print("Error message: ", ierLN3)
        print("Num Iterations: ", itsLN3)
    if lab:
        x0 = [1.0,0.0]
        tol = 1e-10
        Nmax = 100

        ts = time.perf_counter()
        [xstar,ier,its] = Newton(x0,tol,Nmax)
        te = time.perf_counter()
        print("Newton:")
        print("xstar: ", xstar)
        print("Error message: ", ier)
        print("Num iterations: ", its)
        print("Time: ", te-ts)
        print()

        ts2 = time.perf_counter()
        for i in range(60):
            [xstar,ier,its] = LazyNewton(x0,tol,Nmax)
        te2 = time.perf_counter()
        print("Lazy Newton:")
        print("xstar: ", xstar)
        print("Error message: ", ier)
        print("Num iterations: ", its)
        print("Time: ", (te2-ts2)/60)
        print()

        ts3 = time.perf_counter()
        for i in range(60):
            [xstar,ier,its] = SlackerNewton(x0,tol,Nmax)
        te3 = time.perf_counter()
        print("Slacker Newton:")
        print("xstar: ", xstar)
        print("Error message: ", ier)
        print("Num iterations: ", its)
        print("Time: ", (te3-ts3)/60)

# Function for system
def evalF(x):
    # input x - x1,x2
    # Prelab function
    F = np.zeros(2)
    # F[0] = x[0]**2 + x[1]**2 - 2
    # F[1] = np.e**(x[0]-1) + x[1]**2 - 2
    # In lab function:
    F[0] = 4*x[0]**2 + x[1]**2 - 4
    F[1] = x[0] + x[1] - np.sin(x[0]-x[1])
    return F

def evalJ(x): 
# Jacobian of the vector function you want to find the roots of
    # Prelba jacobina
    # J = np.array([[2*x[0],2*x[1]], 
    #     [(x[0]-1)*np.e**(x[0]-1),2*x[1]]])
    # In lab jacobian
    J = np.array([[8*x[0],2*x[1]], 
            [1-np.cos(x[0]-x[1]),1+np.cos(x[0]-x[1])]])
    return J

def SlackerNewton(x0,tol,Nmax):
    ''' Slacker Newton = update Jacobian only after certain "event" occurs'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):
       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       if (norm(x1-x0) < tol):
           xstar = x1
           ier = 0
           return[xstar, ier,its]
       # Gets exponent of tolerance, updates everytime the difference is greater than half of tolerance
       # tol 1e-10, check if difference greater than 1e-5
       if (norm(x1-x0) > 10**(np.floor(np.log10(tol))/2)):
           J = evalJ(x1)
           Jinv = inv(J)
       x0 = x1   
    xstar = x1
    ier = 1
    return[xstar,ier,its] 

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

def LazyNewton(x0,tol,Nmax):
    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):
       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its] 
       x0 = x1   
    xstar = x1
    ier = 1
    return[xstar,ier,its]  

driver()