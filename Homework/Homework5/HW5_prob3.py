## APPM4600 Homework 5, Problem 3
# Cooper Wark

import numpy as np
from numpy.linalg import norm

def driver():
    x0 = [1.0,1.0,1.0] # x,y,z initial guesses
    tol = 1e-14
    Nmax = 100

    [xstar,hist,ier,nit] = itScheme(x0,tol,Nmax)
    print("Sol Approx: ")
    print("x star: ", xstar[0])
    print("y star: ", xstar[1])
    print("z star: ", xstar[2])
    print("Error Message: ", ier)
    print("Num Its: ", nit)
    print()

    error(hist,xstar)

def itScheme(x0,tol,Nmax):
    x1 = np.zeros(3)
    hist = [x0.copy()]
    for nit in range(Nmax):
        df = evaldf(x0)
        x1[0] = x0[0] - df[0]
        x1[1] = x0[1] - df[1]
        x1[2] = x0[2] - df[2]
        hist.append(x1.copy())
        if norm(x1-x0) < tol:
            ier = 0
            return [x1,np.array(hist),ier,nit]
        # Update for next it
        x0[0] = x1[0]
        x0[1] = x1[1]
        x0[2] = x1[2]
    ier = 1
    return [x1,np.array(hist),ier,nit]

def evaldf(x):
    # Function: x^2 + 4y^2 + 4z^2 - 16 = 0
    f = x[0]**2 + 4*x[1]**2 + 4*x[2]**2 - 16

    # d = f/(fx^2 + fy^2 + fz^2)
    fx = 2*x[0] # 2x
    fy = 8*x[1] # 8y
    fz = 8*x[2] # 8z
    d = f/(fx**2 + fy**2 + fz**2)

    return [d*fx, d*fy, d*fz]

def error(history,xstar):
    for i in range(len(history)):
        print(norm(history[i,:]-xstar))

driver()