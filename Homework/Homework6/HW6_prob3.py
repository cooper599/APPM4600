## APPM4600 Homework 6 Prob 3
# Cooper Wark

import numpy as np
from numpy.linalg import norm
from numpy.linalg import inv

def driver():
    tol = 1e-6
    tol2 = 5e-2 # 5*10^-2
    Nmax = 100
    # x0 = [0.0, 0.0, 0.0]
    x0 = [0.5, 0.5, 0.5]

    [xstar,ier,its] = Newton(x0,tol,Nmax)
    print("Newton: ")
    printRes(xstar,ier,its)

    [x,g1,ier,its] = SteepestDescent(x0,tol,Nmax)
    print("Steepest Descent: ")
    printRes(x,ier,its)

    print("Combined Steepest + Newton: ")
    [x,g1,ier,its1] = SteepestDescent(x0,tol2,Nmax)
    [xstar,ier,its2] = Newton(x,tol,Nmax)
    printRes(xstar,ier,its1+its2)

def printRes(xstar,ier,its):
    print("x approx: ", xstar[0])
    print("y approx: ", xstar[1])
    print("z approx: ", xstar[2])
    print("Error: ", ier)
    print("Num Its: ", its)
    print()

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

def SteepestDescent(x,tol,Nmax):
    
    for its in range(Nmax):
        g1 = evalg(x)
        z = eval_gradg(x)
        z0 = norm(z)

        if z0 == 0:
            print("zero gradient")
        z = z/z0
        alpha1 = 0
        alpha3 = 1
        dif_vec = x - alpha3*z
        g3 = evalg(dif_vec)

        while g3>=g1:
            alpha3 = alpha3/2
            dif_vec = x - alpha3*z
            g3 = evalg(dif_vec)
            
        if alpha3<tol:
            print("no likely improvement")
            ier = 0
            return [x,g1,ier,its]
        
        alpha2 = alpha3/2
        dif_vec = x - alpha2*z
        g2 = evalg(dif_vec)

        h1 = (g2 - g1)/alpha2
        h2 = (g3-g2)/(alpha3-alpha2)
        h3 = (h2-h1)/alpha3

        alpha0 = 0.5*(alpha2 - h1/h3)
        dif_vec = x - alpha0*z
        g0 = evalg(dif_vec)

        if g0<=g3:
            alpha = alpha0
            gval = g0

        else:
            alpha = alpha3
            gval =g3

        x = x - alpha*z

        if abs(gval - g1)<tol:
            ier = 0
            return [x,gval,ier,its]

    print('max iterations exceeded')    
    ier = 1        
    return [x,g1,ier,its]

def evalg(x):
# The least squares function you are minimizing

    F = evalF(x)
    g = F[0]**2 + F[1]**2 + F[2]**2
    return g

def eval_gradg(x):
# Grad g
    F = evalF(x)
    J = evalJ(x)
    
    gradg = np.transpose(J).dot(F)
    return gradg

def evalF(x): 
# vector function that you want to find the roots of
    F = np.zeros(3)
    F[0] = x[0] + np.cos(x[0]*x[1]*x[2]) - 1
    F[1] = (1-x[0])**0.25 + x[1] + 0.05*x[2]**2 - 0.15*x[2] - 1
    F[2] = -x[0]**2 - 0.1*x[1]**2 + 0.01*x[1] + x[2] - 1
    return F

def evalJ(x): 
# Jacobian of the vector function you want to find the roots of
    J = np.array([[1,-x[0]*x[2]*np.sin(x[0]*x[1]*x[2]),-x[0]*x[1]*np.sin(x[0]*x[1]*x[2])], 
        [(-1/4)*(1-x[0])**(-3/4),1,0.1*x[2]-0.15],
        [-2*x[0],0.2*x[1]+0.01,1]])
    return J

driver()