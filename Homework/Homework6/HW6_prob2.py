## APPM4600 Homework 6 Prob 2
# Cooper Wark

import numpy as np
from numpy.linalg import norm
from numpy.linalg import inv

def driver():
    tol = 1e-12
    Nmax = 100

    x01 = [1.0, 1.0]
    x02 = [1.0, -1.0]
    x03 = [0.0, 0.0]

    print("Newton: ")
    [xstar1,ier1,its1] = Newton(x01,tol,Nmax)
    [xstar2,ier2,its2] = Newton(x02,tol,Nmax)
    # [xstar3,ier3,its3] = Newton(x03,tol,Nmax)
    printRes(xstar1,ier1,its1)
    printRes(xstar2,ier2,its2)
    # printRes(xstar3,ier3,its3)

    print("Lazy Newton: ")
    [xstar1,ier1,its1] = LazyNewton(x01,tol,Nmax)
    [xstar2,ier2,its2] = LazyNewton(x02,tol,Nmax)
    # [xstar3,ier3,its3] = LazyNewton(x03,tol,Nmax)
    printRes(xstar1,ier1,its1)
    printRes(xstar2,ier2,its2)
    # printRes(xstar3,ier3,its3)

    print("Broyden: ")
    [xstar1,ier1,its1] = Broyden(x01,tol,Nmax)
    [xstar2,ier2,its2] = Broyden(x02,tol,Nmax)
    # [xstar3,ier3,its3] = Broyden(x03,tol,Nmax)
    printRes(xstar1,ier1,its1)
    printRes(xstar2,ier2,its2)
    # printRes(xstar3,ier3,its3)    

def Broyden(x0,tol,Nmax):
    '''tol = desired accuracy
    Nmax = max number of iterations'''

    '''Sherman-Morrison 
   (A+xy^T)^{-1} = A^{-1}-1/p*(A^{-1}xy^TA^{-1})
    where p = 1+y^TA^{-1}Ax'''

    '''In Newton
    x_k+1 = xk -(G(x_k))^{-1}*F(x_k)'''


    '''In Broyden 
    x = [F(xk)-F(xk-1)-\hat{G}_k-1(xk-xk-1)
    y = x_k-x_k-1/||x_k-x_k-1||^2'''

    ''' implemented as in equation (10.16) on page 650 of text'''
    
    '''initialize with 1 newton step'''
    
    A0 = evalJ(x0)

    v = evalF(x0)
    A = np.linalg.inv(A0)

    s = -A.dot(v)
    xk = x0+s
    for  its in range(Nmax):
       '''(save v from previous step)'''
       w = v
       ''' create new v'''
       v = evalF(xk)
       '''y_k = F(xk)-F(xk-1)'''
       y = v-w;                   
       '''-A_{k-1}^{-1}y_k'''
       z = -A.dot(y)
       ''' p = s_k^tA_{k-1}^{-1}y_k'''
       p = -np.dot(s,z)                 
       u = np.dot(s,A) 
       ''' A = A_k^{-1} via Morrison formula'''
       tmp = s+z
       tmp2 = np.outer(tmp,u)
       A = A+1./p*tmp2
       ''' -A_k^{-1}F(x_k)'''
       s = -A.dot(v)
       xk = xk+s
       if (norm(s)<tol):
          alpha = xk
          ier = 0
          return[alpha,ier,its]
    alpha = xk
    ier = 1
    return[alpha,ier,its]

def printRes(xstar,ier,its):
    print("x approx: ", xstar[0])
    print("y approx: ", xstar[1])
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

def evalF(x): 
# vector function that you want to find the roots of
    F = np.zeros(2)
    F[0] = x[0]**2 + x[1]**2 - 4
    F[1] = np.e**x[0] + x[1] - 1
    return F

def evalJ(x): 
# Jacobian of the vector function you want to find the roots of
    J = np.array([[2*x[0],2*x[1]], 
        [np.e**x[0],1]])
    return J

driver()