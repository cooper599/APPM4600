## APPM 4600 Homework 4, Problem 4
# Cooper Wark

import numpy as np
import matplotlib.pyplot as plt

def driver():
    f = lambda x: x**6 - x - 1
    fp = lambda x: 6*x**5 - 1

    tol = 1.e-14
    Nmax = 100

    p0 = 2
    (p_newt,pstar,info,it) = newton(f,fp,p0,tol,Nmax)
    print("root approx: ", pstar)
    print("number it: ", it)
    alpha = pstar # use newton's approx of root for actual root

    p1 = 1
    (p_sec,pstar,info,it) = secant(f,p0,p1,tol,Nmax)
    print("root approx: ", pstar)
    print("number it: ", it)

    # Error table calculations
    [newt_ek,newt_ekp1] = calcError(p_newt,alpha)
    [sec_ek,sec_ekp1] = calcError(p_sec,alpha)

    plt.loglog(newt_ek,newt_ekp1,label="Newton's Method")
    plt.loglog(sec_ek,sec_ekp1,label="Secant Method")
    plt.xlabel("log(e_k)")
    plt.ylabel("log(e_k+1)")
    plt.title("log(e_k+1) vs log(e_k) for Newton's Method and Secant Method")
    plt.legend()
    plt.show()

def calcError(p,alpha):
    # Calculates e_k and e_k+1 array
    # p - root approximations, alpha - closest numerical approx or actual root found
    all_err = abs(p-alpha)
    val_err = all_err > 1.e-16 # floor to prevent log(0)
    err = all_err[val_err]
    ek = err[:-1] # first to second last
    ekp1 = err[1:] # second to last 
    return [ek, ekp1]


def secant(f,p0,p1,tol,Nmax):
  """
  Secant iteration.
  
  Inputs:
    f - function
    p0,p1   - initial guesses for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
  """
  p = np.zeros(Nmax+2);
  p[0] = p0
  p[1] = p1
  f0 = f(p0)
  f1 = f(p1)

  for it in range(1,Nmax+1):
      if abs(f(p0)-f(p1)) == 0:
          pstar = p0
          info = 1
          it = 0
          return [p[:it+2],pstar,info,it] # p[:it+2] returns only non zero iterations in p array
      p2 = p1-f1*(p1-p0)/(f1-f0)
      p[it+1] = p2
      if (abs(p2-p1) < tol):
          pstar = p2
          info = 0
          return [p[:it+2],pstar,info,it]
      p0 = p1
      f0 = f1
      p1 = p2
      f1 = f(p2)
  pstar = p2
  info = 1
  return [p[:it+2],pstar,info,it]

def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
  """
  p = np.zeros(Nmax+1);
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p[:it+2],pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p[:it+2],pstar,info,it]

driver()