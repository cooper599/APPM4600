## APPM 4600 Homework 4, Problem 3
# Cooper Wark

import numpy as np

def driver():
    f = lambda x: np.e**(3*x) - 27*(x**6) + 27*(x**4)*np.e**x - 9*(x**2)*np.e**(2*x)
    fp = lambda x: 3*np.e**(3*x) - 162*(x**5)+108*(x**3)*np.e**x + 27*(x**4)*np.e**x - 18*x*np.e**(2*x) - 18*(x**2)*np.e**(2*x)
    
    p0 = 4
    Nmax = 100
    tol = 1.e-12

    (p,pstar,info,it) = newton(f,fp,p0,tol, Nmax)
    print('the approximate root is', '%16.16e' % pstar)
    print('the error message reads:', '%d' % info)
    print('Number of iterations:', '%d' % it)

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
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

driver()