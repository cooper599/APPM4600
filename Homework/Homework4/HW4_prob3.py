## APPM 4600 Homework 4, Problem 3
# Cooper Wark

import numpy as np

def driver():
    f = lambda x: np.e**(3*x) - 27*(x**6) + 27*(x**4)*np.e**x - 9*(x**2)*np.e**(2*x)
    fp = lambda x: 3*np.e**(3*x) - 162*(x**5)+108*(x**3)*np.e**x + 27*(x**4)*np.e**x - 18*x*np.e**(2*x) - 18*(x**2)*np.e**(2*x)
    fpp = lambda x: 9*np.e**(3*x)-810*x**4+324*(x**2)*np.e**x+108*(x**3)*np.e**x+108*(x**3)*np.e**x+27*(x**4)*np.e**x-18*np.e**(2*x)-36*x*np.e**(2*x)-36*x*np.e**(2*x)-36*(x**2)*np.e**(2*x)

    p0 = 4
    Nmax = 100
    tol = 1.e-12

    alpha = 3.7330793798506652 # root approx to check order of convergence

    print('--------------------')
    (p,pstar,info,it) = newton(f,fp,p0,tol, Nmax)
    print('the approximate root is', '%16.16e' % pstar)
    print('the error message reads:', '%d' % info)
    print('Number of iterations:', '%d' % it)
    k = 15
    e_curr = abs(p[k+1] - alpha)
    e_prev = abs(p[k] - alpha)
    estimated_p = np.log(e_curr) / np.log(e_prev)
    print(f'Estimated order of convergence p: {estimated_p:.2f}')
    print('--------------------')

    m = 3
    (p,pstar,info,it) = newtonMult(f,fp,p0,m,tol, Nmax)
    print('the approximate root is', '%16.16e' % pstar)
    print('the error message reads:', '%d' % info)
    print('Number of iterations:', '%d' % it)
    k = 1
    e_curr = abs(p[k+1] - alpha)
    e_prev = abs(p[k] - alpha)
    estimated_p = np.log(e_curr) / np.log(e_prev)
    print(f'Estimated order of convergence p: {estimated_p:.2f}')
    print('--------------------')

    (p,pstar,info,it) = newton_gofx(f,fp,fpp,p0,tol, Nmax)
    print('the approximate root is', '%16.16e' % pstar)
    print('the error message reads:', '%d' % info)
    print('Number of iterations:', '%d' % it)
    k = 2
    e_curr = abs(p[k+1] - alpha)
    e_prev = abs(p[k] - alpha)
    estimated_p = np.log(e_curr) / np.log(e_prev)
    print(f'Estimated order of convergence p: {estimated_p:.2f}')
    print('--------------------')
    

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

# Literally just added m in front f(x)/f'(x) 
def newtonMult(f,fp,p0,m,tol,Nmax):
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
      p1 = p0-m*f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

def newton_gofx(f,fp,fpp,p0,tol,Nmax):
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
      p1 = p0-(f(p0)*fp(p0))/((fp(p0)**2)-f(p0)*fpp(p0))
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