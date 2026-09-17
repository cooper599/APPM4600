## APPM 4600 Lab 3, Newtons + Other method
# Core function Newtons + bisection from class code
# Cooper Wark
import numpy as np

def driver():
# use routines    
    f = lambda x: np.e**(x**2+7*x-30) - 1
    fp = lambda x: np.e**(x**2+7*x-30) * (2*x+7)
    fpp = lambda x: np.e**(x**2+7*x-30) * (2*x+7)**2 + 2*np.e**(x**2+7*x-30)
    a = 2
    b = 4.5
    x0 = 4.5
    tol = 1e-10
    Nmax = 100
    [astar, ier, count] = bisection(f,a,b,tol)
    [p1,pstar1,ier1,it1] = newton(f,fp,x0,tol,Nmax)
    [p2,pstar2,ier2,it2] = BiNewton(f,fp,fpp,a,b,tol,Nmax)
    print("     Bisection Method       ")
    print("Root approx: ", astar)
    print("Iterations: ", count)

    print("     Newtons Method       ")
    print("Root approx: ", pstar1)
    print("Iterations: ", it1)

    print("     Hybrid Method       ")
    print("Root approx: ", pstar2)
    print("Iterations: ", it2)    

# Newtons Method
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

# regular bisection
def bisection(f,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 
    count = 0
    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier, count]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier, count]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier, count]

    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier, count]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    return [astar, ier, count]

# Combined function
def BiNewton(f,fp,fpp,a,b,tol,Nmax):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 

    p = np.zeros(Nmax+1)
    fa = f(a)
    fb = f(b)
    count = 0
    if (fa*fb>0):
       ier = 1
       pstar = a
       return [p, pstar, ier, count]

#   verify end points are not a root 
    if (fa == 0):
      pstar = a
      ier =0
      return [p, pstar, ier, count]

    if (fb ==0):
      pstar = b
      ier = 0
      return [p, pstar, ier, count]

    d = 0.5*(a+b)
    while (abs(d-a) > tol):
      fd = f(d)
      if (fd ==0):
        pstar = d
        ier = 0
        return [p, pstar, ier, count]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)

      # Added this exit scenario for when d = astar is within basin of convergence, doesn't need to keep looping to within tolerance
      if abs(1-(fp(d)**2-f(d)*fpp(d))/fp(d)**2) < 1:
         ier = 0
         pstar = d
        #  return [astar, ier] 
         break # exits out of while
      count = count + 1
      
    # astar is p0
    p[0] = pstar
    for it in range(Nmax):
        p1 = pstar-f(pstar)/fp(pstar)
        p[it+1] = p1
        if (abs(p1-pstar) < tol):
            pstar = p1
            ier = 0
            return [p,pstar,ier,it + count]
        pstar = p1
    pstar = p1
    ier = 1
    return [p,pstar,ier,it + count]

driver()
