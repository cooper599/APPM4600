## APPM 4600 Homework 3 prob 2
# Cooper Wark

import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt

def driver():
    Ti = 20
    T0 = -15
    alpha = 0.138*(10**-6) # m^2/s
    tol = 1e-13 # 10^-13
    t = 60*24*60*60 # day -> hour -> minute -> seconds 

    f = lambda x: erf(x/(2*np.sqrt(alpha*t))) - 3/7
    fp = lambda x: (1/np.sqrt(np.pi*alpha*t)) * np.e**(-(x/(2*np.sqrt(alpha*t)))**2)

    # xbar_test = [-1,0,1,2]
    # print(f(xbar_test)) f(1) > 0
    plot = False

    xbar = np.linspace(0,1,100)

    if plot:
        plt.plot(xbar,f(xbar))
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("f(x) for [0 to 1]")
        plt.axhline(0,linestyle="dashed", color = "black")
        plt.show()

    # part b
    a0 = 0
    b0 = 1
    [astar, ier] = bisection(f,a0,b0,tol)
    print("Biscetion Root Approx: ", astar)

    # part c
    Nmax = 100
    x0 = 0.01 # m
    [p,pstar,info,it] = newton(f,fp,x0,tol,Nmax)
    print("Newton Root Approx x = 0.01m : ", pstar)
    x0 = 1.0 # m
    [p1,pstar1,info1,it1] = newton(f,fp,x0,tol,Nmax)
    print("Newton Root Approx x = 1 m : ", pstar1)



# % Newton
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

# Bisection
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

    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
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
    print('count = ', count)
    return [astar, ier]

driver()