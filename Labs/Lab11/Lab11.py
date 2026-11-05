## APPM 4600 Lab 11
# Cooepr Wark

import numpy as np
import matplotlib.pyplot as plt

def driver():
    method = [eval_composite_trap,eval_composite_simpsons,eval_gauss_quad]

    for i in range(len(method)):
        method_comparison(method[i])

def method_comparison(method):
    a = 0.1
    b = 2.0
    f = lambda x: np.sin(1/x); labl = "sin(1/x)"
    # absolute tolerance for adaptive quad 
    tol = 1e-3
    # machine eps in numpy
    eps = np.finfo(float).eps
    # number of nodes and weights, per subinterval 
    Ms = np.arange(2,15); nM = len(Ms)
    # storage for error
    err_old = np.zeros((nM,))
    err_new = np.zeros((nM,))
    I_true = 1.145580834099501
    for iM in range(nM):
        M = Ms[iM]; 
        # non adaptive routine 
        # Note: the _,_ are dummy vars/Python convention 
        # to store uneeded returns from the routine
        I_old,_,_ = method(M,a,b,f)
        # adaptive routine
        I_new,X,nsplit = adaptive_quad(a,b,f,tol,M,method)
        err_old[iM] = np.abs(I_old-I_true)/I_true
        err_new[iM] = np.abs(I_new-I_true)/I_true 
        # clean the error for nice plots
        if err_old[iM] < eps:
            err_old[iM] = eps 
        if err_new[iM] < eps:
            err_new[iM] = eps
        # save grids for M = 2
        if M == 2:
            mesh = X

    # plot the old and new error for each f and M
    fig,ax = plt.subplots(1,2)
    ax[0].semilogy(Ms,err_old,'ro--')
    ax[0].set_ylim([1e-16,2]);
    ax[0].set_xlabel('$M$')
    ax[0].set_title('Non-adaptive')
    ax[0].set_ylabel('Relative Error');
    ax[1].semilogy(Ms,err_new,'ro--',label=labl)
    ax[1].set_ylim([1e-16,2]);
    ax[1].set_xlabel('$M$')
    ax[1].set_title('Adaptive')
    ax[1].legend()
    plt.show()

    # plot the adaptive mesh for M=2
    fig,ax = plt.subplots(1,1)
    ax.semilogy(mesh,f(mesh),'ro',label=labl)
    ax.legend()
    plt.show()

    print("Nsplit: ",nsplit)

def CompGauss(a,b,n,f):
    x,w = lgwt(n,a,b)
    I_gauss = np.sum(f(x)*w)
    return I_gauss

def lgwt(N,a,b):
    """
     Computes the Legendre-Gauss nodes and weights on an interval
     [a,b] with truncation order N.

     Suppose you have a continuous function f(x) which is defined on [a,b]
     which you can evaluate at any x in [a,b]. Simply evaluate it at all of
     the values contained in the x vector to obtain a vector f. Then compute
     the definite integral using np.sum(f*w)

     Written by Greg von Winckel - 02/25/2004
     translated to Python - 10/30/2022
    """
    N = N-1
    N1 = N+1
    N2 = N+2
    eps = np.finfo(float).eps
    xu = np.linspace(-1,1,N1)
    # Initial guess
    y = np.cos((2*np.arange(0,N1)+1)*np.pi/(2*N+2))+(0.27/N1)*np.sin(np.pi*xu*N/N2)
    # Legendre-Gauss Vandermonde Matrix
    L = np.zeros((N1,N2))
    # Compute the zeros of the N+1 Legendre Polynomial
    # using the recursion relation and the Newton-Raphson method
    y0 = 2.
    one = np.ones((N1,))
    # Iterate until new points are uniformly within epsilon of old points
    while np.max(np.abs(y-y0)) > eps:
        L[:,0] = one
        L[:,1] = y
        for k in range(2,N1+1):
            L[:,k] = ((2*k-1)*y*L[:,k-1]-(k-1)*L[:,k-2])/k
        lp = N2*(L[:,N1-1]-y*L[:,N2-1])/(1-y**2)
        y0 = y
        y = y0-L[:,N2-1]/lp
    # Linear map from [-1,1] to [a,b]
    x = (a*(1-y)+b*(1+y))/2
    # Compute the weights
    w = (b-a)/((1-y**2)*lp**2)*(N2/N1)**2
    return x,w

def eval_composite_trap(n,a,b,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    I_trap = h*f(xnode[0])*1/2
    for j in range(1,n):
        I_trap = I_trap+h*f(xnode[j])
    I_trap= I_trap + 1/2*h*f(xnode[n])
    x = xnode
    w = None
    return I_trap,x,w

def eval_composite_simpsons(n,a,b,f):
    h = (b-a)/n
    sum_even = 0    
    sum_odd = 0    
    xnode = a+np.arange(0,n+1)*h
    for j in range(1,int(n/2) -1 + 1): # +1 for indexing bc it stops 1 early
        sum_even += f(xnode[2*j])
    for j2 in range(0,int(n/2)): # convert to int
        sum_odd += f(xnode[2*j2+1])
    I_simp = (h/3)*(f(xnode[0]) + 2*sum_even + 4*sum_odd + f(xnode[n]))
    x = xnode
    w = None
    return I_simp,x,w

def eval_gauss_quad(M,a,b,f):
  """
  Non-adaptive numerical integrator for int_a^b f(x)w(x)dx
  Input:
    M - number of quadrature nodes
    a,b - interval [a,b]
    f - function to integrate
  Output:
    I_hat - approx integral
    x - quadrature nodes
    w - quadrature weights
  Currently uses Gauss-Legendre rule
  """
  x,w = lgwt(M,a,b)
  I_gauss = np.sum(f(x)*w)
  return I_gauss,x,w

def adaptive_quad(a,b,f,tol,M,method):
  """
  Adaptive numerical integrator for \int_a^b f(x)dx
  
  Input:
  a,b - interval [a,b]
  f - function to integrate
  tol - absolute accuracy goal
  M - number of quadrature nodes per bisected interval
  method - function handle for integrating on subinterval
         - eg) eval_gauss_quad, eval_composite_simpsons etc.
  
  Output: I - the approximate integral
          X - final adapted grid nodes
          nsplit - number of interval splits
  """
  # 1/2^50 ~ 1e-15
  maxit = 50
  left_p = np.zeros((maxit,))
  right_p = np.zeros((maxit,))
  s = np.zeros((maxit,1))
  left_p[0] = a; right_p[0] = b;
  # initial approx and grid
  s[0],x,_ = method(M,a,b,f);
  # save grid
  X = []
  X.append(x)
  j = 1;
  I = 0;
  nsplit = 1;
  while j < maxit:
    # get midpoint to split interval into left and right
    c = 0.5*(left_p[j-1]+right_p[j-1]);
    # compute integral on left and right spilt intervals
    s1,x,_ = method(M,left_p[j-1],c,f); X.append(x)
    s2,x,_ = method(M,c,right_p[j-1],f); X.append(x)
    if np.max(np.abs(s1+s2-s[j-1])) > tol:
      left_p[j] = left_p[j-1]
      right_p[j] = 0.5*(left_p[j-1]+right_p[j-1])
      s[j] = s1
      left_p[j-1] = 0.5*(left_p[j-1]+right_p[j-1])
      s[j-1] = s2
      j = j+1
      nsplit = nsplit+1
    else:
      I = I+s1+s2
      j = j-1
      if j == 0:
        j = maxit
  return I,np.unique(X),nsplit

driver()