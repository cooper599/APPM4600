## APPM 4600 Lab 9
# Cooper Wark

import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import math
from scipy.integrate import quad

def driver():
    #  orign fnct
    # f = lambda x: math.exp(x)
    '''
    error close to 10^-1 at end points at -0.75, 0, 0.75 drops down to 10^-4 making parabolas between
    '''

    # function num 2
    f = lambda x: 1/(1. + x**2)
    '''
    error close to 10^-1 at end point again, drops down to 10^-3 at 4 locations instead of 3
    '''

# Interval of interest    
    a = -1
    b = 1
# weight function    
    w = lambda x: 1.
    # 3.3 excersise weight fnct
    # w = lambda x: 1/np.sqrt(1.-x**2)

# order of approximation
    n = 2

#  Number of points you want to sample in [a,b]
    N = 1000
    xeval = np.linspace(a,b,N+1)
    pval = np.zeros(N+1)

    for kk in range(N+1):
      pval[kk] = eval_legendre_expansion(f,a,b,w,n,xeval[kk])
      
    ''' create vector with exact values'''
    fex = np.zeros(N+1)
    for kk in range(N+1):
        fex[kk] = f(xeval[kk])
        
    plt.figure()    
    plt.plot(xeval,fex,'ro-', label= 'f(x)')
    plt.plot(xeval,pval,'bs--',label= 'Expansion') 
    plt.legend()
    plt.show()    
    
    err = abs(pval-fex)
    plt.semilogy(xeval,err,'ro--',label='error')
    plt.legend()
    plt.show()

def eval_legendre_expansion(f,a,b,w,n,x): 

#   This subroutine evaluates the Legendre expansion

#  Evaluate all the Legendre polynomials at x that are needed
# by calling your code from prelab 

# orig
  p = eval_legendre(n,x)
# chebshev
#   p = eval_chebyshev(n,x)

  # initialize the sum to 0 
  pval = 0.0    
  for j in range(0,n+1):
      # make a function handle for evaluating phi_j(x)
      phi_j = lambda x: eval_legendre(j,x)[j]
    #   phi_j = lambda x: eval_chebyshev(j,x)[j]

      # make a function handle for evaluating phi_j^2(x)*w(x)
      phi_j_sq = lambda x: phi_j(x)**2 * w(x)
      # use the quad function from scipy to evaluate normalizations
      norm_fac,err = quad(phi_j_sq,a,b)
      # make a function handle for phi_j(x)*f(x)*w(x)/norm_fac
      func_j = lambda x: phi_j(x)*f(x)*w(x)/norm_fac
      # use the quad function from scipy to evaluate coeffs
      aj,err = quad(func_j,a,b)
      # accumulate into pval
      pval = pval+aj*p[j] 
       
  return pval

def eval_legendre(n,x):
    phi = np.zeros(n+1)
    phi[0] = 1.0
    if n >= 1:
        phi[1] = x
    for i in range (2,n+1):
        k = i-1 # i starts at 2, k starts at 1. Sets k to n since i to n+1, k = n+1-1 = n
        phi[i] = 1.0/(k+1) * ((2*k+1.0)*x*phi[k]-k*phi[k-1])
    return phi

def eval_chebyshev(n,x):
    T = np.zeros(n+1)
    T[0] = 1.0
    if n >= 1:
        T[1] = x
    for i in range (2,n+1):
        k = i-1 # i starts at 2, k starts at 1. Sets k to n since i to n+1, k = n+1-1 = n
        T[i] = 2.0*x*T[k] - T[k-1]
    return T

driver()