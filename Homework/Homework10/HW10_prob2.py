## APPM 4600 Homeowrk 10 Prob 2
# Cooepr Wark

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def driver():
    a = -5.0
    b = 5.0
    # I_ex = -12.070346316389632 # Exact value of integral
    tol = 1e-6
    f = lambda x: 1/(1+x**2)

    ntrap = 1291
    nsimp = 108

    Iquad6, _, info6 = quad(f, a, b, epsabs=1e-6, epsrel=1e-6, full_output=True)
    Iquad4, _, info4 = quad(f, a, b, epsabs=1e-4, epsrel=1e-4, full_output=True)  
    neval6 = info6['neval']
    neval4 = info4['neval']

    # Found from exact integral
    Iex = 2.746801533890032

    Itrap = CompTrap(a,b,ntrap,f)
    Isimp = CompSimp(a,b,nsimp,f)

    print("Trapezoid Estimate: ", Itrap, "Absolute Error: ", calcAbsErr(Itrap,Iex), "Neval: ", ntrap)
    print("Simpon's Estimate: ", Isimp, "Absolute Error: ", calcAbsErr(Isimp,Iex), "Neval: ", nsimp)
    print("Quad (1e-6) Estimate: ", Iquad6, "Absolute Error: ", calcAbsErr(Iquad6, Iex), "Neval: ", neval6)
    print("Quad (1e-4) Estimate: ", Iquad4, "Absolute Error: ", calcAbsErr(Iquad4, Iex), "Neval: ", neval4)


    
def CompTrap(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    I_trap = h*f(xnode[0])*1/2
    for j in range(1,n):
        I_trap = I_trap+h*f(xnode[j])
    I_trap= I_trap + 1/2*h*f(xnode[n])
    return I_trap

def CompSimp(a,b,n,f):
    h = (b-a)/n
    sum_even = 0    
    sum_odd = 0    
    xnode = a+np.arange(0,n+1)*h
    for j in range(1,int(n/2) -1 + 1): # +1 for indexing bc it stops 1 early
        sum_even += f(xnode[2*j])
    for j2 in range(0,int(n/2)): # convert to int
        sum_odd += f(xnode[2*j2+1])
    I_simp = (h/3)*(f(xnode[0]) + 2*sum_even + 4*sum_odd + f(xnode[n]))
    return I_simp

def calcRelErr(I_est,I_ex):
    rerr = np.abs(I_est-I_ex)/np.abs(I_ex)
    return rerr

def calcAbsErr(I_est,I_ex):
    aerr = np.abs(I_est-I_ex)
    return aerr

driver()