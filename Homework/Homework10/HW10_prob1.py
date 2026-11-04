## APPM 4600 Homeowrk 10 Prob 1
# Cooepr Wark

import numpy as np
import matplotlib.pyplot as plt

def driver():
    a = 0
    b = np.pi
    I_ex = -12.070346316389632 # Exact value of integral
    tol = 1e-10
    n_trap = 2


    I_est_trap = 0
    nit = 0
    rerr_trap = calcRelErr(I_est_trap,I_ex)
    while rerr_trap > tol:
        # Update n at beginning and nit so check at end doesn't update n before exiting
        n_trap = 2*n_trap
        nit += 1
        I_est_trap = CompTrap(a,b,n_trap,f)
        rerr_trap = calcRelErr(I_est_trap,I_ex)

    rerr_trap = calcRelErr(I_est_trap,I_ex)
    print("I estimate trapezoidal:", I_est_trap)
    print("Relative Error: ", rerr_trap)
    print("n trap (how many factors of 2)", n_trap)

def CompTrap(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h
    I_trap = h*f(xnode[0])*1/2
    for j in range(1,n):
        I_trap = I_trap+h*f(xnode[j])
    I_trap= I_trap + 1/2*h*f(xnode[n])
    return I_trap

def calcRelErr(I_est,I_ex):
    rerr = np.abs(I_est-I_ex)/np.abs(I_ex)
    return rerr

def f(x):
    f = np.e**x * np.cos(x)
    return f

driver()