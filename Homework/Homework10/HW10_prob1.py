## APPM 4600 Homeowrk 10 Prob 1
# Cooepr Wark

import numpy as np
import matplotlib.pyplot as plt

def driver():
    a = 0
    b = np.pi
    I_ex = -12.070346316389632 # Exact value of integral
    tol = 1e-10
    # initial n start at lowest possible
    n_trap = 2
    n_simp = 2
    # stuff for plotting
    h_trap = []
    h_simp = []
    err_trap = []
    err_simp = []

    I_est_trap = 0
    rerr_trap = calcRelErr(I_est_trap,I_ex)
    while rerr_trap > tol:
        # Update n at beginning and nit so check at end doesn't update n before exiting
        n_trap = 2*n_trap
        I_est_trap = CompTrap(a,b,n_trap,f)
        rerr_trap = calcRelErr(I_est_trap,I_ex)
        h_trap.append((b-a)/n_trap)
        err_trap.append(np.abs(I_ex-I_est_trap))

    rerr_trap = calcRelErr(I_est_trap,I_ex)
    print("I estimate trapezoidal:", I_est_trap)
    print("Relative Error: ", rerr_trap)
    print("n trap (how many factors of 2)", n_trap)
    print()

    I_est_simp = 0
    rerr_simp = calcRelErr(I_est_simp,I_ex)
    while rerr_simp > tol:
        # Update n at beginning and nit so check at end doesn't update n before exiting
        n_simp = 2*n_simp
        I_est_simp = CompSimp(a,b,n_simp,f)
        rerr_simp = calcRelErr(I_est_simp,I_ex)
        h_simp.append((b-a)/n_simp)
        err_simp.append(np.abs(I_ex-I_est_simp))

    rerr_simp = calcRelErr(I_est_simp,I_ex)
    print("I estimate simpson:", I_est_simp)
    print("Relative Error: ", rerr_simp)
    print("n simp (how many factors of 2)", n_simp)

    plt.loglog(h_trap,err_trap,label="Trapezoidal Error")
    plt.loglog(h_simp,err_simp,label="Simpson's Error")
    plt.xlabel("h")
    plt.ylabel("Absolute Error")
    plt.legend()
    plt.show()

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

def f(x):
    f = np.e**x * np.cos(x)
    return f

driver()