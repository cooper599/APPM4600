## APPM 4600 Lab 9
# Cooper Wark

import numpy as np

def driver():
    phi = eval_legendre(2,0.5)
    for i in range(2+1):
        print(phi[i])

def eval_legendre(n,x):
    phi = np.zeros(n+1)
    phi[0] = 1.0
    phi[1] = x
    for i in range (2,n+1):
        k = i-1 # i starts at 2, k starts at 1. Sets k to n since i to n+1, k = n+1-1 = n
        phi[i] = 1.0/(k+1) * ((2*k+1)*x*phi[k]-k*phi[k-1])
    return phi

driver()