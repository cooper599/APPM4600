## APPM4600 Homework 7, Prob 2
# Cooper Wark

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

def driver():
    N = [19,45,99]
    for i in range(len(N)):
        interp_baycentric(N[i])

def interp_baycentric(N):
    h = 2.0 / (N-1)
    xi = -1.0 + np.arange(N) * h
    yi = 1.0 / (1.0 + (10 * xi)**2)

    w = np.ones(N)
    for i in range(N):
        diff = xi[i]-xi
        diff[i] = 1.0
        w[i] = 1.0/np.prod(diff)

    xfine = np.linspace(-1,1,1001)
    ffine = 1.0 / (1.0 + (10*xfine)**2)

    pfine = np.zeros(len(xfine))
    for j in range(len(xfine)):
        x = xfine[j]
        matches = np.abs(x-xi) < 1e-14 # if lands really clsoe to node
        if np.any(matches):
            pfine[j] = yi[matches][0]
        else:
            terms = w/(x-xi)
            pfine[j] = np.sum(terms * yi) / np.sum(terms)

    plt.plot(xi,yi,'o',label="data points")
    plt.plot(xfine,ffine,label="f(x)")
    plt.plot(xfine,pfine,label="p(x) baycentric approx")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.ylim(-1,2)
    plt.title(f"f(x) vs p(x) using baycentric approx for N = {N}")
    plt.legend()
    plt.show()

driver()