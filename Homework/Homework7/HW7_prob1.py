## APPM4600 Homework 7, Prob 1
# Cooper Wark

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

def driver():
    N = 19
    interp(N)

def interp(N):
    h = 2.0 / (N-1)
    xi = -1.0 + np.arange(N) * h
    # j = np.arange(1, N + 1)
    # xi = np.cos((2 * j - 1) * np.pi / (2 * N))
    yi = 1.0 / (1.0 + (10 * xi)**2)

    V = np.vander(xi,N)
    c = inv(V).dot(yi)

    xfine = np.linspace(-1,1,1001)
    ffine = 1.0 / (1.0 + (10*xfine)**2)
    pfine = np.zeros(len(xfine))
    # Make polynomail line byline
    for i in range(N):
        pfine += c[i] * (xfine**(N-1-i))

    plt.plot(xi,yi,'o',label="data points")
    plt.plot(xfine,ffine,label="f(x)")
    plt.plot(xfine,pfine,label="p(x) approximation")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("f(x) vs monomial p(x) approx using Vandermonde")
    plt.legend()
    plt.show()


driver()