## APPM4600 Homework 7
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
    plt.legend()
    plt.show()


driver()