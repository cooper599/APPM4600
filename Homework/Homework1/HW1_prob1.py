# APPm 4600 Homework 1 Problem 1
# Cooper Wark
import numpy as np
import matplotlib.pyplot as plt

pExp = lambda x: x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 + 5376*x**3-4608*x**2+2304*x-512
pPar = lambda x: (x-2)**9

# start, stop (exclusive), step
x = np.arange(1.920,2.081,0.001)

plt.plot(x,pExp(x),label="Coefficient Expansion")
plt.plot(x,pPar(x), label="Expression")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Plot of p(x) using coef expansion and expression")
plt.legend()
plt.show()

# Expression works better, the coefficient expansion leads to jagged ups and downs
# Floating point arithemetic rounding on the multiplication causing errors in combination with exponent
