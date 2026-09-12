## APPM 4600 Homework 2, Prob 3
# Cooper Wark
import numpy as np
import math

x = 9.999999995000000*(10**-10)

y = math.e**x

alg_y = y-1

print(alg_y)

print(x**2/6)

P = lambda x: x + x**2/2
print(f"Taylor approximation of e^x - 1: {P(x):.16e}")
rel_err = abs(P(x)-alg_y)/P(x)
print("Relative function errors b/w alg and Taylor approx: ", rel_err)