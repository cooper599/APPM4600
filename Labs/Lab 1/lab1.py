# Lab 1 code 3.1.1
x = [1,2,3]
print(x)
print(3*x) # Outputs 1,2,3,1,2,3,1,2,3

import numpy as np

y = np.array([1,2,3])
print(y)
print(3*y) # Outputs 3 6 9

## Section 3.12
print('this is 3y', 3*y)

## Section 3.1.3
import matplotlib.pyplot as plt

X = np.linspace(0,2*np.pi,100)
Ya = np.sin(X)
Yb = np.cos(X)

plt.plot(X,Ya)
plt.plot(X,Yb)
plt.xlabel("x")
plt.ylabel('y')
plt.show()
print(len(X)) # size 100

## Section 3.2, 0 to 10 1 step size
x = np.linspace(0,10,11)
y = np.arange(0,11,1) 

print("The first 3 entries of x are",x[0:3])

## Section 3.3, plotting
w = 10**(-np.linspace(1,10,10))
x = np.arange(1,len(w)+1,1)
print(w)
print(x)

plt.semilogy(x,w,label="w")
plt.xlabel('x')
plt.ylabel('y')

s = 3*w
plt.semilogy(x,s,label="s")
plt.legend()
plt.show()