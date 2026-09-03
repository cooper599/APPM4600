# APPM 4600 Homework 1 Problem 5
# Cooper Wark

import numpy as np
import matplotlib.pyplot as plt

cosEquation = lambda x, delta: np.cos(x+delta) - np.cos(x)
cosExpression = lambda x, delta: -2*np.sin(x+(delta/2))*np.sin(delta/2)

delta_arr = np.logspace(-16,0,17)

x_arr = [np.pi, 10**6]

plt.semilogx(delta_arr,cosEquation(x_arr[0],delta_arr)-cosExpression(x_arr[0],delta_arr),label="x=pi")
plt.semilogx(delta_arr,cosEquation(x_arr[1],delta_arr)-cosExpression(x_arr[1],delta_arr),label="x=10^6")
plt.xlabel("deltas")
plt.ylabel("Difference between equation and manipulated expression")
plt.title("Difference between equation and manupulated expression for range of deltas")
plt.legend()
plt.show()