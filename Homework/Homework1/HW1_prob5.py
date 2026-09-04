# APPM 4600 Homework 1 Problem 5
# Cooper Wark

import numpy as np
import matplotlib.pyplot as plt

cosEquation = lambda x, delta: np.cos(x+delta) - np.cos(x)
cosExpression = lambda x, delta: -2*np.sin(x+(delta/2))*np.sin(delta/2)

delta_arr = np.logspace(-16,0,17)

x_arr = [np.pi, 10**6]

diff_pi = cosEquation(x_arr[0],delta_arr)-cosExpression(x_arr[0],delta_arr)
diff_106 = cosEquation(x_arr[1],delta_arr)-cosExpression(x_arr[1],delta_arr)

print("Difference when x = pi:", diff_pi)
print("Difference when x = 10^6:", diff_106)
                                                        
plt.loglog(abs(delta_arr),abs(diff_pi),label="x=pi")
plt.loglog(abs(delta_arr),abs(diff_106),label="x=10^6")
plt.xlabel("deltas")
plt.ylabel("Difference between equation and manipulated expression")
plt.title("Difference between equation and manupulated expression for range of deltas")
plt.legend()
plt.show()

# Part c, using same x_arr and delta_arr
cos_Taylor = lambda x, delta: -delta*np.sin(x)-(delta**2)/2 * np.cos(x)

diff_myalg_pi = cos_Taylor(x_arr[0],delta_arr)-cosExpression(x_arr[0],delta_arr)
diff_myalg_106 = cos_Taylor(x_arr[1],delta_arr)-cosExpression(x_arr[1],delta_arr)

plt.loglog(abs(delta_arr),abs(diff_myalg_pi),label="x=pi")
plt.loglog(abs(delta_arr),abs(diff_myalg_106),label="x=10^6")
plt.xlabel("deltas")
plt.ylabel("Difference between equation and manipulated expression")
plt.title("Difference between Taylor Algorithm and manipulated expression for range of deltas")
plt.legend()
plt.show()

print("Taylor exp Difference when x = pi:", diff_myalg_pi)
print("Taylor exp Difference when x = 10^6:", diff_myalg_106)

print("Difference when x = 10^6:", diff_106)
print("Taylor exp Difference when x = 10^6:", diff_myalg_106)