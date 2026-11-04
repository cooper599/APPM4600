## APPM4600 Homework 9 Prob 4
# Cooper Wark

import numpy as np
import matplotlib.pyplot as plt

def driver():  
    N = 1000
    x = np.linspace(0,5,N)

    # actual
    fx = np.sin(x)
    #maclauren
    f_max = x - x**3.0/6.0 + x**5.0/120.0
    # M= N = 3
    f_M3_N3 = (x-7.0/60.0 * x**3.0)/(1.0+1.0/20.0 * x**2.0)
    # M = 4, N = 2
    f_M4_N2 = (x-7.0/60.0 * x**3.0)/(1.0+1.0/20.0*x**2.0)
    # M = 2, N = 4
    f_M2_N4 = x / (1.0+1.0/6.0*x**2.0+7.0/360.0*x**4.0)

    # Erros
    err_mac = np.abs(f_max-fx)
    err_M3_N3 = np.abs(f_M3_N3-fx)
    err_M4_N2 = err_M3_N3
    err_M2_N4 = np.abs(f_M2_N4-fx)

    # Plotting
    plt.plot(x,err_mac,label="Maclauren Error")
    plt.plot(x,err_M3_N3,label = "M3 N3 Approx") # same as below
    plt.plot(x,err_M4_N2,label="M4 N2 Approx") # same as above
    plt.plot(x,err_M2_N4,label="M2 N4 Approx")
    plt.xlabel("x")
    plt.ylabel("error")
    plt.title("6th degree Padé and Maclauren error for sin(x) approximation")
    plt.legend()
    plt.show()

driver()