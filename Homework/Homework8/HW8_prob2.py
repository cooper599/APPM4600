##APPM4600 Homework 8 Problem 2
# Cooper Wark

import numpy as np
import matplotlib.pyplot as plt

def driver():
    # Fncts
    f = lambda x: 1.0/(1.0+x**2)
    fp = lambda x: -2*x/(1.0+x**2)**2

    # Parms
    N = 5 # num nodes
    a = -5
    b = -5

    # Equispaceaced nodes
    xint = np.linspace(a,b,N+1)
    # Interpdata
    yint = np.zeros(N+1)
    ypint = np.zeros(N+1)
    for i in range(N+1):
       yint[i] = f(xint[i])
       ypint[i] = fp(xint[i])

    # Poitns for lagrange and hermite
    Neval = 1000
    xeval = np.linspace(a,b,Neval+1)
    # 4 interp methods
    yeval_Lagrange = np.zeros(Neval+1)
    yeval_Hermite = np.zeros(Neval+1)
    yeval_NatCubSpl = np.zeros(Neval+1)
    yeval_ClaCubSpl = np.zeros(Neval+1)
    fex = np.zeros(Neval + 1)
    for ii in range(Neval + 1):
        yeval_Lagrange[ii] = eval_lagrange(xeval[ii],xint,yint,N)
        yeval_Hermite[ii] = eval_hermite(xeval[ii],xint,yint,ypint,N)
    #     yeval_NatCubSpl[ii] = 
    #     yeval_ClaCubSpl[ii] = 
        fex[ii] = f(xeval[ii])


    plt.figure()
    plt.plot(xeval,fex,label="Function")
    plt.plot(xeval,yeval_Lagrange,label="Lagrange Interp.")
    plt.plot(xeval,yeval_Hermite,label="Hermite Interp.")

    plt.xlabel("x")
    plt.ylabel("f(x)/y")
    plt.title("Comparison of Interpolation Methods")
    plt.legend()
    plt.show()

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N+1)
    
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)

def eval_hermite(xeval,xint,yint,ypint,N):
    ''' Evaluate all Lagrange polynomials'''
    lj = np.ones(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])
    ''' Construct the l_j'(x_j)'''
    lpj = np.zeros(N+1)
#    lpj2 = np.ones(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
#              lpj2[count] = lpj2[count]*(xint[count] - xint[jj])
              lpj[count] = lpj[count]+ 1./(xint[count] - xint[jj])  
    yeval = 0.
    for jj in range(N+1):
       Qj = (1.-2.*(xeval-xint[jj])*lpj[jj])*lj[jj]**2
       Rj = (xeval-xint[jj])*lj[jj]**2
#       if (jj == 0):
#         print(Qj)    
#         print(Rj)
#         print(Qj)
#         print(xeval)
 #        return
       yeval = yeval + yint[jj]*Qj+ypint[jj]*Rj  
    return(yeval)

driver()