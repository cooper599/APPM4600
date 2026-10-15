## APPM4600 Lab 8 
# Cooper Wark

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv 
from numpy.linalg import norm

def driver():
    ''' original fucniiotn and interal'''
    # f = lambda x: np.exp(x)
    # a = 0
    # b = 1
    f = lambda x: 1/(1+(10*x)**2)
    a = -1
    b = 1
    
    ''' create points you want to evaluate at'''
    Neval = 100
    xeval =  np.linspace(a,b,Neval)
    ''' number of intervals'''
    Nint = 10
    '''evaluate the linear spline'''
    yeval = eval_lin_spline(xeval,Neval,a,b,f,Nint)
    ''' evaluate f at the evaluation points'''
    fex = f(xeval)
      
    plt.plot(xeval,fex,'ro-',label='exact function')
    plt.plot(xeval,yeval,'bs--',label='natural spline') 
    plt.legend()
    plt.show()
     
    err = abs(yeval-fex)
    plt.figure() 
    plt.semilogy(xeval,err,'ro--',label='absolute error')
    plt.legend()
    plt.show()

    # Experiments 1) }}f0S}}_inf / ||f||_inf
    # 2) ||f-S||_1 / ||f||_1
    # 3) ||f-S||_2 / ||f||_2
    dif = fex - yeval
    err_inf = norm(fex - yeval,np.inf)
    err_1 = norm(fex - yeval,1)
    err_2 = norm(fex - yeval,2)

    f_inf = norm(fex,np.inf)
    f_1 = norm(fex,1)
    f_2 = norm(fex,2)

    rel_inf = err_inf/f_inf
    rel_1 = err_1/f_1
    rel_2 = err_2/f_2

    print("Location global max err: ", rel_inf)
    print("Location 1 norm err: ", rel_1)
    print("Location 2 norm err: ", rel_2)



def  eval_lin_spline(xeval,Neval,a,b,f,Nint):
    '''create the intervals for piecewise approximations'''
    xint = np.linspace(a,b,Nint+1)
    '''create vector to store the evaluation of the linear splines'''
    yeval = np.zeros(Neval) 
    for jint in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        '''let n denote the length of ind'''
        '''temporarily store your info for creating a line in the interval of 
         interest'''
        a1= xint[jint]
        fa1 = f(a1)
        b1 = xint[jint+1]
        fb1 = f(b1)

        if jint == Nint - 1:
            ind = np.where((xeval >= a1) & (xeval <= b1))[0]
        else:
            ind = np.where((xeval >= a1) & (xeval < b1))[0]
        n = len(ind)

        
        for kk in range(n):
            '''use your line evaluator to evaluate the lines at each of the points 
            in the interval'''
            '''yeval(ind(kk)) = call your line evaluator at xeval(ind(kk)) with 
            the points (a1,fa1) and (b1,fb1)'''
            alpha = xeval[ind[kk]]
            yeval[ind[kk]] = line([a1,fa1],[b1,fb1],alpha)
    return yeval

def line(x0,x1,alpha):
    # x0 = [x0,y0], x1 = [x1,y1]
    lhat_0 = (alpha-x1[0])/(x0[0]-x1[0])
    lhat_1 = (alpha-x0[0])/(x1[0]-x0[0])
    p1 = x0[1]*lhat_0 + x1[1]*lhat_1
    return p1

driver()