## APPM 4600, Lab 2 Code
# Cooper Wark

############################################# 
"""
Copyright (C) 2025  Adrianna M. Gillman

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
############################################# 

import numpy as np

def driver():

# Prelab stuff SECTION 2
     f1 = lambda x: (10/(x+4))**(1/2)
     p = 1.3652300134140976 # actual fixed point

     Nmax = 100
     tol = 1e-10

     p0 = 1.5
     [xstar,ier, x_arr, count] = fixedpt(f1,p0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('f1(xstar):',f1(xstar))
     print('Error message reads:',ier)
     print("Num iterations: ", count)

     [alpha,const] = computeConvergence(x_arr, p)
     print("Alpha = ", alpha) # alpha = 0.999836
     print("Const = ", const)
     '''
     Part 2 questions; num iteration = 12 using abs tolerance
     alpha = 1.0000003318338138
     Constant = 0.1272300010325035
     '''

     ## SECTION 3 STUFF
     # Section 3.2
     convertedVec = seqToVec(x_arr, tol, Nmax) # Aitkins acceleration conversion
     [alpha,const] = computeConvergence(convertedVec,p) # Recompute alpha and const
     print("Alpha = ", alpha)

'''
Write a subroutine that takes in a sequence of approximations and returns a vector of the
approximations. It should also have tolerance and max number of iterations as input.
• Apply Aitken’s ∆2 method to the sequence created by the fixed point iteration in the before
lab exercise set. Determine if the convergence is in fact faster than the fixed point iteration.
Can you figure out the order of convergence?
'''
def seqToVec(p_seq, tol, Nmax):
    nmax = len(p_seq) - 2 # max index for the conversion
    p_vec = np.array((nmax,1)) # preallocate
    # loop trhough to create vector of new approximations
    for i in range(nmax):
        print("i",i)
        num = (p_seq[i+1]-p_seq[i])**2
        den = p_seq[i+2] - 2*p_seq[i+1] + p_seq[i]
        p_vec[i] = p_seq[i] - num/den
    return p_vec


# Fixed pt calculation routine
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''
    x_arr = np.zeros((Nmax,1))

    count = 0
    while (count <Nmax):
       x_arr[count] = x0
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier, x_arr, count]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier, x_arr, count]
    

# Function to compute convergence
def computeConvergence(pk_arr, p):
    # modify pk_arr to be list of only the relevant values
    modpk = pk_arr[pk_arr != 0]
    # Convert middle term but not last to integer, should probably have a check
    mid_k = np.floor(len(modpk)/2).astype(int) 

    eps_k = abs(modpk[mid_k]-p)
    eps_km1 = abs(modpk[mid_k-1]-p)
    eps_kp1 = abs(modpk[mid_k+1]-p)
    num = np.log(abs(eps_kp1/eps_k))
    den = np.log(abs(eps_k/eps_km1))
    alpha = num/den
    const = abs(eps_k)/(abs(eps_km1))**alpha
    return [alpha,const]

driver()