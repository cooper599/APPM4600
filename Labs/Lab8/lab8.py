## APPM4600 Lab 8 
# Cooper Wark

import numpy as np
from numpy.linalg import inv

def driver():
    x0 = [0,0]
    x1 = [2,2]
    alpha = 1

def line(x0,x1,alpha):
    # x0 = [x0,y0], x1 = [x1,y1]
    lhat_0 = (alpha-x1[0])/(x0[0]-x1[0])
    lhat_1 = (alpha-x0[0])/(x1[0]-x0[0])
    p1 = x0[1]*lhat_0 + x1[1]*lhat_1
    return p1

driver()