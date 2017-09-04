## Curve Fitting

```python {cmd=true .line-numbers matplotlib=true}}

import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
import pandas as pd
#%matplotlib inline

#functions for fitting
def Quadratic_RPM(x, a, b, c):
    return a + b*x + c*x**2

def Quadratic_RPM_2nd(x, a, b):
    return a + b*x

def Linear_Torque(x, a):
    return a*x

def GetCoef(D,R,fun):
    popt,pcov=opt.curve_fit(fun,D,R)
    return popt

#read data for RPM
df = pd.read_csv("./data.csv")
#print df.head()

DP = df["DP"].values
RPM = np.array([df["RPM1"].values, df["RPM2"].values, df["RPM3"].values])

#max min medium flow rate
Q = np.array([120,100,80])

#coefs for the first fitting for A B C in f(Q)=A+BQ+CQ^2
coef={}
for i in range(3):
    coef[i]=GetCoef(DP,RPM[i],Quadratic_RPM)   



coefA=GetCoef(DP,RPM[i],Linear_Torque)   

print coefA    

coefB={}

coefB=GetCoef(DP,RPM[i],Quadratic_RPM_2nd)   

print coefB    