#!/usr/bin/env python
A=100000.0
Ea=52000.0
R=8.3145

from math import *

def get_residual_activity(T,t):
    k=exp(-Ea/(R*T))
    perc=100.0*exp(-k*t)
    perc=k
    return perc
    

t=15*60
xs=[]
ys=[]
for T in range(20,100,5):
    xs.append(T+273.14)
    ys.append(get_residual_activity(float(T)+273.15,t))
A=A*7000.0
Ea=Ea+25000.0
zs=[]
for T in range(20,100,5):
    #xs.append(T)
    zs.append(get_residual_activity(float(T)+273.15,t))
import pylab
#pylab.plot(xs,ys,'ro-')
pylab.plot(xs,zs,'bo-',label='steep')
pylab.legend()
pylab.xlabel('Temperature (C)')
pylab.ylabel('Residual activity (%)')
pylab.savefig('test.png',dpi=300)
pylab.show()