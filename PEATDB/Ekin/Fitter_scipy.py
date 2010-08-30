#!/usr/bin/env python

import scipy.optimize

def funktion(variables):
    print variables
    
    return variables[0]**2*variables[1]*0.554


def callback(vk):
    print vk
    return

start_values=[44.0,100.0]

junk=scipy.optimize.fmin(funktion,start_values,xtol=1E-10,callback=callback)

print 'junk',junk
print funktion(junk)
