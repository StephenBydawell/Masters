# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 11:21:21 2016

@author: stephen
"""

import numpy as np
from iapws import IAPWS97
import sympy

P = 19.0
T = 813
d = 0.250
t = 0.029
Ai = np.pi*d
Ao = np.pi*(d+2*t)
ri = d/2.0
ro = ri + t
rins = ro + 0.160

steam=IAPWS97(P=19, T=713)               #steam with known P and T

fr = 260

#v = fr/steam.rho

v = 12.5

Re = (steam.rho*v*d)/steam.mu

Nu = 0.023*Re**0.8*steam.Prandt**0.4

h = (steam.k/d)*Nu
k = 26.4
kins = 0.033

print(v)
print(h)

opt = 1

#forward

if opt == 1:

    R1 = 1/(Ai*h)
    R2 = np.log(ro/ri)/(2*np.pi*k)
    R3 = np.log(rins/ro)/(2*np.pi*kins)
    
    Tinf = 30+273
    To = sympy.Symbol('To')
    T1 = sympy.Symbol('T1')
    Ti = sympy.Symbol('Ti')
    
    ho = 1.32*((To-Tinf)/(2*rins))**(1.0/4.0)
    R4 = 1/(ho*2*np.pi*rins)
    
    Ti = sympy.solve((T-Ti)/R1 - (Ti-T1)/R2, Ti)[0]
    T1 = sympy.solve((Ti-T1)/R2 - (T1-To)/R3, T1)[0]
    
    term1 = (T1-To)/R3
    term2 = 2*np.pi*rins*(1.32/(2*rins)**(0.25))*(To-Tinf)**(1.25)
    
    To = sympy.solve(term1 - term2, To)[0]
    ho = ho.subs('To',To)
    T1 = T1.subs('To',To)
    Ti = Ti.subs('T1',T1)
    
    print(To-273)
    print(T1)
    print(Ti)
    print(Ti-T1)
    
    heff = (1/(ho*2*np.pi*rins) + (np.log(rins/ro))/(2*np.pi*kins))**(-1)
    heff = (1/(ho) + (2*np.pi*rins*np.log(rins/ro))/(2*np.pi*kins))**(-1)

#Reverse

if opt == 2:
    
    Tinf = 30+273
    To = sympy.Symbol('To')
    Ti = sympy.Symbol('Ti')
    T1 = Ti - 4
    kins = sympy.Symbol('kins')
    
#    ho = 6.51182459352817
    ho = 1
    R1 = 1/(Ai*h)
    R3 = np.log(rins/ro)/(2*np.pi*kins)
    R2 = np.log(ro/ri)/(2*np.pi*k)
    R4 = 1/(ho*2*np.pi*rins)
    
    Ti = sympy.solve((T-Ti)/R1 - (Ti-T1)/R2, Ti)[0]
    T1 = Ti - 4
    
    kins = sympy.solve((T1-To)/R3 - (To-Tinf)/R4, kins)[0]
    R3 = R3.subs('kins', kins)
    To = sympy.solve((Ti-T1)/R3 - (T1-To)/R4, To)[0]
    kins = kins.subs('To',To)
    
    print(To-273)
    print(T1)
    print(Ti)
    print(kins)

    heff = (1/(ho*2*np.pi*rins) + (np.log(rins/ro))/(2*np.pi*kins))**(-1)
    
    
#    
#    Ti = sympy.solve((T-Ti)/R1 - (Ti-To)/R2, Ti)[0]
#    To = Ti-4
#    ho = sympy.solve((Ti-To)/R2 - (To-Tinf)/R4, ho)[0]
    

else:
    R1 = 1/(Ai*h)
    R2 = np.log(ro/ri)/(2*np.pi*k)
    
    
    Tinf = 30+273
    Ti = sympy.Symbol('Ti')
    To = Ti - 4
    ho = sympy.Symbol('ho')
    
#    R4 = 1/(ho*2*np.pi*ro)
    R4 = 1/(Ao*ho)
    
    Ti = sympy.solve((T-Ti)/R1 - (Ti-To)/R2, Ti)[0]
    To = Ti-4
    ho = sympy.solve((Ti-To)/R2 - (To-Tinf)/R4, ho)[0]
    
#    heff = (1/(ho) + (2*np.pi*rins*np.log(rins/ro))/(2*np.pi*kins))**(-1)
    
#    heff = 0.166081738599762
#    heff = 0.144247170734293

