#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 02:37:37 2019

@author: ayushi
"""
import math
import scipy as sci
import scipy.linalg

#defining constants
pi = 3.14159
eps = 1e-10 #smallest value for eccentricity
deg = pi/180

def coe_from_sv(R, V, mu):

    print(R,V)
    r = sci.linalg.norm(R)
    v = sci.linalg.norm(V)
    vr = sci.dot(R,V)/r

    #print(r, v, vr)

    #angular momentum
    H = sci.cross(R,V)
    h = sci.linalg.norm(H)

    #inclination

    incl = math.acos(H[2]/h)

    N = sci.cross([0,0,1], H)
    n = sci.linalg.norm(N)

    if n!=0:
        RA = math.acos(N[0]/n)
        if N[1] < 0:
            RA = 2*pi - RA
    else:
        RA = 0

    E = ((v**2 - mu/r)*R - r*vr*V)/mu

    e = sci.linalg.norm(E)

    if n!=0:
        if e > eps:
            w = math.acos(sci.dot(N,E)/n/e)
            if E[2] < 0:
                w = 2*pi - w
        else:
            w = 0
    else:
        w = 0

    if e > eps:
        TA = math.acos(sci.dot(E,R)/e/r)
        if vr < 0:
            TA = 2*pi - TA
    else:
        cp = sci.cross(N,R)
        if cp[2] >= 0:
            TA = math.acos(sci.dot(N,R)/n/r)
        else:
            TA = 2*pi - math.acos(sci.dot(N,R)/n/r)


    a = h**2/mu/(1 - e**2)
    coe = sci.array([h, e, RA, incl, w, TA, a], dtype="float64")
    print(coe)
    return coe
'''
if __name__=="__main__":
    mu = 398600
    r = sci.array([-6045 ,-3490, 2500], dtype="float")
    v = sci.array([-3.457, 6.618, 2.533], dtype="float64")
    coe = coe_from_sv(r,v,mu)

    print("\n Angular momentum, h (km^2/s) = ", coe[0])
    print("\n Eccentricity, e =", coe[1])
    print("\n Right ascension, omega (deg) =", coe[2]/deg)
    print("\n Inclination , i(deg) =", coe[3]/deg)
    print("\n Argument of perigee, w (deg) =", coe[4]/deg)
    print("\n True anomaly, theta (deg) =", coe[5]/deg)
    print("\n Semimajor axis, a (km): = ", coe[6])

    if coe[1]<1:
        T = 2*pi/math.sqrt(mu)*coe[6]**1.5
        print("\n Period:")
        print("\n Seconds =", T)
        print("\n Minutes =", T/60)
        print("\n Hours =", T/3600)
        print("\n Days =", T/24/3600)
'''