#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 01:01:21 2019

@author: ayushi
"""
import scipy.optimize
import scipy as sci
import math
import coe_from_sv

#defining constants
pi = 3.14159

def lambert(R1, R2, t, string, mu):
    #print("R1, R2: ", R1, R2)
    r1 = sci.linalg.norm(R1)
    r2 = sci.linalg.norm(R2)

    c12 = sci.cross(R1, R2) #cross product
    r12 = sci.dot(R1, R2)
   #print("cross, dot:", c12, r12)

    theta = math.acos(r12/(r1*r2))
    if (string=="pro" and c12[2] < 0) or (string=="retro" and c12[2] >=0):
        theta = 2*pi - theta

    A = math.sin(theta)*math.sqrt(r1*r2/(1-math.cos(theta)))
   #print("theta, A:", theta, A)

    def C(z):
        if z > 0:
            c = (1 - math.cos(math.sqrt(z)))/z
        elif z < 0:
            c = (math.cosh(math.sqrt(-z)) - 1)/(-z)
        else:
            c = 1/2
        return c

    def S(z):
        if z > 0:
            s = (math.sqrt(z) - math.sin(math.sqrt(z)))/(math.sqrt(z))**3
        elif z < 0:
            s = (math.sinh(math.sqrt(-z)) - math.sqrt(-z))/(math.sqrt(-z))**3
        else:
            s = 1/6
        return s

    def Y(z):
        return r1 + r2 + A*(z*S(z) - 1)/math.sqrt(C(z))

    def F(z):
        s = S(z)
        y = Y(z)
        c = C(z)
        val = A*math.sqrt(y) + s*(y/c)**1.5 - math.sqrt(mu)*t

        return val

    def Fprime(z):
        if z!=0:
            y = Y(z)
            c = C(z)
            s = S(z)
            dum = (y/c)**1.5*(1/2/z*(c - 3*s/2/c) + 3*s**2/4/c) + A/8*(3*s*math.sqrt(y)/c + A*math.sqrt(c/y))
        else:
            y = Y(0)
            dum = math.sqrt(2)*y**1.5/40* + A/8*(math.sqrt(y) + A/math.sqrt(2*y))

        return dum

    root = scipy.optimize.newton(F, 0, Fprime, args=(), tol=1e-8, maxiter=5000)

    if root <0:
        print("Hyperbola\n")
    elif root ==0:
        print("Parabola")
    else:
        print("Ellipse")

    ysol = Y(root)

    f = 1 - ysol/r1
    g = A*math.sqrt(ysol/mu)
    gdot = 1 - ysol/r2
   #print("z, f, g, gdot: ", root, f, g, gdot)
    #velocities
    V1 = (R2 - f*R1)/g
    V2 = (gdot*R2 - R1)/g
   #print("V1, V2: ", V1, V2)

    return sci.array([V1, V2], dtype="float64")


if __name__=="__main__":
    mu = 1.327124e11#398600
    deg = pi/180
    t = 20.346053*24*60*60#3600
    '''
    e = 0.4
    h = 67232
    deg = pi/180
    theta1 = 45*deg
    theta2 = 190.57*deg
    R1 = h**2/mu/(1 + e * math.cos(theta1))*sci.array([math.cos(theta1),math.sin(theta1),0])
    R2 = h**2/mu/(1 + e * math.cos(theta2))*sci.array([math.cos(theta2),math.sin(theta2),0])
    '''
    R1 = sci.array([-165066454.034 , 144849637.245 , -3423998.596])
    R2 = sci.array([-194482766.201 , 118099855.449 , -2847890.309])
    string = "pro"
    vel = lambert(R1, R2, t, string, mu)
    coe = coe_from_sv.coe_from_sv(R1, vel[0], mu)
    print("\n Angular momentum, h (km^2/s) = ", coe[0])
    print("\n Eccentricity, e =", coe[1])
    print("\n Right ascension, omega (deg) =", coe[2]/deg)
    print("\n Inclination , i(deg) =", coe[3]/deg)
    print("\n Argument of perigee, w (deg) =", coe[4]/deg)
    print("\n True anomaly, theta (deg) =", coe[5]/deg)
    print("\n Semimajor axis, a (km): = ", coe[6])
