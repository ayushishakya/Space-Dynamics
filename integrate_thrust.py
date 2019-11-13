#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 03:14:08 2019

@author: ayushi
"""
import scipy as sci
import scipy.integrate
#defining constants
pi = 3.14159
deg = pi/180
mu = 398600
RE = 6378
g0 = 9.807

def integrate_thrust(y0, t0, t_burn):

    def F(w,t,T, Isp):
        x = w[0]
        y = w[1]
        z = w[2]

        vx = w[3]
        vy = w[4]
        vz = w[5]

        m = w[6]

        r = sci.linalg.norm([x, y, z])
        v = sci.linalg.norm([vx, vy, vz])
        ax = -mu*x/r**3 + T/m*vx/v
        ay = -mu*y/r**3 + T/m*vy/v
        az = -mu*z/r**3 + T/m*vz/v
        mdot = -T*1000/g0/Isp

        dfdt = sci.array([vx, vy, vz, ax, ay, az, mdot])
        return dfdt

    time_span=sci.linspace(t0, t_burn, 5000)
    solution = sci.integrate.odeint(F, y0, time_span, args=(T, Isp))
    return solution

if __name__=="__main__":
    r0 = sci.array([RE+480, 0, 0], dtype="float64")
    v0 = sci.array([ 0, 7.7102, 0], dtype="float64")
    t0 = 0
    t_burn = 261.1127
    m0 = 2000
    T = 10
    Isp = 300

    yinit = sci.array([r0, v0, m0])
    yinit = sci.hstack(yinit)
    sol = integrate_thrust(yinit, t0, t_burn)
    print(sol)
'''
    r1 = sol[-1,:3]
    v1 = sol[-1,3:6]
    m1 = sol[-1,6]
    coe = coe_from_sv(r1,v1,mu)
    e = coe[1]
    TA = coe[5]
    a = coe[6]

    print("\nBefore ignition")
    print("\n Mass = ", m0)
    print("\n State vector:")
    print("\n r = ", r0(1), r0(2), r0(3))
    print("\n Radius =", sci.linalg.norm(r0))
    print("\n v =", v0(1), v0(2), v0(3))
    print("\n Speed =", sci.linalg.norm(v0))
    print("\nThrust =", T)
    print("\nBurn time =", t_burn)
    print("\nMass after burn =", m1)
    print("\nEnd-of-burn-state vector:")
    print("\n r =", r1(1), r1(2), r1(3))
    print("\n Radius =", sci.linalg.norm(r1))
    print("\n v =", v1(1), v1(2), v1(3))
    print("\n Speed =", sci.linalgnorm(v1))
    print("\nPost-burn trajectory:")
    print("\n Eccentricity =", e)
    print("\n Semimajor axis =", a)
'''