#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 13:09:59 2019

@author: ayushi
Hohmann transfer
1. co-centric circular orbits
2. co-axial elliptic orbits
"""
import matplotlib.pyplot as plt
import math
import numpy

mu = 3.99e14  #Earth

def hohmann_circle(r1, r2):

    a = (r2 + r1)/2
    vr1 = math.sqrt(mu/r1)
    vr2 = math.sqrt(mu/r2)

    vp = math.sqrt(mu*(2/r1 - 1/a))
    va = math.sqrt(mu*(2/r2 - 1/a))

    v1 = vp - vr1
    v2 = vr2 - va

    delV = [v1, v2]

    e = (vp - va)/(vp + va)
    print("Eccentricity: ", e)
    b = a*math.sqrt(1 - e**2)
    print("Semi minor axis", b)

    theta = [i for i in numpy.arange(0,3.14, 1e-3)]
    xt = [a*math.cos(angle) - r1 for angle in theta]
    yt = [b*math.sin(angle) for angle in theta]

    x1 = [r1*math.cos(angle) for angle in theta]
    y1 = [r1*math.sin(angle) for angle in theta]

    x2 = [r2*math.cos(angle) for angle in theta]
    y2 = [r2*math.sin(angle) for angle in theta]

    plt.plot(xt,yt)
    plt.plot(x1,y1)
    plt.plot(x2,y2)
    plt.xlabel('x - axis')
    plt.ylabel('y - axis')
    plt.title("Hohmann transfer")
    plt.show()

    return delV

def hohmann_ellipse(h1, r1, h2, r2):

    a = (r1 + r2)/2

    v1 = h1/r1
    v2 = h2/r2

    vp = math.sqrt(mu*(2/r1 - 1/a))
    va = math.sqrt(mu*(2/r2 - 1/a))

    vt1 = vp - v1
    vt2 = v2 - va

    delV = [vt1, vt2]

    e = (vp - va)/(vp + va)
    print("Eccentricity: ", e)
    b = a*math.sqrt(1 - e**2)
    print("Semi minor axis", b)


    return delV

if __name__ == '__main__':

    choice = input("1. Hohmann transfer between circular orbits \n2. Hohmann transfer between elliptical orbits \nEnter your choice:  ")

    if choice == 1:

        r1 = 480e3 + 6378e3#input("radius of smaller orbit")
        r2 = 16000e3 + 6378e3#input("radius of larger orbit")

        delV = hohmann_circle(r1, r2)
        print("Velocity boost required for desired trajectory transfer (km/s): ", delV[0]/1e3, delV[1]/1e3)

    elif choice == 2:

        h1 = input("specific angular momentum of inner elliptical orbit")
        rp1 = input("perigee radius of inner elliptic orbit")
        ra1 = input("apogee radius of inner elliptic orbit")

        h2 = input("specific angular momentum of outer elliptical orbit")
        rp2 = input("perigee radius of outer elliptic orbit")
        ra2 = input("apogee radius of outer elliptic orbit")

        print("Hohmann transfer from perigee of inner elliptical orbit to apogee of outer elliptical orbit")
        delV = hohmann_ellipse(h1, rp1, h2, ra2)
        print("Velocity boost required for desired trajectory transfer (km/s): ", delV[0]/1e3, delV[1]/1e3)


        print("Hohmann transfer from apogee of inner elliptical orbit to perigee of outer elliptical orbit")
        delV = hohmann_ellipse(h1, ra1, h2, rp2)
        print("Velocity boost required for desired trajectory transfer (km/s): ", delV[0]/1e3, delV[1]/1e3)
