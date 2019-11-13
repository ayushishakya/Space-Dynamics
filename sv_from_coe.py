#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 04:54:30 2019

@author: ayushi
"""
import math
import scipy as sci
def sv_from_coe(coe,mu):
    h = coe[0]
    e = coe[1]
    RA = coe[2]
    incl = coe[3]
    w = coe[4]
    TA = coe[5]

    rp = (h**2/mu) * (1/(1 + e*math.cos(TA))) * (math.cos(TA)*sci.array([1,0,0]) + math.sin(TA)*sci.array([0,1,0]))
    vp = (mu/h) * (-math.sin(TA)*sci.array([1,0,0]) + (e + math.cos(TA))*sci.array([0,1,0]))

    R3_W = sci.array([[math.cos(RA), math.sin(RA), 0], [-math.sin(RA), math.cos(RA), 0], [0,0,1]])

    R1_i = sci.array([[1, 0, 0], [0, math.cos(incl), math.sin(incl)], [0, -math.sin(incl), math.cos(incl)]])

    R3_w = sci.array([[math.cos(w), math.sin(w), 0], [-math.sin(w), math.cos(w), 0], [0,0,1]])

    Q_pX = (R3_w*R1_i*R3_W).T

    r = Q_pX*rp
    v = Q_pX*vp

    return [r, v]