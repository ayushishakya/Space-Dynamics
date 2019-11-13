#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 04:52:22 2019

@author: ayushi
"""
import math
pi = 3.14159

def kepler_E(e, M):
    error = 1e-8

    if M < pi:
        E = M + e/2
    else:
        E = M - e/2

    ratio = 1

    while abs(ratio) > error:
        ratio = (E - e*math.sin(E) - M)/(1 - e*math.cos(E))
        E = E - ratio
    return E