#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 04:44:44 2019

@author: ayushi
"""
import coe_from_sv
import jdcal
import de421
from jplephem import Ephemeris
import month_planet_names

def planet_elements_and_sv(planet_id, year, month, day, hour, minute, second):
    global mu
    mu = 1.327124e11

    #calculating julian time
    j0 = sum(jdcal.gcal2jd(year, month, day))
    ut = (hour+minute/60+second/3600)/24

    jd = j0+ut

    #determining state vectors
    eph = Ephemeris(de421)
    planet = month_planet_names.planet_name(planet_id)
    if planet=="earth":
        barycenter = eph.position_and_velocity('earthmoon', jd)
        moonvector = eph.position_and_velocity('moon', jd)
        r = barycenter[0] - eph.earth_share*moonvector[0]
        v = barycenter[1] - eph.earth_share*moonvector[1]

    elif planet=="moon":
        barycenter = eph.position_and_velocity('earthmoon', jd)
        moonvector = eph.position_and_velocity('moon', jd)
        r = barycenter[0] - eph.moon_share*moonvector[0]
        v = barycenter[1] - eph.moon_share*moonvector[1]

    else:
        r, v = eph.position_and_velocity(planet, jd)

    r = r.flatten()
    v = v.flatten()/(24*3600)

    coe = coe_from_sv.coe_from_sv(r, v, mu)

    return [coe, r, v, jd]