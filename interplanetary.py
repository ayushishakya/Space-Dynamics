#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 04:41:10 2019

@author: ayushi
"""
import planet_elements_and_sv
import lambert
import month_planet_names
import scipy as sci
import math

mu = 1.327124e20
pi = 3.14159
deg = pi/180

def interplanetary(depart, arrive):

    #DEPART
    planet_id = depart[0]
    year = depart[1]
    month = depart[2]
    day = depart[3]
    hour = depart[4]
    minute = depart[5]
    second = depart[6]

    [coe1, Rp1, Vp1, jd1] = planet_elements_and_sv.planet_elements_and_sv(planet_id, year, month, day, hour, minute, second)

    #ARRIVE
    planet_id = arrive[0]
    year = arrive[1]
    month = arrive[2]
    day = arrive[3]
    hour = arrive[4]
    minute = arrive[5]
    second = arrive[6]

    [coe2, Rp2, Vp2, jd2] = planet_elements_and_sv.planet_elements_and_sv(planet_id, year, month, day, hour, minute, second)
    tof = (jd2-jd1)*86400
    #Patched conic assumption:
    R1 = Rp1
    R2 = Rp2

    #departure and arrival, assuming a prograde trajectory:
    [V1, V2] = lambert.lambert(R1, R2, tof, "pro", mu)
    planet1 = [coe1, Rp1, Vp1, jd1]
    planet2 = [coe2, Rp2, Vp2, jd2]
    trajectory = [V1, V2]

    return [planet1, planet2, trajectory]

if __name__=="__main__":

    #Departure
    planet_id = 2
    year = 1996
    month = 10
    day = 7
    hour = 0
    minute = 0
    second = 0
    depart = [planet_id, year, month, day, hour, minute, second]

    #Arrival
    planet_id = 3
    year = 1997
    month = 9
    day = 12
    hour = 0
    minute = 0
    second = 0
    arrive = [planet_id, year, month, day, hour, minute, second]

    [planet1, planet2, trajectory] = interplanetary(depart, arrive)

    coe = planet1[0]
    R1 = planet1[1]
    Vp1 = planet1[2]
    jd1 = planet1[3]

    coe2 = planet2[0]
    R2 = planet2[1]
    Vp2 = planet2[2]
    jd2 = planet2[3]

    V1 = trajectory[0]
    V2 = trajectory[1]
    tof = (jd2-jd1)*86400

    vinf1 = V1 - Vp1
    vinf2 = V2 - Vp2

'''
    print("\n\n Departure:\n")
    print("\n Planet: ", month_planet_names.planet_name(depart[0]))
    print("\n Year : ", depart[1])
    print("\n Month : ", month_planet_names.month_name(depart[2]))
    print("\n Day : ", depart[3])
    print("\n Hour : ", depart[4])
    print("\n Minute: ", depart[5])
    print("\n Second: ", depart[6])
    print("\n\n Julian day: ", jd1)
    print("\n Planet position vector (km) = ",R1)
    print("\n Magnitude = \n", sci.linalg.norm(R1))
    print("\n Planet velocity (km/s) = ",Vp1)
    print("\n Magnitude = \n", sci.linalg.norm(Vp1) )
    print("\n Spacecraft velocity (km/s) = ", V1)
    print("\n Magnitude = \n", sci.linalg.norm(V1))
    print("\n v-infinity at departure (km/s) = ", vinf1 )
    print("\n Magnitude = \n", sci.linalg.norm(vinf1))
    print("\n\n Time of flight =  days\n", tof)
    print("\n\n Arrival:\n")
    print("\n Planet: ", month_planet_names.planet_name(arrive[0]))
    print("\n Year : ", arrive[1])
    print("\n Month : ", month_planet_names.month_name(arrive[2]))
    print("\n Day : ", arrive[3] )
    print("\n Hour : ", arrive[4])
    print("\n Minute: ", arrive[5] )
    print("\n Second: ", arrive[6] )
    print("\n\n Julian day: \n", jd2)
    print("\n Planet position vector (km) = ", R2)
    print("\n Magnitude = \n", sci.linalg.norm(R1))
    print("\n Planet velocity (km/s) = ", Vp2)
    print("\n Magnitude = \n", sci.linalg.norm(Vp2))
    print("\n Spacecraft Velocity (km/s) = ", V2)
    print("\n Magnitude = \n", sci.linalg.norm(V2))
    print("\n v-infinity at arrival (km/s) = ", vinf2)
    print("\n Magnitude = ", sci.linalg.norm(vinf2))
    print("\n\n\n Orbital elements of flight trajectory:\n")
    print("\n Angular momentum (km^2/s) = ", coe[0] )
    print("\n Eccentricity = ",  coe[1] )
    print("\n Right ascension of the ascending node (deg) = ", coe[2]/deg)
    print("\n Inclination to the ecliptic (deg) = ", coe[3]/deg)
    print("\n Argument of perihelion (deg) = ", coe[4]/deg)
    print("\n True anomaly at departure (deg) = ", coe[5]/deg)
    print("\n True anomaly at arrival (deg) = \n", coe2[5]/deg)
    print("\n Semimajor axis (km) = ",  coe[6])

    if coe[1] < 1:
        print("\n Period (days) = ", 2*pi/math.sqrt(mu)*coe[6]**1.5/24/3600)
'''