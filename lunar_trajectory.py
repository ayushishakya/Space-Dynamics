'''
Co-planar lunar trajectory
We will assume that the geocentric
trajectory is direct and that lunar arrival occurs prior to apogee of the
geocentric orbit.

subscript 2 indicates initial conditions wrt moon
'''

import math

#canonical units for earth satellites
DU = 6378.1e3
TU = 806.80
DUTU = 7.90538e3

pi = pi = 3.14159
mE = 5.974e24 #mass of earth
mM = mE/81.30 #mass of moon
mS = 1.989e30
D = 384400e3/DU #mean distance between earth and moon
mu = 3.986e14 /(DU*(DUTU**2))
mu_m = 4.9048695e12 /(DU*(DUTU**2))



def lunar_trajectory(r0, v0, phi, lambda1):

    #sphere of influence radius
    RS = D*(mM/mE)**(2/5)
    print("Radius of sphere of influence of moon: *(DU)", RS)


    #energy and angular momoentum
    E = (v0**2)/2 - mu/r0
    #print("Energy (*(DUTU**2)): ", E)
    h = r0*v0*math.cos(phi)
    #print("momentum (*DU*DUTU): ", h)

    #radius at lunar arrival
    r1 = math.sqrt(D**2 + RS**2 - 2*D*RS*math.cos(lambda1))
    print("Radius at lunar arrival, ", r1*DU)

    #speed at lunar arrival
    v1 = math.sqrt(2*(E + mu/r1))

    print("Velocity at lunar arrival, ", v1*(DUTU))

    #flight path angle at arrival
    phi1 = math.acos(h/(r1*v1))
    #print("phi1 ", phi1)


    gamma1 = math.asin(RS*math.sin(lambda1)/r1)
    #print("gamma1 ", gamma1)
    #geocentric trajectory
    p = (h**2)/mu
    a = -mu/(2*E)
    e = math.sqrt(1- p/a)
    #print("orbit elements ", p, a , e)

    nu0 = math.acos((p-r0)/(r0*e))
    nu1 = math.acos((p-r1)/(r1*e))
    #print("nu0, nu1 ", nu0, nu1)

    #eccentric anomalies
    E0 = math.acos((e+math.cos(nu0))/(1+ e*math.cos(nu0)))
    E1 = math.acos((e+math.cos(nu1))/(1+ e*math.cos(nu1)))
    #print("Eccentric anomalities: ", E0, E1)

    #time of flight
    tof = math.sqrt((a**3)/mu)*((E1 - e*math.sin(E1)) - (E0 - e*math.sin(E0)))
    print("time of flight ", tof)

    #angular velocity of moon in its orbit
    w_m = 2.649e-3 #based on simplified model of earth moon system

    #phase angle at departure
    gamma0 = nu1 - nu0 - gamma1 - w_m*tof
    print("phase angle of departure ", gamma0)

    #Conditions at patch point
    #trajectory inside moon's sphere of influence
    vm = 1.018e3/DUTU #velocity of moon relative to earth
    r2 = RS
    v2 = math.sqrt(v1**2 + vm**2 - 2*v1*vm*math.cos(phi1-gamma1))
    print("patch position and velocity ", r2*DU, v2*DUTU)

    # direction of initial selenocentric velocity
    epsilon2 = math.asin((vm*math.cos(lambda1) - v1*math.cos(gamma1 + lambda1 - phi1))/v2)
    #print("epsilon2 ", epsilon2)
    #selenocentric arrival orbit

    #energy and momentum relative to moon
    E_m = (v2**2)/2 - mu_m/r2
    h_m = r2*v2*math.cos(epsilon2)
    print("E, h", E_m, h_m*DU*DUTU)

    #geocentric trajectory
    p_m = (h_m**2)/mu_m
    e_m = math.sqrt(1 + (2*E_m*(h_m**2))/(mu_m**2))
    #print("Orbital elements: ", p_m, e_m)

    #positon and velocity at selenocentric orbit
    rp = p_m/(1+e_m)
    vp = math.sqrt(2*(E_m + mu_m/rp))

    print("positon and velocity at selenocentric orbit", rp*DU, vp*DUTU)


if __name__ == '__main__':

    #injecttion conditions based on canonical units based on Earth
    r0 = 1.05
    v0 = 1.372
    phi0 = 0 #flight path angle
    lambda1 = 30*pi/180 #point at which the geocentric trajectory crosses the lunar sphere of influence.

    lunar_trajectory(r0, v0, phi0, lambda1)
