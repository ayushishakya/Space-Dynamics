#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 00:38:50 2019

@author: ayushi
"""
import scipy as sci
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.integrate

#Define universal gravitation constant
G = 6.67408e-11

#Reference quantities
m0 = 1.989e+30 #kg #mass of the sun
r0 = 5.326e+12 #m #distance between stars in Alpha Centauri
v0 = 30000 #m/s #relative velocity of earth around the sun
t0 = 79.91*365*24*3600*0.51 #s #orbital period of Alpha Centauri

#Net constants
k1 = G*t0*m0/(r0**2*v0)
k2 = v0*t0/r0

#Define masses
m1 = 1.1 #Alpha Centauri A
m2 = 0.907 #Alpha Centauri B
m3 = 1.0 #Third Star

#Define initial position vectors
r1 = sci.array([-0.5,0,0],dtype="float64")
r2 = sci.array([0.5,0,0],dtype="float64")
r3 = sci.array([0,1,0],dtype="float64")

#Find Centre of Mass
r_com = (m1*r1+m2*r2+m3*r3)/(m1+m2+m3)

#Define initial velocities
v1 = sci.array([0.01,0.01,0],dtype="float64")
v2 = sci.array([-0.05,0,-0.1],dtype="float64")
v3 = sci.array([0,-0.01,0],dtype="float64")

#velocity of COM
vcom = (m1*v1+m2*v2+m3*v3)/(m1+m2+m3)

def ThreeBodyEquations(w,t,G,m1,m2,m3):
    r1=w[:3]
    r2=w[3:6]
    r3=w[6:9]
    v1=w[9:12]
    v2=w[12:15]
    v3=w[15:18]
    r12=sci.linalg.norm(r2-r1)
    r13=sci.linalg.norm(r3-r1)
    r23=sci.linalg.norm(r3-r2)

    dv1dt=k1*m2*(r2-r1)/r12**3+k1*m3*(r3-r1)/r13**3
    dv2dt=k1*m1*(r1-r2)/r12**3+k1*m3*(r3-r2)/r23**3
    dv3dt=k1*m1*(r1-r3)/r13**3+k1*m2*(r2-r3)/r23**3
    dr1dt=k2*v1
    dr2dt=k2*v2
    dr3dt=k2*v3
    r12_derivs=sci.concatenate((dr1dt,dr2dt))
    r_derivs=sci.concatenate((r12_derivs,dr3dt))
    v12_derivs=sci.concatenate((dv1dt,dv2dt))
    v_derivs=sci.concatenate((v12_derivs,dv3dt))
    derivs=sci.concatenate((r_derivs,v_derivs))
    return derivs

#Package initial parameters
init_params=sci.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
init_params=init_params.flatten() #Flatten to make 1D array
time_span=sci.linspace(0,20,500) #20 orbital periods and 500 points
#Run the ODE solver
import scipy.integrate
three_body_sol=sci.integrate.odeint(ThreeBodyEquations,init_params,time_span,args=(G,m1,m2,m3))

r1_sol=three_body_sol[:,:3]
r2_sol=three_body_sol[:,3:6]
r3_sol=three_body_sol[:,6:9]

#Create figure
fig=plt.figure(figsize=(15,15))

#Create 3D axes
ax=plt.gca(projection="3d")

#Plot the orbits
ax.plot(r1_sol[:,0],r1_sol[:,1],r1_sol[:,2],color="darkblue")
ax.plot(r2_sol[:,0],r2_sol[:,1],r2_sol[:,2],color="tab:red")
ax.plot(r3_sol[:,0],r3_sol[:,1],r3_sol[:,2],color="tab:green")

#Plot the final positions of the stars
ax.scatter(r1_sol[-1,0],r1_sol[-1,1],r1_sol[-1,2],color="darkblue",marker="o",s=100,label="Alpha Centauri A")
ax.scatter(r2_sol[-1,0],r2_sol[-1,1],r2_sol[-1,2],color="tab:red",marker="o",s=100,label="Alpha Centauri B")
ax.scatter(r3_sol[-1,0],r3_sol[-1,1],r3_sol[-1,2],color="tab:green",marker="o",s=100,label="Third Star")

#Add a few more bells and whistles
ax.set_xlabel("x-coordinate",fontsize=14)
ax.set_ylabel("y-coordinate",fontsize=14)
ax.set_zlabel("z-coordinate",fontsize=14)
ax.set_title("Visualization of orbits of stars in a two-body system\n",fontsize=14)
ax.legend(loc="upper left",fontsize=14)

plt.show()