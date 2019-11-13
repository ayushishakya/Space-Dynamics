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

#Define initial position vectors
r1 = sci.array([-0.5,0,0],dtype="float64")
r2 = sci.array([0.5,0,0],dtype="float64")

#Find Centre of Mass
r_com = (m1*r1+m2*r2)/(m1+m2)

#Define initial velocities
v1 = sci.array([0.01,0.01,0],dtype="float64")
v2 = sci.array([-0.05,0,-0.1],dtype="float64")

#Find velocity of COM
v_com = m1*v1+m2*v2/(m1+m2)

#initial state
yinit = sci.array([r1, r2, v1, v2])
yinit = yinit.flatten()

#A function defining the equations of motion
def TwoBodyEquations(w,t,G,m1,m2):
    r1 = w[:3]
    r2 = w[3:6]
    v1 = w[6:9]
    v2 = w[9:12]
    r  = sci.linalg.norm(r2-r1) #Calculate magnitude or norm of vector
    dv1dt = k1*m2*(r2-r1)/r**3
    dv2dt = k1*m1*(r1-r2)/r**3
    dr1dt = k2*v1
    dr2dt = k2*v2
    r_derivs=sci.concatenate((dr1dt,dr2dt))
    derivs=sci.concatenate((r_derivs,dv1dt,dv2dt))
    return derivs

#Package initial parameters
init_params=sci.array([r1,r2,v1,v2]) #create array of initial params
init_params=init_params.flatten() #flatten array to make it 1D
time_span=sci.linspace(0,8,500) #8 orbital periods and 500 points

#Run the ODE solver
ysol=sci.integrate.odeint(TwoBodyEquations,init_params,time_span,args=(G,m1,m2))

r1_sol = ysol[:, :3]
r2_sol = ysol[:, 3:6]
'''
#Create figure
fig=plt.figure(figsize=(15,15))

#Create 3D axes
ax=plt.gca(projection="3d")

#Plot the orbits
ax.plot(r1_sol[:,0],r1_sol[:,1],r1_sol[:,2],color="darkblue")
ax.plot(r2_sol[:,0],r2_sol[:,1],r2_sol[:,2],color="tab:red")

#Plot the final positions of the stars
ax.scatter(r1_sol[-1,0],r1_sol[-1,1],r1_sol[-1,2],color="darkblue",marker="o",s=100,label="Alpha Centauri A")
ax.scatter(r2_sol[-1,0],r2_sol[-1,1],r2_sol[-1,2],color="tab:red",marker="o",s=100,label="Alpha Centauri B")

#Add a few more bells and whistles
ax.set_xlabel("x-coordinate",fontsize=14)
ax.set_ylabel("y-coordinate",fontsize=14)
ax.set_zlabel("z-coordinate",fontsize=14)
ax.set_title("Visualization of orbits of stars in a two-body system\n",fontsize=14)
ax.legend(loc="upper left",fontsize=14)

plt.show()
'''
#Find location of COM
rcom_sol=(m1*r1_sol+m2*r2_sol)/(m1+m2)
#Find location of Alpha Centauri A w.r.t COM
r1com_sol=r1_sol-rcom_sol
#Find location of Alpha Centauri B w.r.t COM
r2com_sol=r2_sol-rcom_sol

#Create figure
fig=plt.figure(figsize=(15,15))

#Create 3D axes
ax=plt.gca(projection="3d")

#Plot the orbits
ax.plot(r1com_sol[:,0],r1com_sol[:,1],r1com_sol[:,2],color="darkblue")
ax.plot(r2com_sol[:,0],r2com_sol[:,1],r2com_sol[:,2],color="tab:red")

#Plot the final positions of the stars
ax.scatter(r1com_sol[-1,0],r1com_sol[-1,1],r1com_sol[-1,2],color="darkblue",marker="o",s=100,label="Alpha Centauri A")
ax.scatter(r2com_sol[-1,0],r2com_sol[-1,1],r2com_sol[-1,2],color="tab:red",marker="o",s=100,label="Alpha Centauri B")

#Add a few more bells and whistles
ax.set_xlabel("x-coordinate",fontsize=14)
ax.set_ylabel("y-coordinate",fontsize=14)
ax.set_zlabel("z-coordinate",fontsize=14)
ax.set_title("Visualization of orbits of stars in a two-body system\n",fontsize=14)
ax.legend(loc="upper left",fontsize=14)

plt.show()