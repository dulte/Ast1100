# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 20:13:18 2016

@author: Daniel
"""

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pylab
import planets as p

"""
def integrate(x,v,a,dt):
    #Euler-Cromer
    m = 3.e-6
    v += a(x)*dt/m
    x += v*dt
    
    return x,v
    
def gravitational_force(x):
    G = 4*np.pi**2   
    ms = 1
    mp = 3.e-6

    try:    
        
        r = norm(x)
    except:
        r = x
    
    return (-G*mp*ms/r**3)*x
    
def loop(T,dt,x0,v0):
    steps = int(float(T)/dt)
    x = np.zeros([steps,2])
    x[0,:] = x0
    
    v = np.zeros([steps,2])
    v[0,:] = v0
    
    for i in range(1,steps):
        x[i,:],v[i,:] = integrate(x[i-1,:],v[i-1,:],gravitational_force,dt)
    
    return x,v
  

"""
    


def main():
    sun_mass = 1
    earth_mass = 3.e-6
    earth_pos0 = np.array([1,0])
    earth_vel0 = np.array([0,2.3*np.pi])
    sun_pos = np.array([0,0])
    
    
    venus_mass = 2.e-6
    venus_pos0 = np.array([.5,0])
    venus_vel0 = np.array([1,3*np.pi])
    
    
    
    N = 2
    dt = 1./20000
    steps = int(float(N)/dt)
    
    
    earth = p.planet(earth_pos0,earth_vel0,earth_mass,N,dt)
    venus = p.planet(venus_pos0,venus_vel0,venus_mass,N,dt)
    
    for i in range(1,steps):
        earth.integrate(i)
        venus.integrate(i)
        
    
    xe,ve,te = earth.return_data()
    xv,vv,tv = venus.return_data()
    
    
    
    """
    pylab.ion()
    
    limits = [-1.1,1.1]
    
    fig = pylab.figure() 
    ax = fig.add_subplot(111 , autoscale_on=False ,  xlim=limits ,ylim=limits )
    
    
    line , = ax.plot(x[0,0],x[0,1],"ro")
    pylab.show()
    
    for i in xrange(len(te)):
        line.set_data(x[i,0],x[i,1])
        
        try:
            pylab.draw()
        except:
            print "her"
            break
    print "der"
           
    """ 

    
    
    plt.plot(sun_pos[0],sun_pos[1],"yo")
    plt.plot(xe[:,0],xe[:,1],"g")
    plt.plot(xv[:,0],xv[:,1],"r")
    plt.axis("equal")
   

main()
    
    
    
    



        
    
    