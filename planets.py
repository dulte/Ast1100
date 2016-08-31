# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 14:21:35 2016

@author: Daniel
"""


import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import pylab


class planet:
    
    
    def __init__(self, x0,v0,m,T,dt):
        self.x0 = x0
        self.v0 = v0
        self.m = m
        self.T = T
        self.dt = dt
        
        self.steps = int(float(self.T)/self.dt)
        
        self.x = np.zeros([self.steps,2])
        self.v = np.zeros([self.steps,2])
        self.t = np.zeros([self.steps])
        
        self.x[0,:] = self.x0
        self.v[0,:] = self.v0
        self.t[0] = 0
        
        self.G = 4*np.pi**2
        
    
    def gravitational_force(self,i):
        ms = 1
        mp = self.m
    
        try:    
            
            r = norm(self.x[i-1,:])
        except:
            r =self.x[i-1,:]
        
        return (-self.G*mp*ms/r**3)*self.x[i-1,:]
        
    def integrate(self,i):
        self.v[i,:] = self.v[i-1,:] + self.gravitational_force(i)*self.dt/self.m
        self.x[i,:] = self.x[i-1,:] + self.v[i,:]*self.dt
        
        self.t[i] += self.dt
        
    def return_data(self):
        return self.x,self.v,self.t
        
        
        
        
        