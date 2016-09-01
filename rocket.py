# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:30:05 2016
@author: Daniel
"""

import random
import numpy as np
import os

class engine:
    def __init__(self, sideSize, holeSize, fuel, temp,dt):
        self.L = float(sideSize)
        self.hole = holeSize
        self.fuel = fuel
        self.T = temp
        self.k = 1.38e-23
        self.dt = dt



        self.m =  2*1.7e-27
        self.N = 1e5

        self.sigma = np.sqrt(self.k*self.T/self.m)

        self.collisions = np.zeros([3,2])
        self.wall_force = np.zeros([3,2])
        self.force_on_top = 0
        self.wall_momentum = np.zeros([3,2])
        self.momentum_gaines = 0
        self.numb_escaping = 0





        self.x = np.zeros([int(self.N),3])
        self.v = np.zeros([int(self.N),3])




        for i in range(int(self.N)):
            for j in range(3):
                self.x[i,j] = -self.L/2+(random.random()*self.L)
                self.v[i,j] = random.gauss(0,self.sigma)#(3./2)*self.k*self.T)



    def integrate(self):


        self.x += self.v*self.dt




    def detect_collision(self,dt):
        self.collisions[:,:] = 0
        self.wall_force[:,:] = 0


        for i in range(int(self.N)):
            for j in range(3):
                if(self.x[i,j] > self.L/2):

                    self.wall_force[j,1] += 2*self.v[i,j]*self.m/dt
                    self.wall_momentum[j,1] += 2*self.v[i,j]*self.m
                    self.v[i,j] *= -1
                    self.x[i,j] = self.L/2

                    self.collisions[j,1] += 1



                elif(self.x[i,j] < -self.L/2):
                    if (j == 1 and not((self.x[i,0] > self.hole/2 or self.x[i,0] < -self.hole/2) and(self.x[i,3] > self.hole/2 or self.x[i,3] < -self.hole/2))):
                        self.momentum_gaines += np.abs(self.v[i,j]*self.m)
                        self.x[i,j] = -3*self.L
                        self.v[i,j] = 0
                        continue
                    else:
                        self.wall_force[j,0] += 2*self.v[i,j]*self.m/dt
                        self.wall_momentum[j,0] += 2*self.v[i,j]*self.m
                        self.v[i,j] *= -1
                        self.x[i,j] = -self.L/2
                        self.collisions[j,0] += 1




        return self.collisions, self.wall_force

    def detect_collision_vec(self):


        for i in range(3):

            if i == 1:
                index_part_through_floor = np.logical_and(self.x[:,i] < -self.L/2., np.logical_and(self.x[:,i-1] < self.hole/2.,self.x[:,i-1] > -self.hole/2.),np.logical_and(self.x[:,i+1] < self.hole/2.,self.x[:,i+1] > -self.hole/2.))
                self.momentum_gaines += np.abs(np.sum(self.v[index_part_through_floor,i])*self.m)
                self.numb_escaping += np.sum(index_part_through_floor)
                self.force_on_top += 2*np.sum(self.v[self.x[:,i] > self.L/2.,i])*self.m/self.dt




            self.v[np.abs(self.x[:,i]) > self.L/2.,i] = -self.v[np.abs(self.x[:,i]) > self.L/2.,i]

            self.x[self.x[:,i] > self.L/2.,i] = self.L/2.
            self.x[self.x[:,i] < -self.L/2.,i] = -self.L/2.

            # if i == 1:
            #     index_part_through_floor = np.logical_and(self.x[:,i] < -self.L/2., np.logical_and(self.x[:,i-1] < self.hole/2.,self.x[:,i-1] > -self.hole/2.),np.logical_and(self.x[:,i+1] < self.hole/2.,self.x[:,i+1] > -self.hole/2.))
            #     self.momentum_gaines += np.abs(np.sum(self.v[index_part_through_floor,i])*self.m)
            #     self.numb_escaping += np.sum(index_part_through_floor)
            #     self.v[index_part_through_floor,i] = np.random.normal(0,self.sigma)
            #     self.x[index_part_through_floor,i] = np.random.uniform(-self.L/2.,self.L/2)
            #
            #     index_not_through_floor = np.logical_and(self.x[:,i] < -self.L/2., np.logical_and(self.x[:,i-1] >= self.hole/2.,self.x[:,i-1] <= -self.hole/2.),np.logical_and(self.x[:,i+1] >= self.hole/2.,self.x[:,i+1] <= -self.hole/2.))
            #     self.v[index_not_through_floor,i] = -self.v[index_not_through_floor,i]
            #     self.x[index_not_through_floor,i] = -self.L/2.
            #
            #     self.force_on_top += 2*np.sum(self.v[self.x[:,i] > self.L/2.,i])*self.m/self.dt
            #     self.v[self.x[:,i] > self.L/2.,i] = -self.v[self.x[:,i] > self.L/2.,i]
            #     self.x[self.x[:,i] > self.L/2.,i] = self.L/2.
            #
            #
            #
            # else:
            #     self.v[np.abs(self.x[:,i]) > self.L/2.,i] = -self.v[np.abs(self.x[:,i]) > self.L/2.,i]
            #
            #     self.x[self.x[:,i] > self.L/2.,i] = self.L/2.
            #     self.x[self.x[:,i] < -self.L/2.,i] = -self.L/2.


    def print_ke(self):
        print "Calculated energy: ", (3./2)*self.k*self.T
        ke = 0
        mid = 0
        for i in range(int(self.N)):
            mid = 0
            for j in range(3):
                mid += self.v[i,j]**2

            ke += 0.5*self.m*mid

        print "Numerical energy: ", ke/self.N

    def print_press(self):
        #print "Trykk: " ,((2*self.momentum_gaines/(self.dt*self.N))/(L**2))
        print "Trykk: " , (self.force_on_top/(self.L**2)) / steps
        print "Utregnet trykk: ", self.N*self.k*self.T/(self.L**3)





interval = 1e-9
steps = 1000
dt = interval/steps
L = 1e-6
H = L/2
e = engine(L,H,10,10000,dt)
xforce = 0
yforce = 0
zforce = 0

xcol = 0
ycol = 0
zcol = 0

xmom,ymom,zmom = 0,0,0

print e.momentum_gaines
e.print_ke()

for i in range(steps):

    e.detect_collision_vec()

    e.integrate()
    print (float(i)/steps)*100, "%            \r",
print ""



print e.momentum_gaines/dt
e.print_press()
e.print_ke()
print e.numb_escaping/steps
print e.numb_escaping
