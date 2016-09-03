# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:30:05 2016
@author: Daniel
"""
from AST1100SolarSystem import AST1100SolarSystem
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import random
import numpy as np
import matplotlib.pyplot as plt
import os

class engine:
    def __init__(self, sideSize, holeSize, fuel, temp,interval,steps,number_of_particles):
        self.L = float(sideSize)
        self.hole = holeSize
        self.fuel = fuel
        self.T = temp
        self.k = 1.38e-23
        self.dt = interval/steps
        self.steps = steps
        self.interval = interval

        self.launcher_mass = 1100



        self.m =  3.3e-27
        self.N = number_of_particles

        self.sigma = np.sqrt(self.k*self.T/self.m)

        self.force_on_top = 0
        self.momentum_gaines = 0
        self.force_gained = 0

        self.numb_escaping = 0
        self.numb_coll = 0
        self.mass_escaped_per_sec = 0

        self.engine_stated = False

        self.v_escape = 12473





        self.x = np.zeros([int(self.N),3])
        self.v = np.zeros([int(self.N),3])




        for i in range(int(self.N)):
            for j in range(3):
                self.x[i,j] = -self.L/2+(random.random()*self.L)
                self.v[i,j] = random.gauss(0,self.sigma)#(3./2)*self.k*self.T)



    def integrate(self):


        self.x += self.v*self.dt


    def detect_collision_vec(self):


        for i in range(3):

            if i == 1:
                index_part_through_floor = np.logical_and(self.x[:,i] < -self.L/2, np.logical_and(np.logical_and(self.x[:,i+1] < self.hole/2, self.x[:,i+1] > -self.hole/2),np.logical_and(self.x[:,i-1] < self.hole/2,self.x[:,i-1] > -self.hole/2)))
                #index_part_through_floor = np.logical_and(self.x[:,i] < -self.L/2., np.logical_and(self.x[:,i-1] < self.hole/2.,self.x[:,i-1] > -self.hole/2.),np.logical_and(self.x[:,i+1] < self.hole/2.,self.x[:,i+1] > -self.hole/2.))
                self.momentum_gaines += np.abs(np.sum(self.v[index_part_through_floor,i])*self.m)
                self.numb_escaping += np.sum(index_part_through_floor)
                self.force_on_top += 2*np.sum(self.v[self.x[:,i] > self.L/2.,i])*self.m/self.dt
                self.numb_coll += np.sum(self.x[:,i] < -self.L/2.)





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

    def startEngine(self):

        for i in range(self.steps):

            self.detect_collision_vec()

            self.integrate()
            print (float(i)/self.steps)*100, "%            \r",
        print ""

        self.force_gained = self.momentum_gaines/(self.dt*self.steps)
        self.mass_escaped_per_sec = self.numb_escaping*self.m/self.interval

        self.engine_stated = True

        print "Engine has stated!"
        print "Data from the engine: "
        print "--------------------------------"
        print "Force: ", self.force_gained
        self.print_press()
        self.print_ke()
        print "Number of particles colliding with floor: ", self.numb_coll
        print "Number of particles escaped: ",self.numb_escaping
        print "Number of particles escaped per sec: ",self.numb_escaping/self.interval
        print "Mass lost per sec: " ,self.mass_escaped_per_sec
        print "--------------------------------"

    # def update_lines(self,num, dataLines, lines) :
    #     for line, data in zip(lines, dataLines) :
    #         line.set_data(data[0:2, num-1:num])
    #         line.set_3d_properties(data[2,num-1:num])
    #     return lines


    def plot_box(self,frames):
        def update_lines(num, dataLines, lines) :
            for line, data in zip(lines, dataLines) :
                line.set_data(data[0:2, num-1:num])
                line.set_3d_properties(data[2,num-1:num])
            return lines


        # Attach 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)

        m = frames   # number of frames in the animation
        n = self.N   # number of particles you want to animate
        N = self.steps # number of time steps in your data

        data = np.zeros([n, 3, N]) # this is the positions of the particles
        # to be animated. In this code it should be an array of shape (n, 3, N)

        # ...Fill out the data array here! data[p, :, i] is particle p, time step i.

        for i in range(N):

            for p in range(self.N):
                data[p,:,i] = self.x[p,:]
            self.detect_collision_vec()
            self.integrate()
            print (float(i)/self.steps)*100, "%            \r",
        print ""


        # creates animation data for all your different particles
        lines = [i for i in range(n)]
        for i in range(n):
            lines[i] = [ax.plot(data[i][0,0:1],
            data[i][1,0:1], data[i][2,0:1], 'o')[0]]

        # Set the axes properties
        ax.set_xlim3d([-self.L/2., self.L/2.])
        ax.set_xlabel('X')

        ax.set_ylim3d([-self.L/2., self.L/2.])
        ax.set_ylabel('Y')

        ax.set_zlim3d([-self.L/2., self.L/2.])
        ax.set_zlabel('Z')

        ax.set_title('3D random particles')

        # Creating the Animation object
        ani = [i for i in range(n)]
        for i in range(n):
            ani[i] = animation.FuncAnimation(fig,
            update_lines, m, fargs=([data[i]], lines[i]),
            interval=50, blit=False)
        plt.show()




    def launch(self, numberOfBoxes, burnTime, launch_step):

        launch_dt = float(burnTime)/launch_step

        total_mass_escaping_per_sec = numberOfBoxes*self.mass_escaped_per_sec
        mass_fuel = numberOfBoxes*self.N*self.m
        time = 0
        total_force = numberOfBoxes*self.force_gained
        fuel_mass = []
        fuel_mass.append(0)
        acceleration = 0

        reached_escape_vel = False


        launcher_vel = []
        launcher_vel.append(0)

        while (time < burnTime and launcher_vel[-1] < self.v_escape):

            acceleration = total_force/(self.launcher_mass)# + mass_fuel)

            launcher_vel.append(launcher_vel[-1] + acceleration*launch_dt)

            fuel_mass.append(fuel_mass[-1] + total_mass_escaping_per_sec*launch_dt)

            time += launch_dt

            print (float(time)/burnTime)*100, "%            \r",
        print ""

        if (launcher_vel[-1] > self.v_escape):
            reached_escape_vel = True

        self.checkFuelAmount(numberOfBoxes,fuel_mass[-1])

        if(reached_escape_vel):
            #print "Reached escape velocity with: ", fuel[-1], " to spare!"
            print "Fuel needed: ", fuel_mass[-1]
            print "After ", time, " seconds!"
        else:
            print "Did not reach escape velocity"
            print "Only reach ", launcher_vel[-1], " m/s"

        plt.plot([self.v_escape for i in range(len(launcher_vel))])
        plt.plot(launcher_vel)
        plt.plot(fuel_mass)

        plt.show()

        #return launcher_vel, fuel_mass, time, reached_escape_vel

    def checkFuelAmount(self, numberOfBoxes,fuel):

        seed = int("01595")
        myStarSystem = AST1100SolarSystem(seed)
        myStarSystem.massNeededCheck(numberOfBoxes,self.v_escape,self.force_gained, self.numb_escaping/interval, fuel)










interval = 1e-9
steps = 1000
frames = 1000
L = 1e-6
H = L/2
T = 10000
numb_part = 1.e5#25

e = engine(L,H,4130,T,interval,steps,numb_part)


#e.plot_box(frames)

e.startEngine()
e.launch(8.6552e12, 20*60,10000)
