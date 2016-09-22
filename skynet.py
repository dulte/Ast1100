import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn


class skynet:

    def __init__(self):



        with open("positionsHomePlanet.npy", "rb") as npy:
            self.planetPositions = np.load(npy)

        with open("himmelkule.npy") as him:
            self.himmelkule = np.load(him)

        self.time_planet_simulated = 20
        self.steps =  self.planetPositions[0,0,:].size
        self.time = np.linspace(0,self.time_planet_simulated,self.steps)
        self.planetPosFunction = inter.interp1d(self.time, self.planetPositions)

        self.seed = 75041
        self.system = AST1100SolarSystem(self.seed)
        self.numberOfPlanets = self.system.numberOfPlanets




    def find_distance_planets(self,sat_pos,time):
        distances = np.zeros(self.numberOfPlanets + 1)
        r = self.planetPosFunction(time) - sat_pos[:,np.newaxis]
        distances[:-1] = norm(r,axis = 0)
        distances[-1] = norm(sat_pos)
        r_with_sun = np.zeros([2,self.numberOfPlanets + 1])
        r_with_sun[:,:-1] = r
        r_with_sun[:,-1] = sat_pos
        return distances,r_with_sun


    def find_position(self, guess, weight,distances, real_r, time):

        x0 = guess[0]
        y0 = guess[1]



        regression_steps = 1e4

        for i in range(int(regression_steps)):
            dist,r = self.find_distance_planets(guess,time)

            dx = (4./dist.size)*sum((r[0] - real_r[0]))
            dy = (4./dist.size)*sum((r[1] - real_r[1]))

            guess[0] += weight*dx
            guess[1] += weight*dy

        return guess



    def test_dist(self,number_of_tests):

        epsilon = 1e-5
        correct = 0
        wrong = 0

        for i in range(number_of_tests):
            time = np.random.uniform(0,18)
            point = np.array([np.random.uniform(-50,50),np.random.uniform(-50,50)])

            guess =  np.array([np.random.uniform(-50,50),np.random.uniform(-50,50)])
            dist, r = self.find_distance_planets(point,time)
            guess = self.find_position(guess,0.01,dist,r,time)

            if np.max(abs(guess-point)) < epsilon:
                correct += 1
            else:
                wrong += 1


        print "Correct: ", correct
        print "Wrong: ",wrong


    def get_image(self,phi_midpoint):


        theta_midpoint = np.pi/2.
        ypix = 480
        xpix = 640
        field_of_view = 70.
        max_xy = (2*np.sin(field_of_view/2.))/(1+np.cos(field_of_view/2.))

        pic = np.zeros((ypix,xpix,3), dtype=np.uint8)
        x = np.linspace(-max_xy,max_xy,xpix)
        y = np.linspace(-max_xy,max_xy,ypix)

        print x.size


        for i in xrange(xpix):
            for j in xrange(ypix):
                rho = norm(np.array([x[i],y[j]]))
                c = 2*np.arctan2(rho,2.)


                theta = (np.pi/2) - np.arcsin(np.cos(c)*np.cos(theta_midpoint) + y[j]*np.sin(c)*np.sin(theta_midpoint)/rho)
                phi = phi_midpoint + np.arctan2(x[i]*np.sin(c),(rho*np.sin(theta_midpoint)*np.cos(c) - y[j]*np.cos(theta_midpoint)*np.sin(c)))

                temp = self.himmelkule[self.system.ang2pix(theta,phi)]

                rbg = np.array([temp[2],temp[3],temp[4]])

                pic[j,i,:] = rbg

        img = Image.fromarray(pic)

        img.save("star.png")



















skynet = skynet()
#dist, r =  skynet.find_distance_planets(np.array([-1.2,0.8]),2)
#skynet.find_position(np.array([0.0001,-0.0001]), 0.01,dist,r,2)
#skynet.test_dist(5)
skynet.get_image(43.5)
