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
        return distances


    def find_position(self, guess, weight,distances, time):

        x0 = guess[0]
        y0 = guess[1]

        h = 1e-8

        regression_steps = 1e5

        for i in range(int(regression_steps)):
            dist= self.find_distance_planets(guess,time)


            dx = (8./dist.size)*np.sum(dist[:-1]*(distances[:-1]-dist[:-1])*(self.planetPosFunction(time)[0,:]-guess[0]))
            dy = (8./dist.size)*np.sum(dist[:-1]*(distances[:-1]-dist[:-1])*(self.planetPosFunction(time)[1,:]-guess[1]))


            guess[0] -= weight*dx
            guess[1] -= weight*dy

        return guess





    def test_dist(self,number_of_tests):

        epsilon = 1e-5
        correct = 0
        wrong = 0

        for i in range(number_of_tests):
            time = np.random.uniform(0,18)
            point = np.array([np.random.uniform(-50,50),np.random.uniform(-50,50)])

            guess =  np.array([np.random.uniform(-50,50),np.random.uniform(-50,50)])
            dist = self.find_distance_planets(point,time)
            guess = self.find_position(guess,0.00001,dist,time)

            if np.max(abs(guess-point)) < epsilon:
                correct += 1
            else:
                wrong += 1


        print "Correct: ", correct
        print "Wrong: ",wrong


    def make_image(self,phi_midpoint, name = ""):


        folder = "pic/"
        phi_midpoint = phi_midpoint*np.pi/180
        theta_midpoint = np.pi/2.
        ypix = 480
        xpix = 640
        field_of_view = 1.22173#70.
        max_xy = (2*np.sin(field_of_view/2.))/(1+np.cos(field_of_view/2.))

        pic = np.zeros((ypix,xpix,3), dtype=np.uint8)
        x = np.linspace(-max_xy,max_xy,xpix)
        y = np.linspace(-max_xy,max_xy,ypix)


        for i in xrange(xpix):
            for j in xrange(ypix):
                rho = norm(np.array([x[i],y[j]]))
                c = 2*np.arctan2(rho,2.)


                theta = (np.pi/2.) - np.arcsin(np.cos(c)*np.cos(theta_midpoint) + y[j]*np.sin(c)*np.sin(theta_midpoint)/rho)
                phi = phi_midpoint + np.arctan2(x[i]*np.sin(c),(rho*np.sin(theta_midpoint)*np.cos(c) - y[j]*np.cos(theta_midpoint)*np.sin(c)))

                temp = self.himmelkule[self.system.ang2pix(theta,phi)]

                rbg = np.array([temp[2],temp[3],temp[4]])

                pic[j,i,:] = rbg



        if name != "":
            img = Image.fromarray(pic)
            img.save(folder + name + ".png")
            return

        return pic


    def make_sky(self):
        phi = 0
        ypix = 480
        xpix = 640

        degrees = 360

        sky = np.zeros((degrees,ypix,xpix,3), dtype=np.uint8)

        for i in range(degrees):
            #name = "degree" + str(i)
            sky[i,:,:,:] = self.make_image(i)

            print i


        np.save("sky.npy",sky)


    def find_angle(self,pic):
        with open("sky.npy", "rb") as infile:
            sky = np.load(infile)


        degrees = 360

        least_error = 1e7
        least_error_deg = 0

        for i in range(degrees):

            pic_from_file = sky[i,:,:,:]

            error = np.sum((pic_from_file-pic)**2)

            if (error < least_error):
                least_error = error
                least_error_deg = i


        return least_error_deg



    def test_find_angle(self,number_of_tests = 1):

        correct = 0
        wrong = 0

        eps = 1

        for i in range(number_of_tests):
            test_deg = (np.random.uniform(0,359))

            pic = self.make_image(test_deg)

            calculated_deg = self.find_angle(pic)

            if (abs(calculated_deg - test_deg) < eps):
                correct += 1
            else:
                wrong += 1


        print "Correct: ",correct
        print "Wrong: ",wrong


    def find_velocity(self,lambda1,lambda2):
        c = 63239.7263
        h_alpha = 656.3

        star1_phi = 77.518724
        star1_shift = -0.017002884383
        star1_v = (star1_shift/h_alpha)*c

        star2_phi = 325.121916
        star2_shift = 0.017144686369
        star2_v = (star2_shift/h_alpha)*c


        v1 = (star1_v - (lambda1/h_alpha)*c)
        v2 = (star2_v - (lambda2/h_alpha)*c)

        div_factor = 1.0/np.sin(star2_phi-star1_phi)

        vx = div_factor*(np.sin(star2_phi)*v1 - np.sin(star1_phi)*v2)
        vy = div_factor*(-np.cos(star2_phi)*v1 + np.cos(star1_phi)*v2)

        return vx,vy


    def test_find_velocity(self):


        phi1 = 77.518724
        phi2 = 325.121916
        vs_1 = -1.63836317948
        vs_2 = 1.65202692896
        eps = 1e-6




        vtx,vty = self.find_velocity(-0.017002884383, 0.017144686369)


        vx,vy = self.find_velocity(0.0, 0.0)
        v1 = np.cos(phi1)*vx + np.sin(phi1)*vy
        v2 = np.cos(phi2)*vx + np.sin(phi2)*vy


        success = (abs(vtx - 0.0) < eps and abs(vtx - 0.0)<eps) and (abs(v1 - vs_1) < eps and abs(v2 - vs_2) <eps)

        assert success, "There is something wrong with the velocity finder"



    def give_orientation(self,lambda1, lambda2, picname,dist,time):

        pos = self.find_position([1,1],0.00001,dist,time)
        vel = self.find_velocity(lambda1,lambda2)
        pic = np.load(picname)
        angle = self.find_angle(pic)


        print "Position: ",pos
        print "Velocity: ",vel
        print "Angle: ",angle











skynet = skynet()
#dist =  skynet.find_distance_planets(np.array([1.2,100.8]),2)
#print skynet.find_position(np.array([-1.0001,1.0001]), 0.000001,dist,2)
skynet.test_dist(10)
#skynet.test_find_angle(5)
#skynet.test_find_velocity()

#skynet.get_image(43.5)
