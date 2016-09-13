import numpy as np
import matplotlib.pyplot as plt
from rocket import engine
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter


class satelite:
    def __init__(self,delta_v,x0,v0):
        interval = 1e-9
        steps = 1000
        frames = 1000
        L = 1e-6
        H = L/2
        T = 10000
        numb_part = 1e5



        self.G = 4*np.pi**2



        self.seed = 69558
        self.force_per_box = 0
        self.force_engine = 0
        self.mass_lost_per_box_per_sec = 0
        self.mass_lost_per_sec = 0
        self.spare_delta_v = 250
        self.fuel = 0
        self.delta_v = delta_v + self.spare_delta_v
        self.velocity = np.zeros(2)
        self.position = np.zeros(2)
        self.sim_time = 20
        self.steps_per_year = 20000
        self.dt = 1./self.steps_per_year

        self.writeingFreq = 10

        self.pos_over_time = np.zeros((2,self.sim_time*self.steps_per_year/self.writeingFreq))

        self.time_planet_simulated = 20
        # self.steps_per_year = 2000
        # self.dt = 1./(self.time_planet_simulated*self.steps_per_year)

        #self.time = np.linspace(0,self.time_planet_simulated,self.steps_per_year*self.steps_per_year/10)
        #print self.time.size
        self.numb_boxes = 1.333e13

        self.v0 = np.array(v0)
        self.x0 = np.array(x0)

        self.velocity = self.v0
        self.position = self.x0

        self.pos_over_time[:,0] = self.position

        self.e = engine(L,H,T,interval,steps,numb_part)
        self.system = AST1100SolarSystem(self.seed)

        self.numberOfPlanets = self.system.numberOfPlanets

        self.planetMasses = np.array(self.system.mass)
        self.starMass = self.system.starMass



        with open("positionsHomePlanet.npy", "rb") as npy:
            self.planetPositions = np.load(npy)

        # for p in range(self.numberOfPlanets):
        #     plt.plot(self.planetPositions[0,p,-1],self.planetPositions[1,p,-1],"o")
        #     plt.plot(self.planetPositions[0,p,0],self.planetPositions[1,p,0],"x")
        #     plt.plot(self.planetPositions[0,p,:],self.planetPositions[1,p,:])
        #
        #
        # plt.show()

        self.steps =  self.planetPositions[0,0,:].size

        self.time = np.linspace(0,self.time_planet_simulated,self.steps)

        self.planetPosFunction = inter.interp1d(self.time, self.planetPositions)





    def main_sequence(self):

        self.e.startEngine()
        self.force_per_box = self.e.force_gained
        self.mass_lost_per_box_per_sec = self.e.mass_escaped_per_sec
        self.force_engine = self.force_per_box*self.numb_boxes
        self.mass_lost_per_sec = self.mass_lost_per_box_per_sec*self.numb_boxes
        print "Engine ignited"
        self.fuel_cal(self.delta_v)
        print "Lauching with ", self.fuel, " kg fuel."
        print "---------------------"


    def fuel_cal(self,delta_v):
        self.fuel = self.e.launcher_mass*(np.exp(self.e.mass_escaped_per_sec*delta_v/self.e.force_gained) - 1)


    def boost(self,delta_v):
        pass


    def acceleration(self,planetPos):
        acceleration = np.zeros(2)
        r = planetPos.T - self.position

        acceleration = -np.sum((self.G*self.planetMasses)/(np.sum(r**2,axis = 1)**(3./2))*r.T,axis = 1) - (self.G*self.starMass/(np.sum((self.position**2))**(3./2)))*self.position.T
        return acceleration
    def main_loop(self):

        self.velocity += 0.5*self.acceleration(self.planetPositions[:,:,0])*self.dt



        for i in xrange(self.steps_per_year*self.sim_time-1):

            self.position += self.velocity*self.dt
            self.velocity += self.acceleration(self.planetPositions[:,:,i])*self.dt

            if ((i+1)%self.writeingFreq == 0 ):
                self.pos_over_time[:,i/self.writeingFreq+1] = self.position

        plt.plot(self.pos_over_time[0,:],self.pos_over_time[1,:])
        for p in range(self.numberOfPlanets):
            plt.plot(self.planetPositions[0,p,-1],self.planetPositions[1,p,-1],"o")
            plt.plot(self.planetPositions[0,p,0],self.planetPositions[1,p,0],"x")
            plt.plot(self.planetPositions[0,p,:],self.planetPositions[1,p,:])


        

        plt.axis("equal")
        plt.show()


x0 = [0.,2.7]
v0 = [.6*np.pi*1.7,0.]

sat = satelite(20000,x0,v0)
#sat.main_sequence()
sat.main_loop()
