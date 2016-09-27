import numpy as np
import matplotlib.pyplot as plt
from rocket import engine
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn


class satelite:
    def __init__(self,delta_v,destination):
        interval = 1e-9
        steps = 1000
        frames = 1000
        L = 1e-6
        H = L/2
        T = 10000
        numb_part = 1e5



        self.G = 4*np.pi**2
        self.km_to_au = 6.685e-9



        #self.seed = 69558
        self.seed = 75041

        self.e = engine(L,H,T,interval,steps,numb_part)
        self.system = AST1100SolarSystem(self.seed)

        self.force_per_box = 0
        self.force_engine = 0
        self.mass_lost_per_box_per_sec = 0
        self.mass_lost_per_sec = 0
        self.spare_delta_v = 250
        self.fuel = 0
        self.delta_v = delta_v + self.spare_delta_v
        self.velocity = np.zeros(2)
        self.position = np.zeros(2)
        self.sim_time = 5
        self.steps_per_year = 365*24*60*2#300000
        self.dt = 1./self.steps_per_year

        self.writeingFreq = 100

        self.pos_over_time = np.zeros ((2,self.sim_time*self.steps_per_year))#np.zeros((2,self.sim_time*self.steps_per_year/self.writeingFreq))

        self.time_planet_simulated = 20

        self.use_time0 = False

        self.numb_boxes = 1.333e13








        self.numberOfPlanets = self.system.numberOfPlanets
        self.destination_planet = destination
        self.home_planet = 0

        self.escape_velocity = self.calc_escape()



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

        self.calc_trad_parameters(self.destination_planet)
        self.angle_between_planets -=0.02
        self.calc_time_to_burn(self.angle_between_planets)
        #self.time_to_launch -= 0.4
        self.launch_angle = self.calc_tangental_angle(self.time_to_launch)


        if (self.use_time0):
            self.time_to_launch = 0

        print "e home: ", self.system.e[0]
        print "e destinaton: ", self.system.e[self.destination_planet]


        #v0 = [self.system.vx0[self.home_planet] + -np.cos(self.launch_angle)*(1.*self.calc_influence() +.956*self.dv_mainburn),0.*self.system.vy0[self.home_planet] + np.sin(self.launch_angle)*(1.*self.calc_influence() +.956*self.dv_mainburn)]#[0. + .2*self.escape_velocity,0*self.system.vy0[self.home_planet] + 0.3*self.escape_velocity]
        """0*0.767"""
        # self.position = np.array([self.system.x0[0],self.system.y0[0]])
        # v0 = [0,self.dv_mainburn]

        print "Angle: ",self.calc_tangental_angle(self.time_to_launch)
        print "Sin, ",np.sin(self.calc_tangental_angle(self.time_to_launch))
        print "Cos: ",np.cos(self.calc_tangental_angle(self.time_to_launch))

        self.position = np.array([self.planetPosFunction(self.time_to_launch)[0,0],self.planetPosFunction(self.time_to_launch)[1,0]]) +np.array([-np.sin(self.calc_tangental_angle(self.time_to_launch)),np.cos(self.calc_tangental_angle(self.time_to_launch))])*(self.system.radius[self.home_planet]*self.km_to_au)
        #v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn)]
        #For default tid: v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn + .8929*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn + .8929*self.calc_influence())]
        v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn + .9025*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn + .9025*self.calc_influence())]
        #v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.0025*self.dv_mainburn + 1.*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.0025*self.dv_mainburn + 1.*self.calc_influence())]
        #v0 = [0.,0.]
        print v0
        print "Total speed: ",norm(np.array(v0))

        # if (self.use_time0):
        #     v0 = [self.system.vx0[self.home_planet] ,0.*self.system.vy0[self.home_planet] + 0.*self.calc_influence() +1*self.dv_mainburn] #.956*self.dv_mainburn]#[0. + .2*self.escape_velocity,0*self.system.vy0[self.home_planet] + 0.3*self.escape_velocity]
        #     self.position = np.array([self.system.x0[0] +0*self.system.radius[self.home_planet]*self.km_to_au ,self.system.y0[0] + 1*self.system.radius[self.home_planet]*self.km_to_au])
        # else:
        #     v0 = [0*self.calc_vel(self.home_planet, self.time_to_launch)[0]+ -np.sin(self.launch_angle)*(0.*self.calc_influence() +1.048*self.dv_mainburn),0.*self.system.vy0[self.home_planet] + np.cos(self.launch_angle)*(0.*self.calc_influence() +1.048*self.dv_mainburn)]#[0. + .2*self.escape_velocity,0*self.system.vy0[self.home_planet] + 0.3*self.escape_velocity]
        #     self.position = self.planetPosFunction(self.time_to_launch)[:,0] + np.array([-np.sin(self.launch_angle)*(self.system.radius[self.home_planet]*self.km_to_au),np.cos(self.launch_angle)*(self.system.radius[self.home_planet]*self.km_to_au)])
        #
        #print norm(np.array(v0))




        self.velocity = np.array(v0)
        self.pos_over_time[:,0] = self.position





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
        # acceleration = np.zeros(2)
        # r = planetPos.T - self.position
        #
        # to_low = np.abs(r) < self.system.radius[self.home_planet]*self.km_to_au
        # #r[to_low] = self.system.radius[self.home_planet]*self.km_to_au
        #
        # acceleration = np.sum((self.G*self.planetMasses)/(norm(r,axis = 1)**3)*r.T,axis = 1) #-(self.G*self.starMass/(norm(self.position,axis = 0)**3))*self.position
        # return acceleration
        # return -(self.G*self.starMass/(norm(self.position,axis = 0)**3))*self.position


        # acceleration = np.zeros(2)
        #
        # rx = (planetPos[0,:self.destination_planet +1] - self.position[0])
        # ry = (planetPos[1,:self.destination_planet +1] - self.position[1])

        # r = (planetPos[:,:self.destination_planet + 1] - self.position)
        #
        #
        # return np.sum((self.G*self.planetMasses[:self.destination_planet +1]/(norm(r)**3))*r) - (self.G*self.starMass/(norm(self.position,axis = 0)**3))*self.position
        # some = norm(np.array([rx,ry]),axis = 0)
        # acceleration[0] = np.sum((self.G*self.planetMasses[:self.destination_planet +1]/(some**3)*rx))
        # acceleration[1] = np.sum((self.G*self.planetMasses[:self.destination_planet +1]/(some**3)*ry))
        # return acceleration - (self.G*self.starMass/(norm(self.position,axis = 0)**3))*self.position
        #return -(self.G*self.starMass/(norm(self.position,axis = 0)**3))*self.position


        #acceleration = np.zeros(2)

        r = planetPos[:,:2] - self.position[:,np.newaxis]

        #r = planetPos[:,i] - self.position
        #acceleration = np.sum(self.G*self.planetMasses/(norm(r,axis = 0)**3)*r,axis = 1)
        #print acceleration
        return np.sum(self.G*self.planetMasses[:2]/(norm(r,axis = 0)**3)*r,axis = 1) - (self.G*self.starMass/(norm(self.position,axis = 0)**3))*self.position

    def main_loop(self):


        print "---------------------------------"
        print "Beginning to calculate the orbit: "

        if (self.use_time0):
            time = [0.]
        else:
            time =[self.time_to_launch]

        time_save = []

        self.velocity += 0.5*self.acceleration(self.planetPosFunction(time[-1]))*self.dt
        #print self.planetPosFunction(time[-1])

        min_dist_to_planet = 1000
        time_closest_incounter = 0

        dt = self.dt

        start = tid.clock()

        couter = 0
        r_dest = 100

        for i in xrange(1,self.steps_per_year*self.sim_time):
        #while (time[-1] < self.steps_per_year*self.sim_time):


            r_dest = abs(self.calc_dist_to_planet(self.destination_planet,time[-1],self.position))
            r_home = abs(self.calc_dist_to_planet(self.home_planet,time[-1],self.position))

            if (r_dest < 0.005 or r_home < 0.005):
                dt = self.dt/60
            else:
                dt = self.dt



            self.position += self.velocity*dt

            time.append(time[-1] + dt)
            self.velocity += self.acceleration(self.planetPosFunction(time[-1]))*dt

            self.pos_over_time[:,i] = self.position
            time_save.append(time[-1])

            #self.velocity += self.acceleration(self.planetPosFunction(0))*self.dt





            if (r_dest < min_dist_to_planet):
                min_dist_to_planet = r_dest# abs(self.calc_dist_to_planet(self.destination_planet,time[-1],self.position))
                time_closest_incounter = time[-1]


            couter += 1

            # if ((i+1)%self.writeingFreq == 0 ):
            #     self.pos_over_time[:,i/self.writeingFreq+1] = self.position
            #     time_save.append(time[-1])



            #print (float(i)/(self.steps_per_year*self.sim_time))*100, "%            \r",

        print ""

        print "It took", (tid.clock()-start), " sec"


        print "Shortest distance to planet: ", min_dist_to_planet
        print "At time: ",time_closest_incounter
        plt.plot(self.pos_over_time[0,0],self.pos_over_time[1,0],"x")
        plt.plot(self.pos_over_time[0,:],self.pos_over_time[1,:])
        plt.plot(self.pos_over_time[0,-1],self.pos_over_time[1,-1],"o")

        plt.plot(self.planetPosFunction(self.time_to_launch)[0,0],self.planetPosFunction(self.time_to_launch)[1,0],"rv")
        plt.plot(self.planetPosFunction(self.time_to_launch)[0,0],self.planetPosFunction(self.time_to_launch)[1,0],"gv")
        plt.plot(self.planetPosFunction(self.time_to_launch)[0,0] - np.sin(self.calc_tangental_angle(self.time_to_launch))*0.01,self.planetPosFunction(self.time_to_launch)[1,0] + np.cos(self.calc_tangental_angle(self.time_to_launch))*0.01,"yv")

        pos_planets =  self.planetPosFunction(np.array(time_save))

        for p in xrange(self.numberOfPlanets):
            plt.plot(pos_planets[0,p,:],pos_planets[1,p,:])


        plt.plot(self.planetPosFunction(time_closest_incounter)[0,self.destination_planet],self.planetPosFunction(time_closest_incounter)[1,self.destination_planet],"*")





        plt.axis("equal")
        plt.show()


    def calc_escape(self):

        return np.sqrt(2*self.G*self.system.mass[self.home_planet]/(self.system.radius[self.home_planet]*self.km_to_au))

    def calc_influence(self):
        k = 10
        r = norm(np.array([self.system.x0[0],self.system.y0[0]]))
        r_soi = r*(self.planetMasses[0]/self.starMass)**(2./5)

        r_inf = r/np.sqrt(k*self.starMass/self.planetMasses[0])
        return np.sqrt(-(2*self.G*self.planetMasses[0])/(r_soi) + 2*self.G*self.planetMasses[0]/(self.system.radius[0]*self.km_to_au))

    def calc_trad_parameters(self,planet_number):

        time_mod = -0.

        r1 = norm(np.array([self.system.x0[0],self.system.y0[0]]))
        r2 = norm(np.array([self.system.x0[planet_number],self.system.y0[planet_number]]))


        self.dv_mainburn = np.sqrt(2*(self.G*self.starMass/r1 - self.G*self.starMass/(r1+r2)))

        self.time_to_encounter = time_mod + np.pi*np.sqrt((r1+r2)**3/(8*self.G*self.starMass))

        self.omega_destination = 2*np.pi/(np.sqrt((4*np.pi**2 * self.system.a[planet_number]**3)/(self.G*(self.starMass+self.planetMasses[planet_number]))))
        self.omega_home = 2*np.pi/(np.sqrt((4*np.pi**2 * self.system.a[self.home_planet]**3)/(self.G*(self.starMass+self.planetMasses[self.home_planet]))))
        self.angle_between_planets = np.pi - self.omega_destination*self.time_to_encounter

        print "Burn required: ",self.dv_mainburn
        print "Trip takes: ", self.time_to_encounter
        print "Planet moves at angular v: ", self.omega_destination
        print "Home planets anguar v: ", self.omega_home
        print "Angle between planet when burn: ",self.angle_between_planets

    def calc_time_to_burn(self,angle):

        print "Beginning to calculate burn time. This may take a couple of minutes..."

        min_angle = 10000.
        time_for_min_angle = 0
        #time = 0.0
        time = 6.
        eps = 1.e-4

        sim_time = 20.
        steps = 20000
        dt = 1./steps
        for i in xrange(1,int(steps*sim_time)-1):
            angle_home = np.arctan2(self.planetPosFunction(time)[1,self.home_planet],self.planetPosFunction(time)[0,self.home_planet])
            angle_destination = np.arctan2(self.planetPosFunction(time)[1,self.destination_planet],self.planetPosFunction(time)[0,self.destination_planet])

            time += dt

            if abs((angle_destination - angle_home)-angle) < eps:
                break

            print (float(i)/(steps*sim_time))*100, "%            \r",
        print ""


        self.time_to_launch = time
        print "Smallest difference in angle: ",abs((angle_destination - angle_home)-angle)
        print "At time: ",self.time_to_launch


    def calc_dist_to_planet(self,planet_number,time,pos):

        r = norm(self.planetPosFunction(time)[:,planet_number] - pos)
        return abs(r)

    def calc_tangental_angle(self,time):
        return  np.arctan2(self.planetPosFunction(time)[1,self.home_planet],self.planetPosFunction(time)[0,self.home_planet])


    def calc_vel(self,planet_number,time):
        dx = (self.planetPosFunction(time+self.dt)[0,planet_number] - self.planetPosFunction(time-self.dt)[0,planet_number])/(2*self.dt)
        dy = (self.planetPosFunction(time+self.dt)[1,planet_number] - self.planetPosFunction(time-self.dt)[1,planet_number])/(2*self.dt)

        return np.array([dx,dy])

    def calc_injection_burn(self,time):
        r =   self.position - self.planetPosFunction(time)[:,self.destination_planet]

        theta = np.arctan2(r[1],r[0])

        orb_vel = np.sqrt(self.G*self.planetMasses[self.destination_planet]/norm(r))

        return np.array([-orb_vel*np.sin(theta),orb_vel*np.cos(theta)])-self.velocity

    def calc_fuel(self,dv):
        n_e = 2.120316e-13
        f_b = 1.72650965887e-09

        return 1100*(np.exp(dv*n_e/f_b)-1)





destinaton_planet = 1
sat = satelite(20000,destinaton_planet)
#sat.main_sequence()
print sat.calc_escape()
sat.main_loop()
