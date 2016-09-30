import numpy as np
import matplotlib.pyplot as plt
from rocket import engine
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn
import sys


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
        self.sim_time = 2.15
        self.steps_per_year = 300000#365*24*60#300000
        self.dt = 1./self.steps_per_year

        """
        e6 | d
        1    .003
        1.2| .05
        1.5| .12
        2.0| .183
        3.0| .228
        4.0| .2486
        """
        """
        e6 | d
        .4 | 0.4547
        .6 | 0.4292
        1.0| 0.3957
        1.5| 0.3714
        2.0| 0.3564
        """


        self.writeingFreq = 100

        self.pos_over_time = np.zeros ((2,int(self.sim_time*self.steps_per_year/self.writeingFreq)))#np.zeros((2,self.sim_time*self.steps_per_year/self.writeingFreq))

        self.time_planet_simulated = 20

        self.use_time0 = False

        self.numb_boxes = 1.333e13

        self.save = False








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
        #print "Sin, ",np.sin(self.calc_tangental_angle(self.time_to_launch))
        #print "Cos: ",np.cos(self.calc_tangental_angle(self.time_to_launch))
        self.optimal_dist = self.calc_optimal_dist()
        print "Optimal radius of orbit: ",self.optimal_dist


        self.position = np.array([self.planetPosFunction(self.time_to_launch)[0,0],self.planetPosFunction(self.time_to_launch)[1,0]]) +np.array([-np.sin(self.calc_tangental_angle(self.time_to_launch)),np.cos(self.calc_tangental_angle(self.time_to_launch))])*(self.system.radius[self.home_planet]*self.km_to_au)
        #v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn)]
        #For default tid: v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn + .8929*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn + .8929*self.calc_influence())]
        #v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn + .9025*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.004*self.dv_mainburn + .9025*self.calc_influence())]
        #v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.0025*self.dv_mainburn + 1.*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.0025*self.dv_mainburn + 1.*self.calc_influence())]
        #v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.0*self.dv_mainburn + .89622*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.0*self.dv_mainburn + .89622*self.calc_influence())]


        #real v0: (if does not work, check the angle(-=0.2) and time (2.15))
        v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.0*self.dv_mainburn + .894*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.0*self.dv_mainburn + .894*self.calc_influence())]



        #v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(np.sqrt(self.dv_mainburn**2 + 1.0031*self.calc_influence()**2)),np.cos(self.calc_tangental_angle(self.time_to_launch))*(np.sqrt(self.dv_mainburn**2 + 1.0031*self.calc_influence()**2))]


        self.v0_save = v0
        self.planet_vel = self.calc_vel(0,self.time_to_launch)
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




        self.velocity = np.array(v0) + self.planet_vel
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


        r = planetPos[:,:] - self.position[:,np.newaxis]

        return np.sum(self.G*self.planetMasses[:]/(norm(r,axis = 0)**3)*r,axis = 1) - (self.G*self.starMass/(norm(self.position,axis = 0)**3))*self.position

    def main_loop(self):


        print "---------------------------------"
        print "Beginning to calculate the orbit: "

        #time = np.zeros(self.steps_per_year*self.sim_time)
        if (self.use_time0):
            time = 0.
        else:
            time =self.time_to_launch

        time_save = [time]

        #print self.planetPosFunction(time[-1])

        r_dest = (self.calc_dist_to_planet(self.destination_planet,time,self.position))
        r_home = (self.calc_dist_to_planet(self.home_planet,time,self.position))

        min_dist_to_planet = 200
        time_closest_incounter = 0
        index_closest_incounter = 0

        dt_mod = 1000
        dt = self.dt/dt_mod

        start = tid.clock()

        couter = 0
        r_dest = 100

        lauched = False
        close_to_dest = False
        inject_burn = False
        correct_burn = False
        check_feq = 1000



        self.velocity -= 0.5*self.acceleration(self.planetPosFunction(time))*self.dt/dt_mod

        for i in xrange(1,int(self.steps_per_year*self.sim_time)):
        #while (time[-1] < self.steps_per_year*self.sim_time):

            if (i%check_feq == 0):
                r_dest = (self.calc_dist_to_planet(self.destination_planet,time,self.position))
                r_home = (self.calc_dist_to_planet(self.home_planet,time,self.position))

                if r_home > 0.005 and not lauched:

                    print "You have lauched"
                    lauched = True
                elif r_dest < 0.0005 and not close_to_dest:
                    print "Getting close to destinaton"
                    dt_mod = 1000
                    close_to_dest = True
                elif r_dest > 0.0005 and close_to_dest:
                    print "Getting away from the planet"
                    close_to_dest = False
                elif r_dest < 0.003 and not correct_burn:
                    print "Doing a corretion burn at time ",time
                    self.correction_burn = self.calc_correction_burn(time,.5)#.61)
                    print "With velocity ", self.correction_burn
                    self.velocity += self.correction_burn
                    correct_burn = True
                    check_feq = 10
                elif r_dest < self.optimal_dist and not inject_burn:
                    self.injection_burn = self.calc_injection_burn(time)
                    self.velocity += self.injection_burn
                    print "Burning for orbit, at time ", time
                    print "With velocity ", self.injection_burn
                    inject_burn = True


            if (close_to_dest or not lauched):
                dt = self.dt/dt_mod
                for j in xrange(dt_mod):
                    self.velocity += self.acceleration(self.planetPosFunction(time + j*dt))*dt
                    self.position += self.velocity*dt

            else:
                dt = self.dt
                self.velocity += self.acceleration(self.planetPosFunction(time + dt))*dt
                self.position += self.velocity*dt


            time += self.dt
            #self.pos_over_time[:,i] = self.position


            #self.velocity += self.acceleration(self.planetPosFunction(0))*self.dt





            if (r_dest < min_dist_to_planet):
                min_dist_to_planet = r_dest# abs(self.calc_dist_to_planet(self.destination_planet,time[-1],self.position))
                time_closest_incounter = time
                index_closest_incounter = float(i)/self.writeingFreq


            #couter += 1

            if ((i)%self.writeingFreq == 0 ):
                self.pos_over_time[:,i/self.writeingFreq] = self.position
                time_save.append(time)
            if ((i)%10000 == 0 ):
                if r_home > r_dest:
                    print r_dest
                else:
                    print r_home



            #print (float(i)/(self.steps_per_year*self.sim_time))*100, "%            \r",




        print ""


        print "It took", (tid.clock()-start), " sec"
        #time_save = time



        print "Shortest distance to planet: ", min_dist_to_planet
        print "At time: ",time_closest_incounter
        print "Planet position is: ", self.planetPosFunction(time_closest_incounter)[:,self.destination_planet]
        print "Satelite position is: ",self.pos_over_time[:,int(index_closest_incounter)]

        print "-----------------------"
        print "dv and Fuel"
        print "Hohmann burn: ", norm(self.dv_mainburn)
        print "Main Burn: ", norm(self.v0_save)
        print "Correction Burn: ",norm(self.correction_burn)
        print "Injection Burn: ", norm(self.injection_burn)
        print "Total: ", norm(self.v0_save) + norm(self.correction_burn) +norm(self.injection_burn)
        print "Total fuel: ", self.calc_fuel(norm(self.v0_save)+norm(self.correction_burn)+norm(self.injection_burn))

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



        r_relative_from_dest = self.planetPosFunction(np.array(time_save))[:,self.destination_planet] - self.pos_over_time




        plt.axis("equal")
        plt.show()

        plt.plot(r_relative_from_dest[0],r_relative_from_dest[1])
        plt.plot(0,0,"ro")
        plt.axis("equal")
        plt.show()

        if self.save:
            np.save("posOverTime.npy",self.pos_over_time)
            np.save("time.npy",time_save)


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


        self.v_mainburn = np.sqrt(2*(self.G*self.starMass/r1 - self.G*self.starMass/(r1+r2)))
        self.dv_mainburn = np.sqrt(self.G*self.starMass/r1)*(np.sqrt(2*r2/(r1+r2)) - 1)
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

        return np.array([-orb_vel*np.sin(theta),orb_vel*np.cos(theta)])-self.velocity+self.calc_vel(self.destination_planet,time)

    def calc_fuel(self,dv):
        n_e = 2.120316e-13
        f_b = 1.72650965887e-09

        return 1100*(np.exp((dv*4744)*n_e/f_b)-1)

    def calc_optimal_dist(self):
        return self.system.a[self.destination_planet]*np.sqrt(self.planetMasses[self.destination_planet]/(10*self.starMass))

    def calc_correction_burn(self,time,factor):
        r = self.position - self.planetPosFunction(time+0.001)[:,self.destination_planet]
        normal_r = 1.0/(norm(r)) * r
        print "Length normal r: ", norm(normal_r)
        return -factor*normal_r



destinaton_planet = 1
sat = satelite(20000,destinaton_planet)
#sat.main_sequence()
print sat.calc_escape()
sat.main_loop()
