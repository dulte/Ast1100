import numpy as np
#from PIL import Image
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn

class orbit:
    def __init__(self,seed,filename,pos,vel):

        self.x0 = pos
        self.v0 = vel
        sys = AST1100SolarSystem(seed)
        self.sys=sys
        planetID = 1
        self.filename = filename
        self.mu = 28.019072
        self.solar_to_kg = 1.988435e30
        self.km_to_au = 6.685e-9
        self.mh = 1.660539040e-27
        self.k = 1.38064852e-23
        self.G = 6.674e-11
        self.g = self.G*sys.mass[planetID]*self.solar_to_kg/(sys.radius[planetID]*1000)**2
        self.laucherMass = 1100.0

        self.gamma = 1.4

        self.tempStar = sys.temperature
        self.tempPlanet = np.sqrt(sys.starRadius*self.km_to_au/(2.0*sys.a[planetID]))*self.tempStar
        self.rho = sys.rho0[planetID]
        self.laucherArea = 1



        self.mass = sys.mass[planetID]*self.solar_to_kg

        print self.tempPlanet - 273.15,self.g



    def isothermPress(self,h):

        h0 = self.k*self.tempPlanet/(2*self.mh*self.mu*self.g)
        rho1 = self.rho*(0.5)**(1.0/(self.gamma+1))
        hightShift = (self.gamma)/(2*(self.gamma -1))*h0*2

        return rho1*np.exp(-(h-hightShift)/h0)

    def adiabaticPress(self,h):
        h0 = self.k*self.tempPlanet/(self.mh*self.mu*self.g)
        gamma = 1.4

        return self.rho*(1-(gamma-1)/gamma *h/h0)**(1.0/(gamma+1))

    def adiabaticTemp(self,h):
        h0 = self.k*self.tempPlanet/(self.mh*self.mu*self.g)
        return self.tempPlanet*(1-(self.gamma -1)/self.gamma *h/h0)

    def atmosPress(self,h):
        temp = self.adiabaticTemp(h)
        return np.where(temp < self.tempPlanet/2.0,self.isothermPress(h),self.adiabaticPress(h))

    def dragForce(self,rho,A,v):
        return .5*rho*A*v**2

    def gravityForce(self,r):
        return self.G*self.mass*self.laucherMass/(r**2)

    def findSafeHeight(self):
        for h in range(80000,400000,1):
            v = np.sqrt(self.G*self.mass/h)
            fd = self.dragForce(self.atmosPress(h),self.laucherArea,v)
            fg = self.gravityForce(h)
            if fg/fd > 1000:
                return h

    def hohmann1(self,r1,r2,v):
        gravPara = self.mass*self.G
        burn1 = np.sqrt(gravPara/r1)*(-1+np.sqrt(2*r2/(r1+r2)))
        time = np.pi*np.sqrt((r1+r2)**3 /(8*gravPara))
        v_norm = v/norm(v)
        return burn1*v_norm,time

    def circle(self,r1,v):
        gravPara = self.mass*self.G
        #burn2 = np.sqrt(gravPara/r2)*(1-np.sqrt(2*r1/(r1+r2)))
        burn2 = self.orbVel(norm(r1))
        theta = np.arctan2(r1[1],r1[0])
        v_norm = np.array([-np.sin(theta),np.cos(theta),0])
        return burn2*v_norm - v

    def hohmann2(self,r1,r2,v):
        gravPara = self.mass*self.G
        r = norm(r1)
        print norm(r1)
        burn2 = np.sqrt(gravPara/r2)*(1-np.sqrt(2*r/(r+r2)))
        #burn2 = self.orbVel(norm(r1))
        theta = np.arctan2(r1[1],r1[0])
        v_norm = v/norm(v)#np.array([-np.sin(theta),np.cos(theta),0])
        print norm(burn2)
        print self.orbVel(norm(r1))
        return burn2*v_norm - v

    def lauch(self):
        self.sys.landOnPlanet(1,filename)

    def orbVel(self,r):
        return np.sqrt(self.G*self.mass/r)



    def a(self,x):
        return -self.G*self.mass/(norm(x)**3)*x

    def simOrbit(self,time):
        dt = 1e-2

        pos_temp = np.zeros(3)
        pos_temp = np.copy(self.x0)
        writingFreq = 10000.0

        pos = np.zeros((int(time/(writingFreq*dt)),3))
        pos[0,:] = np.copy(self.x0)

        t = np.zeros(int(time/(dt)))

        vel = np.zeros(3)
        vel = np.copy(self.v0)
        vel += self.circle(pos_temp,vel)


        dv1,timeBurn = self.hohmann1(norm(self.x0),3000000,self.v0)

        hohmann2Burned = False

        vel += dv1

        print norm(self.x0)

        vel += 0.5*self.a(pos[0,:])*dt

        for i in range(1,int(time/dt)):

            pos_temp += vel*dt
            vel += self.a(pos_temp)*dt
            t[i] = t[i-1] + dt

            if i*dt >= timeBurn and not hohmann2Burned:
                print "Sirukerer"
                vel -= 0.5*self.a(pos[0,:])*dt
                vel += self.circle(pos_temp,vel)

                #vel += self.hohmann2(self.x0,3000000,vel)
                #vel = self.orbVel(norm(pos_temp))*vel/norm(vel)
                hohmann2Burned = True
                vel += 0.5*self.a(pos[0,:])*dt


            if i%writingFreq == 0:
                pos[int(i/writingFreq),:] = pos_temp

                print (float(i)/int(time/(dt)))*100, "%            \r",
        print ""

        # fig = plt.figure()
        # ax = fig.gca(projection = "3d")
        # ax.plot(pos[:,0],pos[:,1],pos[:,2])
        # plt.title("Orbit")
        # plt.xlabel("x")
        # plt.ylabel("y")
        # plt.show()
        plt.plot(pos[:,0],pos[:,1])
        plt.axis("equal")
        plt.show()




seed = 75041

startPos = (np.array([ 50046912.943 ,  1922577.0044 ,  0]))
startV = np.array([-43.5115671208 ,  667.405339613 ,  0])

filename = "landing.txt"

orb = orbit(seed,filename,startPos,startV)

h = np.linspace(0,200000,1000001)
# print "Safe hight: ",orb.findSafeHeight()
# plt.plot(h,orb.adiabaticPress(h))
# plt.show()
# plt.plot(h,orb.atmosPress(h))
# plt.show()
orb.simOrbit(100000)
#orb.lauch()
