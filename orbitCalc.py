import numpy as np
#from PIL import Image
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
from mpl_toolkits.mplot3d import Axes3D
import seaborn

class orbit:
    def __init__(self,seed,filename,pos,vel):

        self.x0 = pos
        self.v0 = vel
        sys = AST1100SolarSystem(seed)
        self.sys=sys
        planetID = 1
        self.filename = filename
        self.mu = 31.01404
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
        self.laucherArea = 6



        self.mass = sys.mass[planetID]*self.solar_to_kg



        print "temp: ",self.tempPlanet
        print "mass: ",self.mass
        print "rho0: ",self.rho
        print "radius: ",self.sys.radius[planetID]



    def isothermDensity(self,h):

        h0 = self.k*self.tempPlanet/(2*self.mh*self.mu*self.g)
        rho1 = self.rho*(0.5)**(1.0/(self.gamma-1))
        hightShift = (self.gamma)/(2*(self.gamma-1))*h0*2

        return rho1*np.exp(-(h-hightShift)/h0)

    def adiabaticDensity(self,h):
        h0 = self.k*self.tempPlanet/(self.mh*self.mu*self.g)
        gamma = 1.4

        return self.rho*(1-(gamma-1)/gamma *h/h0)**(1.0/(gamma-1))

    def adiabaticTemp(self,h):
        h0 = self.k*self.tempPlanet/(self.mh*self.mu*self.g)
        return self.tempPlanet*(1-(self.gamma -1)/self.gamma *h/h0)

    def atmosDensity(self,h):
        temp = self.adiabaticTemp(h)
        return np.where(temp < self.tempPlanet/2.0,self.isothermDensity(h),self.adiabaticDensity(h))

    def dragForce(self,rho,A,v):
        return .5*rho*A*v**2

    def gravityForce(self,r):
        return self.G*self.mass*self.laucherMass/(r**2)

    def findSafeHeight(self):
        for h in range(80000,8000000,10):
            if self.atmosDensity(h) < 1e-12:
                print "Safe Hight in func[km]: ",h/1000.
                return h
            # v = np.sqrt(self.G*self.mass/h)
            # fd = self.dragForce(self.atmosDensity(h),self.laucherArea,v)
            # fg = self.gravityForce(h)
            # if fg/fd > 1000:
            #     print "Force due to gravity: ",fg
            #     print "Force due to drag: ", fd
            #     return h

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

    def lauch(self):
        self.sys.landOnPlanet(1,filename)

    def orbVel(self,r):
        return np.sqrt(self.G*self.mass/r)



    def a(self,x,v):
        return -self.G*self.mass/(norm(x)**3)*x - 0.5*self.atmosDensity(norm(x)-self.sys.radius[1]*1000)*norm(v)*v*self.laucherArea/1100.

    def simOrbit(self,time,finalR,startTime):
        dt = 1e-1

        pos_temp = np.zeros(3)
        pos_temp = self.x0
        writingFreq = 1000.0

        pos = np.zeros((int(time/(writingFreq*dt)),3))
        pos[0,:] = self.x0

        #t = np.zeros(int(time/(writingFreq*dt)))

        vel = np.zeros(3)
        vel = self.v0
        firstCirc = self.circle(pos_temp,vel)
        vel += firstCirc#self.circle(pos_temp,vel)

        hight = finalR #4000000
        dv1,timeBurn = self.hohmann1(norm(self.x0),hight,self.v0)

        totalFirstBoost = dv1 + firstCirc
        print "Safe Hight: ",hight
        print "First hohmann: ", totalFirstBoost
        print "Sirkulering etter %g sek" %(timeBurn + startTime)
        hohmann2Burned = False

        vel += dv1

        vel += 0.5*self.a(pos[0,:],vel)*dt
        for i in range(1,int(time/dt)):

            pos_temp += vel*dt
            vel += self.a(pos_temp,vel)*dt
            #t[i] = t[i-1] + dt

            if i*dt >= timeBurn and not hohmann2Burned:
                print "Sirukerer"
                vel -= 0.5*self.a(pos[0,:],vel)*dt
                circBurn = self.circle(pos_temp,vel)
                vel += circBurn#self.circle(pos_temp,vel)#self.hohmann2(pos_temp,3000000,vel)
                print "Sirkulerings burn:",circBurn
                #vel = self.orbVel(norm(pos_temp))*vel/norm(vel)
                hohmann2Burned = True
                vel += 0.5*self.a(pos[0,:],vel)*dt


            if norm(pos_temp) < self.sys.radius[1]*1000:
                print "Crash"
                break


            if i%writingFreq == 0:
                pos[int(i/writingFreq),:] = pos_temp
                #print "At radius: ",norm(pos_temp) - self.sys.radius[1]*1000

                print (float(i)/int(time/(dt)))*100, "%            \r",
        print ""


        print "Final r[km]: ", (norm(pos[-1,:]))/1000. - self.sys.radius[1]
        # fig = plt.figure()
        # ax = fig.gca(projection = "3d")
        # # ax.plot(pos[:,0],pos[:,1],pos[:,2])
        # # plt.title("Orbit")
        # # plt.xlabel("x")
        # # plt.ylabel("y")
        # # plt.show()
        #
        # ax = fig.gca(projection='3d')
        # phi = np.linspace(0, 2 * np.pi, 100)
        # theta = np.linspace(0, np.pi, 100)
        # xm = self.sys.radius[1]*1000 * np.outer(np.cos(phi), np.sin(theta))
        # ym = self.sys.radius[1]*1000 * np.outer(np.sin(phi), np.sin(theta))
        # zm = self.sys.radius[1]*1000 * np.outer(np.ones(np.size(phi)), np.cos(theta))
        #
        #
        # ax.plot_surface(xm, ym, zm)


        #ax.plot(pos[:,0],pos[:,1],pos[:,2])
        plt.plot(pos[:,0],pos[:,1])
        #plt.gca().set_aspect('equal', adjustable='box')
        #plt.axis("equal")
        plt.show()

        print "Copy to instuctions: "
        print "----------"
        print "boost %g %g %g %g" %(startTime + 0.001, totalFirstBoost[0],totalFirstBoost[1],totalFirstBoost[2])
        print "boost %g %g %g %g" %(timeBurn, circBurn[0],circBurn[1],circBurn[2])
        print "----------"





seed = 75041

startPos = (np.array([ 50046912.943 ,  1922577.0044 ,  0]))
startV = np.array([-43.5115671208 ,  667.405339613 ,  0])


newPos = (np.array([ 1118641.33046991, -6317229.93852847,0.]))
newVel = (np.array([ 1994.27139912,  -422.71705031 ,0.]))
startTime = 90462


filename = "landing.txt"

orb = orbit(seed,filename,startPos,startV)
#orb = orbit(seed,filename,newPos,newVel)

hight = orb.findSafeHeight() + orb.sys.radius[1]*1000.
newHight = hight/8.0
print "Planet radius: ", orb.sys.radius[1]

#orb.findSafeHeight()
#h = np.linspace(0,210818,1000001)
# plt.plot(h,orb.adiabaticDensity(h))
# plt.show()
# plt.plot(h,orb.atmosDensity(h))
# plt.show()
#orb.simOrbit(160000,hight,0)
#orb.simOrbit(80000,newHight,startTime)
#print norm(np.array([2037990.98306102,  2053420.61273033]))
orb.lauch()
