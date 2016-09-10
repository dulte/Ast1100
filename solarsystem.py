import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
seed = 69558 #daniehei5
system = AST1100SolarSystem(seed,hasMoons=False)


#Variables for simulations
simulationTime = 20#In years
steps = 20000 #Steps per year
dt = 1./steps


times_around = 0

G = 4*np.pi**2

use_many_body = False

#Vaaribles about the planets and star
N = system.numberOfPlanets
postions = np.zeros([2,N + 1,steps*simulationTime]) #Number of datapoints, numb of planets and star, 2 coordinates
velocity = np.zeros([2,N + 1])#,steps*simulationTime])

masses = np.zeros(N+1)
periods = np.zeros(N+1)
acceleration = np.zeros([2,N+1])

time = np.linspace(0,simulationTime,int(steps*simulationTime))
print time

sunMass = system.starMass
masses[N] = sunMass
print "mass: ", sunMass
print "Temp: ",system.temperature
print "Radius: ",system.starRadius
print "___"
print len(system.mass)

#Initializing the arrays
for i in range(N):
    masses[i] = system.mass[i]
    print "Planet %s " %str(i+1)
    print "Planet mass: ",masses[i]
    print "Planet radius: ",system.radius[i]

    postions[0,i,0] = system.x0[i]
    postions[1,i,0] = system.y0[i]



    velocity[0,i] = system.vx0[i]
    velocity[1,i] = system.vy0[i]

    periods[i] = (2*np.pi/system.period[i])*360.
    print "Planet dtheta: ",periods[i]

    print "----"
if(use_many_body):

    for step in range(steps*simulationTime-1):

        for i in range(N+1):
            acceleration[i,:] = 0
            for j in range(N+1):
                r = np.sqrt((postions[step,i,0]-postions[step,j,0])**2 + (postions[step,i,1]-postions[step,j,1])**2)
                if (r>1.e-4):
                    acceleration[i,0] = -G*masses[j]/(r**3)*(postions[step,i,0]-postions[step,j,0])

                    acceleration[i,1] = -G*masses[j]/(r**3)*(postions[step,i,1]-postions[step,j,1])

        velocity[step+1,:,0] = velocity[step,:,0] + acceleration[:,0]*dt
        postions[step+1,:,0] = postions[step,:,0] + velocity[step+1,:,0]*dt

        velocity[step+1,:,1] = velocity[step,:,1] + acceleration[:,1]*dt
        postions[step+1,:,1] = postions[step,:,1] + velocity[step+1,:,1]*dt

        if (np.arctan(postions[step,0,1]/postions[step,0,0]) < 0 and np.arctan(postions[step+1,0,1]/postions[step+1,0,0]) > 0 ):
            times_around += 1

        print (float(step)/(steps*simulationTime))*100, "%            \r",
    print ""
else:
    print velocity[:,:N]
    for p in xrange(N):
        acceleration[:,p] = -G*sunMass/(norm(postions[:,p,0])**3)*postions[:,p,0]
    velocity[:,:N] += 0.5*acceleration[:,:N]*dt


    for step in xrange(steps*simulationTime-1):


        postions[:,:N,step+1] = postions[:,:N,step] + velocity[:,:N]*dt
        for p in xrange(N):
            acceleration[:,p] = -G*sunMass/(norm(postions[:,p,step+1])**3)*postions[:,p,step+1]


        velocity[:,:N] += acceleration[:,:N]*dt




        #velocity[:,:N,step+1] = temp_vel[:,:N] + 0.5*acceleration[:,:N]*dt



        if (np.arctan(postions[1,0,step]/postions[0,0,step]) < 0 and np.arctan(postions[1,0,step+1]/postions[0,0,step+1]) > 0 ):
            times_around += 1

        print (float(step)/(steps*simulationTime))*100, "%            \r",
    print ""
print times_around

system.checkPlanetPositions(postions,simulationTime,20000)
#system.orbitXml(postions[:,:N,:],time)
for p in range(N):
    plt.plot(postions[0,p,-1],postions[1,p,-1],"o")
    plt.plot(postions[0,p,0],postions[1,p,0],"x")
    plt.plot(postions[0,p,:],postions[1,p,:])

plt.plot(postions[0,N,-1],postions[1,N,-1],"y*")
#plt.plot(postions[0,p,0],postions[0,p,1],"x")
plt.plot(postions[0,N,:],postions[1,N,:])

plt.show()
