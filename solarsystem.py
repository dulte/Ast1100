import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem

system = AST1100SolarSystem(1595,hasMoons=False)


#Variables for simulations
simulationTime = 1#In years
steps = 20000 #Steps per year
dt = 1./steps
time = 0

times_around = 0

G = 4*np.pi**2

use_many_body = False

#Vaaribles about the planets and star
N = system.numberOfPlanets

masses = []
postions = np.zeros([steps*simulationTime,N + 1,2]) #Number of datapoints, numb of planets and star, 2 coordinates
velocity = np.zeros([steps*simulationTime,N + 1,2])

temp_pos = np.zeros([N + 1,2])
temp_vel = np.zeros([N + 1,2])

masses = np.zeros(N+1)
periods = np.zeros(N+1)
acceleration = np.zeros([N+1,2])

sunMass = system.starMass
masses[N] = sunMass
print "mass: ", sunMass
print "Temp: ",system.temperature
print "Radius: ",system.starRadius
print "___"


#Initializing the arrays
for i in range(N):
    masses[i] = system.mass[i]
    print masses[i]
    print system.period[i]
    print "----"
    postions[0,i,0] = system.x0[i]
    postions[0,i,1] = system.y0[i]



    velocity[0,i,0] = system.vx0[i]
    velocity[0,i,1] = system.vy0[i]




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
    for step in range(steps*simulationTime-1):
        postions[step+1,:N,:] = postions[step,:N,:] + velocity[step,:N,:]*dt + 0.5*acceleration[:N,:]*dt*dt
        temp_vel[:N,:] = velocity[step,:N,:] + 0.5*acceleration[:N,:]*dt

        for p in range(N):
            acceleration[p,:] = -G*sunMass/(norm(postions[step+1,p,:])**3)*postions[step+1,p,:]
        #acceleration[:N,:] = -G*sunMass/(norm(postions[step,:N,:])**3)*postions[step,:N,:]


        #velocity[step+1,:N,:] = velocity[step,:N,:] + acceleration[:N,:]*dt
        #postions[step+1,:N,:] = postions[step,:N,:] + velocity[step+1,:N,:]*dt
        #print postions[step+1,:,:]

        velocity[step+1,:N,:] = temp_vel[:N,:] + 0.5*acceleration[:N,:]*dt

        if (np.arctan(postions[step,0,1]/postions[step,0,0]) < 0 and np.arctan(postions[step+1,0,1]/postions[step+1,0,0]) > 0 ):
            times_around += 1

        print (float(step)/(steps*simulationTime))*100, "%            \r",
    print ""
print times_around

for p in range(N-1):
    plt.plot(postions[-1,p,0],postions[-1,p,1],"o")
    plt.plot(postions[0,p,0],postions[0,p,1],"x")
    plt.plot(postions[:,p,0],postions[:,p,1])

plt.plot(postions[-1,N,0],postions[-1,N,1],"y*")
#plt.plot(postions[0,p,0],postions[0,p,1],"x")
plt.plot(postions[:,N,0],postions[:,N,1])

plt.show()
