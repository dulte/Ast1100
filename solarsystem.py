import numpy as np
import sys
from numpy.linalg import norm
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
seed = 69558 #daniehei5
system = AST1100SolarSystem(seed,hasMoons=False)


#Variables for simulations
simulationTime = 20#In years
steps = 20000 #Steps per year
dt = 1./steps

writeingFreq = 10

if (float(steps*simulationTime)/writeingFreq - steps*simulationTime/writeingFreq) != 0.0:
    print "change writeingFreq"
    sys.exit()


times_around = 0

G = 4*np.pi**2

use_many_body = True

#Vaaribles about the planets and star
N = system.numberOfPlanets

try:
    postions = np.zeros([2,N + 1,steps*simulationTime/writeingFreq]) #Number of datapoints, numb of planets and star, 2 coordinates
except:
    print "change writeingFreq"
temp_pos = np.zeros([2,N+1])
velocity = np.zeros([2,N + 1])#,steps*simulationTime])
sun_vel = np.zeros([2,steps*simulationTime/writeingFreq])
sun_pos = np.zeros([2,steps*simulationTime/writeingFreq])

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
    temp_pos[0,i] = system.x0[i]
    temp_pos[1,i] = system.y0[i]



    velocity[0,i] = system.vx0[i]
    velocity[1,i] = system.vy0[i]

    periods[i] = (2*np.pi/system.period[i])*360.
    print "Planet dtheta: ",periods[i]

    print "----"







if(use_many_body):

    mom = np.zeros(2)
    mom[0] = np.sum(velocity[0,:]*masses)
    mom[1] = np.sum(velocity[1,:]*masses)

    velocity[:,-1] = -mom/sunMass



    for i in xrange(N+1):
        acceleration[:,i] = 0
        rx = (temp_pos[0,:] - float(temp_pos[0,i]))**2
        ry = (temp_pos[1,:] - float(temp_pos[1,i]))**2
        other_p = np.power(rx+ry,3./2) > 1.e-5
        some = np.maximum(np.power(rx+ry,3./2),1.e-5)
        acceleration[0,i] = np.sum(other_p*(G*masses[:]/some*(temp_pos[0,:] - float(temp_pos[0,i]))))
        acceleration[1,i] = np.sum(other_p*(G*masses[:]/some*(temp_pos[1,:] - float(temp_pos[1,i]))))


    #
    # for i in range(N+1):
    #     acceleration[:,i] = 0
    #     for j in range(N+1):
    #         r = norm(temp_pos[:,i]-temp_pos[:,j])
    #         if (r>1.e-4):
    #             acceleration[:,i] += -G*masses[j]/(r**3)*(temp_pos[:,i]-temp_pos[:,j])










    velocity += 0.5*acceleration*dt

    for step in xrange(steps*simulationTime-1):

        temp_pos += velocity*dt

        for i in xrange(N+1):
            acceleration[:,i] = 0
            rx = (temp_pos[0,:] - float(temp_pos[0,i]))**2
            ry = (temp_pos[1,:] - float(temp_pos[1,i]))**2
            other_p = np.power(rx+ry,3./2) > 1.e-5
            some = np.maximum(np.power(rx+ry,3./2),1.e-5)
            acceleration[0,i] = np.sum(other_p*(G*masses[:]/some*(temp_pos[0,:] - float(temp_pos[0,i]))))
            acceleration[1,i] = np.sum(other_p*(G*masses[:]/some*(temp_pos[1,:] - float(temp_pos[1,i]))))

        # for i in range(N+1):
        #     acceleration[:,i] = 0
        #     for j in range(N+1):
        #         r = norm(temp_pos[:,i]-temp_pos[:,j])
        #         if (r>1.e-4):
        #             acceleration[:,i] = -G*masses[j]/(r**3)*(temp_pos[:,i]-temp_pos[:,j])

        velocity += acceleration*dt

        if ((step+1)%writeingFreq == 0):
            postions[:,:,step/writeingFreq + 1] = temp_pos[:,:]
            sun_pos[:,step/writeingFreq +1] = temp_pos[:,-1]
            sun_vel[:,step/writeingFreq +1] = velocity[:,-1]

        print (float(step)/(steps*simulationTime))*100, "%            \r",
    print ""









else:
    acceleration = -G*sunMass/(norm(temp_pos[:,:N],axis = 0)**3)*temp_pos[:,:N]
    velocity[:,:N] += 0.5*acceleration[:,:N]*dt


    for step in xrange(steps*simulationTime-1):


        temp_pos[:,:N] += velocity[:,:N]*dt

        acceleration = -G*sunMass/(norm(temp_pos[:,:N],axis = 0)**3)*temp_pos[:,:N]


        velocity[:,:N] += acceleration[:,:N]*dt


        if ((step+1)%writeingFreq == 0):
            postions[:,:,step/writeingFreq + 1] = temp_pos[:,:]

        print (float(step)/(steps*simulationTime))*100, "%            \r",
    print ""

#postions[:,:,-1] = temp_pos

#system.checkPlanetPositions(postions,simulationTime,steps/writeingFreq)
#system.orbitXml(postions[:,:N,:],time)
for p in range(N):
    plt.plot(postions[0,p,-1],postions[1,p,-1],"o")
    plt.plot(postions[0,p,0],postions[1,p,0],"x")
    plt.plot(postions[0,p,:],postions[1,p,:])

plt.plot(postions[0,N,-1],postions[1,N,-1],"y*")
#plt.plot(postions[0,p,0],postions[0,p,1],"x")
plt.plot(postions[0,N,:],postions[1,N,:])

plt.show()
plt.plot(sun_vel[0])
plt.show()
