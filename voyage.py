import numpy as np
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn
import sys





def calc_num_boxes(F,v):
    return (v*1100)/(10*60*F)

def calc_escape_m_per_sec(system):
    G = 4*np.pi**2
    km_to_au = 6.685e-9
    return np.sqrt(2*G*system.mass[0]/(system.radius[0]*km_to_au))*4743.7173611


seed = 75041
sys = AST1100SolarSystem(seed)
force_per_box = 1.72650965887e-09
numb_escaping_per_sec = 6.4252e+13
mass_lost_per_sec = 2.120316e-13

km_to_au = 6.685e-9

with open("planetPositions.npy", "rb") as npy:
    [planetPositions,t] = np.load(npy)

steps =  planetPositions[0,0,:].size
print steps
#t = np.linspace(0,20,steps)

planetPosFunction = inter.interp1d(t, planetPositions)
numberOfPlanets = 9
# for p in range(numberOfPlanets):
#     plt.plot(planetPositions[0,p,-1],planetPositions[1,p,-1],"o")
#     plt.plot(planetPositions[0,p,0],planetPositions[1,p,0],"x")
#     plt.plot(planetPositions[0,p,:],planetPositions[1,p,:])
#
# for ti in t[:int(steps/10)]:
#     plt.plot(planetPosFunction(ti)[0,0],planetPosFunction(ti)[1,0])
#
#
#
#
# plt.axis("equal")
# plt.show()




time = 6.71585
startPos = np.array([0.57172522, -2.99619899])#np.array([0.5717085,  -2.99620209])


angle = -1.3822742546
print int(steps*time/20)
print steps
print planetPositions[:,0,int(steps*time/20)]
print planetPosFunction(t[int(steps*time/20)])[:,0]
print planetPosFunction(time)[:,0]
print "Func start pos: ", planetPosFunction(time)[:,0] + np.array([-np.sin(angle),np.cos(angle)])*sys.radius[0]*km_to_au
print t[int(steps*time/20)]
print t[int(steps/time)+1]
print "-----"
print "Array start pos: ",planetPositions[:,0,int(steps*time/20)] + np.array([-np.sin(angle),np.cos(angle)])*sys.radius[0]*km_to_au
#print np.array([0.57165763, -2.99621188])+ np.array([-np.sin(angle),np.cos(angle)])*sys.radius[0]*km_to_au
r = norm(planetPositions[:,0,int(steps*time/20)] - startPos)#norm(planetPosFunction(time)[:,0] - startPos)

startR = r/km_to_au

print sys.radius[0]
print 0.02*sys.radius[0]*km_to_au
print abs(r - sys.radius[0]*km_to_au)
if abs(r - sys.radius[0]*km_to_au) < 0.02*sys.radius[0]*km_to_au:
    print "Check passed!"
escape_velocity = calc_escape_m_per_sec(sys)
print calc_num_boxes(force_per_box,escape_velocity)
sys.sendSatellite("instructions2.txt")
