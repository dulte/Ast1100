import numpy as np
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn
import sys


def calc_r(pos,time,posfunc):
    return norm(pos-posfunc(time)[:,1])

def calc_r_vec(pos,time,posfunc):
    return (pos-posfunc(time)[:,1])

def  calc_correction_burn(pos,tim,posfunc):
    r = (pos-posfunc(time+0.001)[:,1])
    r_norm = 1.0/(norm((pos-posfunc(time)[:,1])))*r
    return 100*r

def calc_planet_vel(time,func):
    h = 1.0/300000
    return (func(time+h)[:,1] - func(time)[:,1])/(h)

def calc_injection_burn(time,func,pos,planet_vel,sat_vel,p_mass):
    G = 4*np.pi**2

    r = pos - func(time)[:,1]

    #print norm(r)
    theta = np.arctan2(r[1],r[0])
    orb_vel = np.sqrt(G*p_mass/norm(r))
    #print "orb_vel: ",orb_vel
    #print "Vel diff",planet_vel-sat_vel

    return np.array([-orb_vel*np.sin(theta),orb_vel*np.cos(theta)]) - sat_vel + planet_vel

def print_theta(func,pos,time):
    theta_p = np.arctan2(func(time)[1,1],func(time)[0,1])
    theta_s = np.arctan2(pos[1],pos[0])

    print "Theta: ",theta_p - theta_s









seed = 75041
sys = AST1100SolarSystem(seed)

km_to_au = 6.685e-9

with open("planetPositions.npy", "rb") as npy:
    [planetPositions,t] = np.load(npy)

steps =  planetPositions[0,0,:].size
print steps
#t = np.linspace(0,20,steps)

planetPosFunction = inter.interp1d(t, planetPositions)

pos = np.array([0.80253929065 , 4.16457524661])
time = 8.81981666666

print calc_r(pos,time,planetPosFunction)
print calc_r_vec(pos,time,planetPosFunction)
print calc_correction_burn(pos,time,planetPosFunction)
print norm(calc_correction_burn(pos,time,planetPosFunction))



pos2 = np.array([0.792907737234 , 4.16748703398])
time2 = 8.82231666666

print "-----"
print calc_r_vec(pos2,time2,planetPosFunction)
print "When burning: ",calc_r(pos2,time2,planetPosFunction)
print "-----"

vel_inj = np.array([-3.77504600841 , 1.09819863799])
pos_inj = np.array([0.792907737234 , 4.16748703398])
time_inj = 8.82231666666

planet_vel = calc_planet_vel(time,planetPosFunction)
p_mass = sys.mass[1]


print "injection burn: " ,calc_injection_burn(time_inj,planetPosFunction,pos_inj,planet_vel,vel_inj,p_mass)
#print norm(calc_injection_burn(time_inj,planetPosFunction,pos_inj,planet_vel,vel_inj,p_mass))
print "-----"

time = 10.82231666666
pos = np.array([-4.1747684284 ,  -1.02819079768])

print "After 2 years: ",calc_r(pos,time,planetPosFunction)
print_theta(planetPosFunction,pos,time)
