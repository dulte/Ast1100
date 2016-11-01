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


seed = 75041
sys = AST1100SolarSystem(seed)

pos = -np.array([-3051899.65683044 , -175360.35571522  ,      0.            ])

phi = np.arctan2(pos[1],pos[0])
theta = np.arccos(pos[2]/float(norm(pos)))

print "r [km]: ", (norm(-pos)/1000. - sys.radius[1])
print "Theta: ",theta
print "Phi: ",phi
