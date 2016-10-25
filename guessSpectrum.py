import numpy as np
#from PIL import Image
import matplotlib.pyplot as plt
#from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn


class spectral:

    def __init__(self, specFile, sigmaFile):

        try:
            self.spectral_lines = np.load(specFile + ".npy")
            self.sigmas = np.load(sigmaFile + ".npy")
            print "Npy files loaded"
        except:
            self.spectral_lines = np.loadtxt(specFile + ".txt")
            self.sigmas = np.loadtxt(sigmaFile + ".txt")
            print "Text files loaded"
            np.save(specFile + ".npy",self.spectral_lines)
            np.save(sigmaFile + ".npy",self.sigmas)
            print "Files saved as npy"



        self.cO2 = np.array([630,690,760])
        self.cH2O = np.array([720,820,940])
        self.cCO2 = np.array([1400,1600])
        self.cCH4 = np.array([1660,2200])
        self.cCO = np.array([2340])
        self.cN2O = np.array([2870])

        self.mO2 =  31.99880
        self.mCO2 = 44.0095
        self.mH2O = 18.01528
        self.mCH4 = 16.0425
        self.mCO = 28.0101
        self.mN2O = 44.01280


    def plot_exact(self):
        plt.plot(self.spectral_lines[:,0],self.spectral_lines[:,1])
        plt.show()


    def exp_expression(self,Lambda,lambdaC,f_min,sigma):
        return np.exp(-(Lambda - lambdaC)**2/(2*sigma**2))

    def model(self,Lambda,lambdaC,f_min,sigma,f_max):
        return (f_min-f_max)*np.exp(-(Lambda - lambdaC)**2/(2*sigma**2))+ f_max

    def estimate(self,lambdaC,f_min,f_max,m):


        Tmin = 150.
        Tmax = 450.


        sigma = np.linspace((1.0/20)*self.sigma(Tmin,lambdaC,m),2*self.sigma(Tmax,lambdaC,m),30)
        Fmin = np.linspace(f_min,1,30)
        lambda_0 = np.linspace(lambdaC-self.find_shift(10000,lambdaC),lambdaC+self.find_shift(10000,lambdaC),300)
        print self.find_shift(10000,lambdaC)

        interesting_Indexes = np.logical_and(self.spectral_lines[:,0] < (lambdaC+0.1),self.spectral_lines[:,0] > (lambdaC-0.1))

        interesting_spec = self.spectral_lines[interesting_Indexes,0]
        interesting_flux = self.spectral_lines[interesting_Indexes,1]
        interesting_sigma = self.sigmas[interesting_Indexes,1]
        print self.calc_temp(sigma[-1],lambdaC,m)

        removedZeros = np.where( self.spectral_lines[:,0]*interesting_Indexes == 0,1e6,self.spectral_lines[:,0]*interesting_Indexes )

        index_first = np.argmin(removedZeros)
        index_last = (np.argmax(self.spectral_lines[:,0]*interesting_Indexes))
        print index_last,index_first
            #self.f_max = f_max
        N_inv = 1.0/self.spectral_lines.shape[0]

        best_parameters = np.zeros(3)
        best_chi = 1e16

        print "Beginning estimation of spectral line at wave lenght %g nm" %lambdaC

        counter = 0
        for f in Fmin:
            for sig in sigma:
                for lam in lambda_0:
                    model = self.model(interesting_spec,lam,f,sig,f_max)
                    chi_squared = np.sum((interesting_flux-model)**2/(interesting_sigma**2))
                    #print chi_squared
                    if (chi_squared < best_chi):
                        best_chi = chi_squared
                        best_parameters = np.array([lam,f,sig])
            print counter
            counter += 1




        return np.append(np.array([lambdaC,index_first,index_last]),best_parameters)

    def find_shift(self,v,lamb):
        c = 299792458
        return lamb*(float(v)/c)

    def find_paramters(self):




        f_min = 0.7
        f_max = 1


        T = 300.

        totalMass = 0
        nMol = 0

        # self.paraCH4 = self.parameter_dump(self.cCH4,self.mCH4,f_min,f_max,T)
        # self.find_gass(self.cCH4,self.mCH4,"noe")
        #self.paraH2O = self.parameter_dump(self.cH2O,self.mH2O,f_min,f_max,T)
        #self.paraN2O = self.parameter_dump(self.cN2O,self.mN2O,f_min,f_max,T)
        #self.paraCO = self.parameter_dump(self.cCO,self.mCO,f_min,f_max,T)


        newMass,newnMol = self.find_gass(self.cH2O,self.mH2O,"noe")
        totalMass += newMass
        nMol += newnMol

        newMass,newnMol = self.find_gass(self.cO2,self.mO2,"noe")
        totalMass += newMass
        nMol += newnMol

        newMass,newnMol = self.find_gass(self.cN2O,self.mN2O,"noe")
        totalMass += newMass
        nMol += newnMol

        newMass,newnMol = self.find_gass(self.cCO,self.mCO,"noe")
        totalMass += newMass
        nMol += newnMol

        newMass,newnMol = self.find_gass(self.cCO2,self.mCO2,"noe")
        totalMass += newMass
        nMol += newnMol

        newMass,newnMol = self.find_gass(self.cCH4,self.mCH4,"noe")
        totalMass += newMass
        nMol += newnMol

        print "Mean molar mass is: ",totalMass/nMol

        #print self.paraO2


    def parameter_calc(self,clamb,m,f_min,f_max,T):
        para = np.zeros((6,clamb.size))
        for i in range(clamb.size):
            para[:,i] = self.estimate(clamb[i],f_min,f_max,m)
            print "Temp: ", self.calc_temp(para[2,i],para[0,i],m)
        return para

    def parameter_dump(self,clamb,m,f_min,f_max,T):
        para = np.zeros((6,clamb.size))
        for i in range(clamb.size):
            para[:,i] = self.estimate(clamb[i],f_min,f_max,m)
            np.savetxt("parameters_%g.txt" % (clamb[i]),para[:,i])



    def sigma(self,T,lamb,m):
        u = 1.66053904e-27
        c = 299792458
        k = 1.38064852e-23
        return (1.0/(np.sqrt(8*np.log(2))))*(2.0*lamb/c)*np.sqrt(2*k*T*np.log(2)/(m*u))




    def plotEstimate(self):

        plt.plot(self.spectral_lines[:,0],self.model(self.spectral_lines[:,0],self.lambdaC,self.f_min,self.sigma,self.f_max))
        plt.show()

    def calc_temp(self,sigma,lamb,m):
        u = 1.66053904e-27
        c = 299792458
        k = 1.38064852e-23
        return (sigma/((1.0/(np.sqrt(8*np.log(2)))))*(1.0/(2.0*lamb/c)))**2*(u*m)*(1.0/(2*k*np.log(2)))


    def find_gass(self,clamb,m,gasname):
        mh = 2.01588
        f_max = 1
        massFraction = []
        for i in range(len(clamb)):
            try:
                parameters = np.loadtxt("parameters_%g.txt" %clamb[i])
                lamb = self.spectral_lines[int(parameters[1]):int(parameters[2])+1,0]
                flux = self.spectral_lines[int(parameters[1]):int(parameters[2])+1,1]
                plt.plot(lamb,flux)
                plt.plot(lamb, self.model(lamb,parameters[3],parameters[4],parameters[5],f_max))
                plt.show()
                ans = raw_input("Er dette ekte (y/n)")
                if (ans == "y"):
                    massFraction.append(m)
                    print "Found something at wavelenght %g nm" %clamb[i]
                    print "Temp: ", self.calc_temp(parameters[5],clamb[i],m)
            except:
                print "Noe skjedde"

        return sum(massFraction),len(massFraction)







sp = spectral("spectrum_seed41_600nm_3000nm","sigma_noise")
sp.find_paramters()
#sp.estimate(200,.001)
#sp.plotEstimate()
