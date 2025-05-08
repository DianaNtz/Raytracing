"""
The code below was written by @author: https://github.com/DianaNtz and is an
implementation of the 2nd order Runge Kutta time integration method for the
geodesic equation of a Schwarzschild black hole.
For more details on the code requirements and licences see
https://github.com/DianaNtz/Raytracing.
"""
#importing libraries
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
#setting up initial values
ntau=10000

M=1

r1n=3.0
r2n=0.0

t1n=0
t2n=0.1

phi1n=0
phi2n=np.sqrt((1-2*M/r1n)/r1n**2*t2n**2)

dtau=(2.0*np.pi/phi2n)/(ntau-1)

taun=0

#setting up functions for time integration
def ft1(t1,t2,r1,r2,phi1,phi2):
    return t2
def fr1(t1,t2,r1,r2,phi1,phi2):
    return r2
def fphi1(t1,t2,r1,r2,phi1,phi2):
    return phi2

def ft2(t1,t2,r1,r2,phi1,phi2):
    return -2*M/(r1*(r1-2*M))*t2*r2

def fr2(t1,t2,r1,r2,phi1,phi2):
    return M/(r1*(r1-2*M))*r2**2-M*(r1-2*M)/r1**3*t2**2+(r1-2*M)*phi2**2

def fphi2(t1,t2,r1,r2,phi1,phi2):

    return -2/r1*phi2*r2

tau=np.zeros(ntau)

r1=np.zeros(ntau)
r2=np.zeros(ntau)

t1=np.zeros(ntau)
t2=np.zeros(ntau)

phi1=np.zeros(ntau)
phi2=np.zeros(ntau)
