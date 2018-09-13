import numpy as np
import random
from AST1100SolarSystem import AST1100SolarSystem
import matplotlib.pyplot as plt

seed = 4252
random.seed(seed)

#Constants
AU = 149597871 #[km]
G = 4 * np.pi**2
sek_yr = 365*24*60*60   #second in year

#AST1100SolarSystem:
sys = AST1100SolarSystem(seed)
massStar = sys.starMass     #[Solar masses]
massPlanets = sys.mass      #[Solar masses]
semimajor = sys.a

#Itterations
time = 300.      #yr (this long to confirm correct orbit over time)
n =  int(time*365) #number of time steps
dt = time/n #time step

#Find 4 most massive planets:
iBP = np.argsort(massPlanets)[-4:]  #Indices for the 4 biggest planets
massBigPlanets = [massPlanets[i] for i in iBP] #mass 4 biggest planets

#Arrays:
pos = np.zeros((2,5,n))    #position star + planets (from origo)
rvec_objekts = np.zeros((2,4)) #vector from star to each planet
vvec = np.zeros_like(pos)   #Velocity vector star + planets
a_planets = np.zeros((2,4)) #Acceleration planets
t = np.zeros(n)             #time

#Initial conditions planets:
pos[0,1::,0] = [sys.x0[i] for i in iBP]
pos[1,1::,0] = [sys.y0[i] for i in iBP]
vvec[0,1::,0] = [sys.vx0[i] for i in iBP]
vvec[1,1::,0] = [sys.vy0[i] for i in iBP]
#reserving pos[:,0,:] and vvec[:,0,:] for star

#Genrating initial velocitys for the star so total momentum is conserved
pTotP = 0
for j in range(1,5):
    pTotP += massBigPlanets[j-1]*vvec[:,j,0]
p_star = -pTotP
vvec[:,0,0] = p_star/massStar

#Initial-CM with star in middle:
M_tot = massStar + np.sum(massBigPlanets)
Rcm = 1/(M_tot) * (np.sum(massBigPlanets*pos[:,1::,0],axis=1))
#Moving pos of objects so CM in origo
for j in range(5):
    pos[:,j,0] += -Rcm #New initial pos when CM in origo
rvec_objekts[:,:] = pos[:,1::,0] #Initial vec between objects when cm in origo

#Integration loop:
for i in range(n-1):
    r_objekts = np.linalg.norm(rvec_objekts[:,:],axis=0) #Norm vec star-planet
    #Acceleration star:
    a_star = G*(np.sum(massBigPlanets*rvec_objekts[:,:]/r_objekts**3,axis=1))
    for k in range(2): #Acceleration for planets:
        a_planets[k,:] = -G*massStar*rvec_objekts[k,:]/r_objekts**3
    vvec[:,0,i+1] = vvec[:,0,i] + a_star * dt   #velocity star
    vvec[:,1::,i+1] = vvec[:,1::,i] + a_planets*dt  #velocity planets
    pos[:,0,i+1] = pos[:,0,i] + vvec[:,0,i+1] * dt  #pos star
    pos[:,1::,i+1] = pos[:,1::,i] + vvec[:,1::,i+1]*dt  #pos planets
    for j in range(4):  #updating vec from star to planets
        rvec_objekts[:,j] = pos[:,j+1,i] - pos[:,0,i]
    t[i+1] = t[i] + dt
print t[-1]
#Radial velocity:
def radialvStar(v,i,time):
    """
    Returns radial velocity with noise
    and timevector for star over time when line of sight along x-axis.
    v: velocity vector (including planets)
    i: inclination
    time: how long time to be calculated out of vector t, rounded up
    """
    indexStop = len(t[np.where(t<time)])-1 #index for t to get min required time
    if t[-1] < time:
        print ('Didnt reach required time, radial velocity shown over %g yrs'
        %t[-1])
    else:
        print ('Radial velocity shown over %g yrs' %t[indexStop])
    v_rad = v[0,0,0:indexStop] * np.sin(i)  #Radial velocity
    noise = np.random.normal(0,np.max(v_rad)/5.,len(v_rad))
    return v_rad + noise, t[0:indexStop]

i = np.pi/2 #inclination
v_rad,t_rad=(radialvStar(vvec,i,
            np.sqrt((np.max([semimajor[i] for i in iBP]))**3)+20))
#covers time for outer planet orbit
print np.sqrt(np.max([semimajor[i] for i in iBP])+20)
print t_rad[-1]

plt.figure(1)
plt.plot(t_rad,v_rad)
plt.xlabel('Time[yr]')
plt.ylabel('Radial velocity[AU/yr]')
plt.title('Radial velocity star')
plt.savefig('RadVelStar')
plt.show()
