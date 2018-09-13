import numpy as np
from AST1100SolarSystem import AST1100SolarSystem
import random
import matplotlib.pyplot as plt
seed = 4252
random.seed(seed)

#Constants
AU = 149597871 #[km]
G = 4 * np.pi**2
sek_yr = 365*24*60*60   #second in year

#Itterations
time = 40.   #[yr]
n = int(time*365*24)       #time steps
dt = time/n


#AST1100SolarSystem info
sys = AST1100SolarSystem(seed,hasMoons = False)

massStar = sys.starMass     #[Solar masses]
radiusStar = sys.starRadius #[km]
tempStar = sys.temperature  #[K]

N = sys.numberOfPlanets     #Number of planets
massPlanets = sys.mass      #[Solar masses]
radiusPlanets = sys.radius  #[km]
majoraxis = sys.a           #major axis [AU]
period = sys.period         #rotational period [earth days]

#Arrays
rvec = np.zeros((2,N,n))
vvec = np.zeros((2,N,n))
t = np.zeros(n)
#test_r = np.zeros((N,n))

#Initial conditions
rvec[:,:,0] = [sys.x0,sys.y0]     #position vector
vvec[:,:,0] = [sys.vx0,sys.vy0]   #Velocity vector

tol = 10**(-5)
time_itterated = 0
for i in range(n-1):  #Itterating time steps
    r = np.linalg.norm(rvec[:,:,i],axis=0)
    a = - G * massStar * rvec[:,:,i] / r**3
    vvec[:,:,i+1] = vvec[:,:,i] + a * dt
    rvec[:,:,i+1] = rvec[:,:,i] + vvec[:,:,i+1] * dt
    t[i+1] = t[i] + dt
    time_itterated += dt

#Prints:
print 'Star: mass, radius, temp', massStar, radiusStar, tempStar
print 'Home: mass, radius', massPlanets[0],radiusPlanets[0]
'''
fig = plt.figure('Planetbaner')
ax = plt.subplot(111)
ax.plot(0,0)
planets_toplot = [0,1,3,4,6,7,8]
for i in planets_toplot:
    ax.plot(rvec[1,i,::100],rvec[0,i,::100], label = 'Planet %i'%i)
ax.plot(rvec[1,2,::100],rvec[0,2,::100], label = 'Planet 2')
ax.plot(rvec[1,5,::100],rvec[0,5,::100], label = 'Planet 5')
ax.set_xlabel('x-position [AU]')
ax.set_ylabel('y-position [AU]')
ax.legend(loc = 'upper center', bbox_to_anchor=(0.5, 1.05),
        ncol=3,fancybox=True,shadow=True)
fig.savefig('Planetbaner_part2.png')
plt.show()
'''
#sys.orbitXml(rvec[:,:,::100],t[::100])
sys.checkPlanetPositions(rvec,int(time),int(n/time))
