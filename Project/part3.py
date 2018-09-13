import numpy as np
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
import matplotlib.pyplot as plt

seed = 4252
np.random.seed(seed)

inFile = open('positionsHomePlanet.npy','rb')

try:
    #Itterations from part 2
    time2 = 40.   #[yr]
    n2 = int(time2*365*24)       #time steps from part 2
    dt2 = time2/n2
    planetPos = np.load(inFile)
    times = np.linspace(0,time2,n2)
finally:
    inFile.close()
posFunc = inter.interp1d(times,planetPos)

#Constants
AU = 149597871. #[km]
G = 4 * np.pi**2

#Itteration constants
stepsperyear= 10**4
time_launch = 45./(365*24*60)   #time launch. 45 min. [yr]
n = int((time2)*stepsperyear) #time steps long travel
n_ = 200000 #time steps launch
dt = time2/n #time step stellar
dtt = 10**(-6) #time step launch

#AST1100SolarSystem info
sys = AST1100SolarSystem(seed,hasMoons = False)
massStar = sys.starMass     #[Solar masses]
radiusStar = sys.starRadius #[km]
tempStar = sys.temperature  #[K]

N = sys.numberOfPlanets     #Number of planets
massPlanets = sys.mass      #[Solar masses]
radiusPlanets = sys.radius  #[km]
majoraxis = sys.a           #major axis [AU]

print radiusPlanets[6]/AU
#Arrays:
r_planets = []
timeplanets = np.linspace(0,time2,500)
for i in range(9):
    r_i = posFunc(timeplanets)[:,i] #position planets 0 through 8.
    r_planets.append(r_i)

#Functions
def velPlanet(N,time):
    """
    Calculates velocity planet by calculating difference in
    position over 3 different time steps periods in x and y direction.
    Then use the mean of the three velocitys to set planet velocity.
    N: planet number
    time: actuall time of interest
    Returns velocity in x and y direction.
    """
    t = time
    v1 = (posFunc(t+times[1])[:,N]-posFunc(t)[:,N])/times[1]
    v2 = (posFunc(t+2*times[1])[:,N]-posFunc(t)[:,N])/(2*times[1])
    v3 = (posFunc(t+3*times[1])[:,N]-posFunc(t)[:,N])/(3*times[1])
    return (v1+v2+v3)/3.

def launch(time,phi,gamma,v_diff):
    """
    Launches satellite from home planet. Only considers gravity from home planet
    and the star. Calculates as long as the satellite is close to home.
    time: time of launch
    phi: Angel from radial direction star-home to launch position
    v_diff: scalar to multiply escape velocity
    gamma: angel velocity from directly upwards direction
    Returns: position satellite at time, velocity satellite at time and time
    when launch complete.
    """
    rvecSatlaunch = np.zeros((2,n_))
    vSat = np.zeros((2,1))
    #Initial values:
    poshome = np.zeros((2,n_))
    poshome[:,0] = posFunc(time)[:,0]    #position home planet at time [AU]
    alpha = np.arctan2(poshome[1,0],poshome[0,0]) #Angle home from x-axis
    radhome = float(radiusPlanets[0]/AU)           #Radius home planet [AU]
    #Angles:
    theta = alpha + phi
    beta = theta + gamma
    #Initial pos sat: Max distance from star when theta = 0,
    #Evt turned along surface with angle phi from radial direction.
    rvecSatlaunch[:,0]=(poshome[:,0]+
        radhome*np.array([np.cos(theta),np.sin(theta)]))
    rvec_sat_home = - poshome[:,0] + rvecSatlaunch[:,0] #vec from home to sat
    r_sat_home = np.linalg.norm(rvec_sat_home) #Initial distance home - sat
    vhome = velPlanet(0,time)   #Velocity home planet at time
    v_esc = np.sqrt(2*G*massPlanets[0]/radhome) #Escape velocity from home planet
    #Initial velocity satellite
    vSat = vhome + (v_esc*v_diff)*np.array([np.cos(beta),np.sin(beta)])
    i = 0
    t = time
    #Itterations:
    while r_sat_home <= 0.1:
        poshome[:,i] = posFunc(t)[:,0]    #position home planet at time [AU]
        rvec_sat_home = -poshome[:,i]+ rvecSatlaunch[:,i] #vec from home to sat
        r_sat_home = np.linalg.norm(rvec_sat_home) #distance home - sat
        a = (-G*massPlanets[0]*(rvec_sat_home)/r_sat_home**3 -
            G*massStar*rvecSatlaunch[:,i]/np.linalg.norm(rvecSatlaunch[:,i])**3)
        vSat = vSat + a*dtt
        rvecSatlaunch[:,i+1] = rvecSatlaunch[:,i] + vSat*dtt
        t += dtt
        poshome[:,i+1] = posFunc(t)[:,0]
        i += 1

    return rvecSatlaunch[:,0:i],t, vSat

def stellar(resultLaunch,boost1,boost2):
    """
    Takes results from launch and continues the satellites journey until it
    has pasted the orbit of an outer planet, or too much time has passed.
    Calculates the effect of gravity from each planet and the star.
    Calculates the position of the planets at the time when it is closest to
    satellite at any time during the journey.
    boost: array with time of boost, and new x and y velocitys
    Returns satellite position vector and points when the planets are closest.
    """
    #Initial setup:
    dt = time2/n                #[yr]
    rvecSat = np.zeros((2,n))
    rvecSat[:,0] = resultLaunch[0][:,-1]
    t = resultLaunch[1]
    vSat = resultLaunch[2]
    i = 0
    posPlanets = posFunc(t)  #inital pos planets
    plpos_closest = posPlanets  #Positions planets when closest to satellite
    satpos_closest = np.zeros((2,9))   #Sat pos when closest to planets
    #Smallest distance from planets to satellite
    r_is_smallest = (np.linalg.norm(
    np.array([rvecSat[:,0]-posPlanets[:,k] for k in range(9)]),axis=1))
    t_close = np.zeros(9)   #Times when closest planets-satellite
    while( i < n-1 and np.linalg.norm(rvecSat[:,i]) < majoraxis[4]
    and t<resultLaunch[1]+5.4):
        if i == boost1[0]:  #Boost when close to target
            print 't boost 1',t
            deltavSat1 = np.array([boost1[1]-vSat[0],boost1[2]-vSat[1]])
            print deltavSat1,np.linalg.norm(deltavSat1)
            print 'vSat before boost1',vSat
            vSat = np.array([boost1[1],boost1[2]])
            print 'vSat after boost1',vSat

        if i == boost2[0]:
            print 't boost 2',t
            v_so = stabilvel(r_is[6])
            posrel = rvecSat[:,i] - posPlanets[:,6]
            vrel = vSat - velPlanet(6,t)
            print 'satpos boost2',rvecSat[:,i]
            print 'r_is',r_is[6]
            print 'v stabil',np.linalg.norm(v_so+velPlanet(6,t))
            print 'v stabil rel planet', np.linalg.norm(v_so)
            print 'vSat before boost2', np.linalg.norm(vSat)
            phi_ = np.arctan2(posrel[1],posrel[0])
            chi_ = np.arctan2(vrel[1],vrel[0])
            vSat1 = np.array([-v_so*np.sin(phi_),v_so*np.cos(phi_)])+velPlanet(6,t)
            print 'Delta v boost 2',vSat1-vSat,np.linalg.norm(vSat1-vSat)
            vSat = vSat1
            print 'vSat after boost2',vSat,np.linalg.norm(vSat1)
            dt = dtt
        if i > boost2[0]: #Only consider gravity form target and star
            posPlanets = posFunc(t)  #pos planets at time
            #Vectors from planets to satellite:
            rvec_is = np.array([rvecSat[:,i]-posPlanets[:,k] for k in range(9)])
            #Length from planets to satellite
            r_is = np.linalg.norm(rvec_is,axis=1)
            a = -G*massStar*rvecSat[:,i]/(np.linalg.norm(rvecSat[:,i])**3)
            a += -G*massPlanets[6]*rvec_is[6,:]/(r_is[6]**3)
            vSat = vSat + a*dt

        else:
            posPlanets = posFunc(t)  #pos planets at time
            #Vectors from planets to satellite:
            rvec_is = np.array([rvecSat[:,i]-posPlanets[:,k] for k in range(9)])
            #Length from planets to satellite
            r_is = np.linalg.norm(rvec_is,axis=1)
            #Acceleration from star:
            a = -G*massStar*rvecSat[:,i]/(np.linalg.norm(rvecSat[:,i])**3)
            for k in range(9):  #Update a to include all planets
                a += -G*massPlanets[k]*rvec_is[k,:]/(r_is[k]**3)
            vSat = vSat + a*dt
        for k in range(9):
            if r_is[k] < r_is_smallest[k]:  #Update closest positions
                plpos_closest[:,k] = posPlanets[:,k]
                satpos_closest[:,k] = rvecSat[:,i]
                r_is_smallest[k] = r_is[k]
                t_close[k] = t
                if k == 6:
                    vSat_at_t = vSat
                    vTarget_at_t = velPlanet(6,t)
                    it = i
        rvecSat[:,i+1] = rvecSat[:,i] + vSat*dt
        i += 1
        t += dt
    print 'Iteration number closest',it
    print 'vSat closest',vSat_at_t, np.linalg.norm(vSat_at_t)
    print 'vtarget closest',vTarget_at_t
    print 'Distance to target min',r_is_smallest[6]
    return rvecSat[:,0:i],plpos_closest,satpos_closest,t_close,r_is_smallest[6]

def stabilorbit(R,m,M,k):
    """
    R: rvecSat
    m: massPlanet
    M: mass star
    k: force scale factor
    """
    R = np.linalg.norm(R)
    r = R*np.sqrt(m/(M*k))
    return r

def stabilvel(r):
    v_stabil = np.sqrt(massPlanets[6]*G/r)
    return v_stabil

def plot():
    """
    Plots the Satellite journey
    """
    fig = plt.figure('Satellite journey')
    ax = plt.subplot(111)
    ax.plot(posSat[0,:],posSat[1,:],label='pos sat',linewidth = 2.0)
    ax.hold('on')
    planets_toplot = [0,1,3,4,6,7,8]
    for k in planets_toplot:
        ax.plot(r_planets[k][0,:],r_planets[k][1,:],
        label='planet %i'%k)
    for k in planets_toplot:
        ax.plot(plclose[0,k],plclose[1,k],'o')
    for k in planets_toplot:
        if all(satclose[:,k])==0:
            ax.plot(satclose[0,k],satclose[1,k])
        else:
            ax.plot(satclose[0,k],satclose[1,k],'o')
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.set_title('Satellite journey')
    fig.savefig('Satellite_journey%.3e.png'%disttarget)

if __name__ == '__main__':
    test = launch(4,0,0,1.25)
    rvecSatlaunch = test[0]
    boost1 = [39850,-2.4,2.4]
    boost2 = [40322]
    posSat,plclose,satclose,t_close,disttarget = stellar(test,boost1,boost2)
    print 'satclose',satclose[:,6],np.linalg.norm(satclose[:,6])
    r_orbit= stabilorbit(satclose[:,6],massPlanets[6],massStar,10)
    v_orbit = stabilvel(r_orbit)
    print 'time closest',t_close[6]
    #Stabil orbit
    print 'distance from target for stable orbit',r_orbit
    print 'pos sat when closest',satclose[:,6]
    print 'pos target when closest',plclose[:,6]
    print 'velocity stable orbit',v_orbit
    plot()
    print 'Done'
    plt.show()
