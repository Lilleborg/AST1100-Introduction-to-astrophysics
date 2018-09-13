from AST1100ShortcutSystem import AST1100SolarSystem
import numpy as np
from funcpart3 import PlanetTemp, FluxSB, Luminosity
import matplotlib.pyplot as plt

seed = 4252

#Constants
AU = 149597871.*10**3 #[m]
pi = np.pi
G = 6.67428e-11         #[m**3/(kg*s)]
solar_mass = 1.98892e30 #[kg]
h = 6.6260693*10**(-34)         #Plancks constant [Js]
c = 299792458.  #Speed of light [m/s]
mu = 30         #mean molecular weight
k = 1.38 * 10**(-23)            #Boltzmanns k
m_H = 1.6817e-27                #Hydrogen mass [kg]
m_sat = 1100                    #Mass of satellite [kg]

#AST1100SolarSystem info
sys = AST1100SolarSystem(seed,hasMoons = True)
massStar = sys.starMass     #[Solar masses]
radiusStar = sys.starRadius #[km]
tempStar = sys.temperature  #[K]
massPlanets = sys.mass      #[Solar masses]
radiusPlanets = sys.radius  #[km]
majoraxis = sys.a           #major axis [AU]
gamma = 1.4                 #Adiabatic index

#In SI-Units target planet:
masstarget = massPlanets[6]*solar_mass  #Mass target planet[kg]
radtarget = radiusPlanets[6]*10**3      #Radius target planet [m]
disttarget = sys.a[6] * AU      #Distance to target planet [m]
periodtarget = sys.period[6]*60*60*24   #Rotational period target [s]

#Atmospere colculations:
rho0 = sys.rho0[6]                      #Surface density
T_0 = PlanetTemp(tempStar,radiusStar*10**3,disttarget) #Surface temperature
P_0 = rho0*k*T_0/(mu*m_H)               #Surface pressure
C = P_0**(1-gamma)*T_0**gamma           #Adiabatic constant
delta_r = 0.05                          #Steplength
h = np.arange(0,50e3,delta_r)           #Height above surface

def adiabaticP(P,height):
    r_center = height + radtarget
    g = gravacc(r_center)
    return -g*mu*m_H*P/(k*C**(1/gamma)*P**((gamma-1)/gamma))

def gravacc(r):
    g = G*masstarget/r**2
    return g

def isothermP(P,T,height):
    r_center = height + radtarget
    g = gravacc(r_center)
    return (-g*mu*m_H*P)/(k*T)

def ForwardEuler(P, P0,h,T0):
    """Solve P'(h)=const*P(h), T = (C/(P(h))**(1-gamma))**(1/gamma) untill
    isoterm part of Atmospere.
    P0: P(0)
    h: array for height above surface
    T0: T(P(0))
    """
    index_isotherm = 0
    n = len(h)
    u = np.zeros(n)   #u[k] is solution at height h[k]
    T = np.zeros(n)   #Tempereture, T[k] temp at height h[k]
    T[0] = T0
    u[0] = P0
    for k in range(n-1):
        if T[k] > T0/float(2):
            T[k+1] = (C/u[k]**(1-gamma))**(1/gamma)
            u[k+1] = u[k] + delta_r*adiabaticP(u[k], h[k])
        else:
            if index_isotherm == 0:
                index_isotherm = k
            T[k+1] = T[k]
            u[k+1] = u[k] + delta_r*isothermP(u[k],T[k],h[k])
    return u,T,index_isotherm

def plotting():
    plt.figure(1)
    plt.subplot(3,1,1)
    plt.plot(h,pressure,'b',label='pressure')
    plt.ylabel('P [N/m^2]')
    plt.plot([h[index_isotherm],h[index_isotherm]],[0,np.max(pressure)],'k')
    plt.legend(loc='best')
    plt.subplot(3,1,2)
    plt.plot(h,T,'r',label='Temperature')
    plt.ylabel('Temperature [K]')
    plt.plot([h[index_isotherm],h[index_isotherm]],[0,np.max(T)],'k')
    plt.legend(loc='best')
    plt.subplot(3,1,3)
    plt.plot(h,rho,label='Density')
    plt.ylabel('tetthet [kg/m^3]')
    plt.plot([h[index_isotherm],h[index_isotherm]],[0,np.max(rho)],'k')
    plt.xlabel('Height')
    plt.legend(loc='best')

if __name__ == '__main__':

    print 'surface temp',T_0
    print 'Surface pressure', P_0
    print 'surface density',rho0
    print 'Radius target', radtarget
    print 'rotational period',periodtarget
    pressure,T,index_isotherm = ForwardEuler(adiabaticP,P_0,h,T_0)
    rho = mu*m_H*pressure/(k*T)
    tol_rho = 0.001*np.max(rho)
    index_atmos = np.where(rho<tol_rho)[0] #Index for height when rho approx = 0
    min_dist = h[index_atmos[0]] + radtarget
    print 'Distance when rho approx equal 0',min_dist
    plotting()
    print 'Done'
    plt.show()
