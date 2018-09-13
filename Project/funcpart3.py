import numpy as np
from AST1100SolarSystem import AST1100SolarSystem

seed = 4252
np.random.seed(seed)
sys = AST1100SolarSystem(seed)

#Constants
AU = 149597871. #[km]
G = 4 * np.pi**2        #Gravitation konstant when astro units
sek_yr = 365.*24*60*60   #second in year
k = 1.38 * 10**(-23)            #Boltzmanns k
solar_mass = 1.9891 * 10**30    #[kg]
h = 6.6260693*10**(-34)         #Plancks constant [Js]
c = 299792458.  #Speed of light [m/s]
nano = 10**(-9)   #nano
pi = np.pi

#AST1100SolarSystem info
starTemp = sys.temperature  #[K]
starRad = sys.starRadius    #[km]
starMass = sys.starMass     #[Solar masses]
dist_to_6 = sys.a[6] * AU * 10**3        #Distance to target planet [m]
rad_6 = sys.radius[6] * 10**3       #Radius target planet [m]

#Arrays:
wavelen = np.linspace(40,3000,10**4) #[nm]

#Functions:
def Blackrad(lamda,temp):
    """
    Plancks radiation law for black body.
    lamda: Array of wavelengths
    temp: Surface temperature of black body
    Returns intensity of different wavelengths
    """
    const = 2*h*c**2
    T = temp
    return const/(lamda*nano)**5 * 1/(np.exp(h*c/(k*T*lamda*nano))-1)

def Wiens(temp):
    """
    Wiens forskyvningslov.
    temp: Surface temp of black body
    Returns the most energetic wavelength
    """
    lamdamax = 2.9*10**(-3)/temp
    return lamdamax

def FluxSB(temp):
    """
    Stefan-Boltzmann law.
    temp: Surface temperature black body
    Returns flux emitted from black body
    """
    zigma = 2*pi**5 * k**4/(15*h**3*c**2)
    return zigma * temp**4

def Fluxrecieved(L,r,R=1):
    """
    Energy per time per area. Optional: Energy/s on absorbtion(shadow)-area
    L: Luminosity from emitter (star / black body object)
    r: distance from emitter
    R: Radius planet, shadow radius
    Returns flux distance r from star, Optional: Energy/s on shadow-area
    """
    area = 4*pi*r**2
    if R!=1:
        shadow = 4*pi*R**2
    else:
        shadow = 1
    return L/area * shadow

def PlanetTemp(temp,starradii,r):
    """
    Estimate surface temp on planet.
    temp: Surface temp star
    starradii: Radius star
    r: distance from star to planet
    Returns surface temperature on planet
    """
    Flux_from_planet = (Luminosity(FluxSB(starTemp),starradii)/(16*pi*r**2))
    zigma = 2*pi**5 * k**4/(15*h**3*c**2)
    T = (Flux_from_planet/zigma)**(1/4.)
    return T

def AreaPanels(F):
    """
    Area of solar panels on lander.
    F: Flux recieved on surface of planet
    Returns required area of panels, to produce 40W at 0.12 efficiency
    """
    A = 40./(F*0.12)
    return A

def Luminosity(flux,radius):
    """
    Total energy per time.
    flux: Flux from body
    radius: Radius of body
    Return Luminosity
    """
    area = 4*pi*radius**2
    return flux*area

def plotting():
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(wavelen,Blackrad(wavelen,starTemp))
    plt.plot(Wiens(starTemp)/nano,Blackrad(Wiens(starTemp)/nano,starTemp),'o',
    label='$\lambda$ most energetic from Wiens')
    plt.xlabel('$\lambda$ [nm]')
    plt.ylabel('intensity [W/m2/sr/Hz]')
    plt.legend()
    plt.show()

if __name__ == '__main__':

    #Prints:
    print 'Wavelength of most energy from star:',Wiens(starTemp)
    print 'Flux from surface of star', FluxSB(starTemp)
    print 'Luminosity star', Luminosity(FluxSB(starTemp),starRad*10**3)
    plotting()
    print 'Surface temperature planet', PlanetTemp(starTemp,starRad*10**3,dist_to_6)
    print 'Area solar panels on lander', AreaPanels(
    Fluxrecieved(Luminosity(FluxSB(starTemp),starRad*10**3),dist_to_6))
