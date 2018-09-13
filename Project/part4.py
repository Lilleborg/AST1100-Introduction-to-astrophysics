import numpy as np
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
import matplotlib.pyplot as plt
from PIL import Image

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

inFile = open('himmelkule.npy','rb')
himmelkulen = np.load(inFile)
inFile.close()

img = Image.open('sample0000.png')
pixels = np.array(img)
print np.shape(pixels)
width = len(pixels[0,:])
hight = len(pixels[1,:])
print width,hight

#Constants
AU = 149597871. #[km]
G = 4 * np.pi**2

#AST1100SolarSystem info
sys = AST1100SolarSystem(seed,hasMoons = False)
massStar = sys.starMass     #[Solar masses]
radiusStar = sys.starRadius #[km]
tempStar = sys.temperature  #[K]

N = sys.numberOfPlanets     #Number of planets
massPlanets = sys.mass      #[Solar masses]
radiusPlanets = sys.radius  #[km]
majoraxis = sys.a           #major axis [AU]
