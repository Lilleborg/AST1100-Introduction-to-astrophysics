from numpy import *
from AST1100SolarSystem import AST1100SolarSystem

seed = 4252
random.seed(seed)

#Constants:
solar_mass = 1.9891 * 10**30
m_H2 = 3.3436*10**(-27)         #Particle mass (H_2)
G = 6.6742 * 10**(-11)
m_satelite = 1100
m_lander = 90
T = 10000                      #temp gas
k = 1.38 * 10**(-23)            #Boltzmanns k
tot_fuel_loss = 1309.82439419
esc_part = 8163                #esc_part in period from box
dp_hole = 5.51748915773e-19     #momentum loss in period from box
nr_boxes = 4.49650353515e+12


#AST1100SolarSystem info:
mys = AST1100SolarSystem(seed)
M_0 = mys.mass[0]      #mass home planet [Solar masses]
R_0 = mys.radius[0]    #Radius home planet [km]

#Itteration constants:
t_box = 10**(-10)   #period from particle box
time = 1200.0         #20 min
n = 1200
dt = time/n
scale = dt/t_box    #Convert values from partclbox to new time steps

#Calculations:
v_esc = sqrt(2*G*M_0*solar_mass/R_0*10**3)#Esc velocity from home planet [m/s]

full_mass = tot_fuel_loss + m_satelite

delta_v = 0
def boost():
    ti = 1200.0         #20 min
    n_ = 12000
    dt_ = ti/n_
    scale = dt_/time    #Convert values from partclbox to new time steps
    full_mass = tot_fuel_loss + m_satelite
    delta_v = np.zeros(n_)
    for i in range(n_-1):
        delta_v[i+1] = delta_v[i]+dp_hole*scale*nr_boxes / full_mass
        full_mass = full_mass - esc_part*m_H2*scale*nr_boxes
    full_mass = full_mass - esc_part*m_H2*scale*nr_boxes
    print 'Velocity reached in 20 min:',delta_v[-1]
    print 'Escape velocity:',v_esc
    print 'Diff velocity reached and escape velocity:',abs(delta_v[-1] - v_esc)
    print 'Total mass after boosting 20 min:',full_mass
    print 'Number of boxes:', nr_boxes
    return delta_v, dt_
for i in range(n):
    delta_v += dp_hole*scale*nr_boxes / full_mass
    full_mass = full_mass - esc_part*m_H2*scale*nr_boxes
print v_esc
print n
print scale
print delta_v
print full_mass
