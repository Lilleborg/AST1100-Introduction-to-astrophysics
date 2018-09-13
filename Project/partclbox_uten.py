import numpy as np
from AST1100SolarSystem import AST1100SolarSystem
from AST1100SolarSystemViewer import AST1100SolarSystemViewer



seed = 4252
np.random.seed(seed)

#Functions:
def plot(pos,N,n):
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as p3
    import matplotlib.animation as animation
    """
    Plots the simulation of the particle Moving
    """
    def update_lines(num, dataLines,lines):
        for line, data in zip(lines,dataLines):
            line.set_data(data[0:2, num-1:num])
            line.set_3d_properties(data[2,num-1:num])
        return lines

    fig = plt.figure()
    ax = p3.Axes3D(fig)

    m = 100 #number of frames
    n = N  #number of particles
    N = n    #number of time steps
    #pos = np.reshape(pos,(n,3,N))
    data = pos
    #np.zeros([n,3,N])

    lines = [i for i in range(n)]
    for i in range(n):
        lines[i] = [ax.plot(data[i][0,0:1],
        data[i][1,0:1],data[i][2,0:1],'o')[0]]

    ax.set_xlim3d([0.0,L])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0,L])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0,L])
    ax.set_zlabel('Z')

    ani = [i for i in range(n)]
    for i in range(n):
        ani[i] = animation.FuncAnimation(fig,update_lines,m,
        fargs=([data[i]],lines[i]),interval=50,blit=False)
    plt.show()

def Kinetic_E(v,m=3.3436*10**(-27)):
    """
    returns kinetic energy for particle with velocity vector v and mass m.
    """
    return 1./2 * m * (np.linalg.norm(v,axis=1))**2
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

#Constants:
N = 10**5                    #nr partc [10**5]
L = 10.**(-6)                   #length box [10**[-6]]
T = 10**4                      #temp gas    [10**4]
k = 1.38 * 10**(-23)            #Boltzmanns k
m_H2 = 3.3436*10**(-27)         #Particle mass (H_2)
sigma_v = np.sqrt(k*T/m_H2)           #Stdr deviation velocity
G = 6.67 * 10**(-11)
m_satelite = 1100
m_lander = 90
solar_mass = 1.9891 * 10**30

#Analytical results:
mean_E_K = 3./2 * k * T           #Theroeticly correct mean kinetic energy
mean_absv = np.sqrt(8*k*T/(np.pi*m_H2))      #Analytic result from 1A.5
P = k * T * (N/L**3)

#Itteration Constants:
time = 10.0**(-9)              #time period
n = 1000                       #nr time steps
dt = time/n                  #time step interval

#AST1100SolarSystem info:
mys = AST1100SolarSystem(seed)
M_0 = mys.mass[0]      #mass home planet [Solar masses]
R_0 = mys.radius[0]    #Radius home planet [km]
mysViewer = AST1100SolarSystemViewer(seed)

#Arrays:
#pos = np.zeros((N,3,n+1)) #position array with time steps for plotting(testing)
L_tol = L/10**6         #Min length toleranse from side of the box to:
                        #prevent the particles to start on the edge
#Sums:
nr_ofhits_z = 0    #nr of hits on z = 0
sum_dp_z = 0        #momentum particles hitting z = 0
esc_part = 0        #particles escaping hole
dp_hole = 0         #momentum through hole

rvec = np.random.uniform(L_tol,L - L_tol,(N,3)) #initial position
vvec = np.random.normal(0,sigma_v,(N,3))    #initial velocity
mean_E_K_calc = np.sum(Kinetic_E(vvec))/N  #mean initial kinetic energy
mean_abs_v_calc = np.sum(np.linalg.norm(vvec,axis=1))/N #mean initial abs vel
#pos[:,:,0] = rvec   #position with time steps for plotting
for j in range(n-1):    #Itterating time steps
    #pressure
    indexpartHitZ = np.nonzero(np.where(rvec[:,2]<=L_tol,1,0)) #Index partcl z=0
    dp_z= 2*m_H2*np.sum(np.absolute(vvec[indexpartHitZ,2]))#momentum left on z=0
    nr_ofhits_z += np.sum(np.where(rvec[:,2]<=L_tol,1,0))
    sum_dp_z += np.sum(dp_z)
    #Escaped momentum
    #Particle hit hole:
    xLow = np.where(rvec[:,0]>L/4,1,0)  #lower x condition
    xHigh = np.where(rvec[:,0]<L-L/4,1,0)   #higher x condition
    yLow = np.where(rvec[:,1]>L/4,1,0)      #lower y condition
    yHigh = np.where(rvec[:,1]<L-L/4,1,0)   #higher y condition
    z = np.where(rvec[:,2]<=L_tol,1,0)      #z condition
    hit_hole = xLow*xHigh*yLow*yHigh*z      #particles hit hole, true=1
    indexHitHole = np.nonzero(hit_hole)     #index for particles hit hole

    dp_hole += np.sum(2 * m_H2 * np.absolute(vvec[indexHitHole,2]))
    esc_part += np.sum(hit_hole)

    #moving escaped particles to the top and middle of box:
    rvec[indexHitHole,0:1] = L/2
    rvec[indexHitHole,2] = L-L_tol

    #Colliding particles:
    for i in range(3):  #itterating x,y,z
        collidMin = np.where(rvec[:,i]<L_tol,1,0)
        collidMax = np.where(rvec[:,i]>L-L_tol,1,0)
        partColliding = collidMax+collidMin
        indexColliding = np.nonzero(partColliding)
        vvec[indexColliding,i] = -vvec[indexColliding,i]
        #moving particles inside box:
        rvec[np.nonzero(collidMin),i] = L_tol
        rvec[np.nonzero(collidMax),i] = L-L_tol
    #pos[:,:,j] = rvec #only included during testing
    rvec = rvec + vvec*dt
#plot(pos,N,n)

#Calculations:
v_esc = np.sqrt(2*G*M_0*solar_mass/(R_0*10**3))#Esc velocity from home planet [m/s]
mean_F_z = sum_dp_z/(time)            #mean force on wall
P_calc = mean_F_z / L**2          #calculated pressure

#Speed gain and fuel loss:
dv = dp_hole / (m_satelite)  #Speed gain from momentum loss in time
nr_boxes = v_esc  /(dv*20*60/time) #Nr box to reach v_esc within 20 min
esc_part_in_20 = esc_part * (20*60/time)
tot_fuel_loss=esc_part_in_20*nr_boxes*m_H2 #mass all particles escaping in 20 min


#Difference and relative error calulations:
diff_mean_abs_v = abs(mean_absv - mean_abs_v_calc)
rel_error_abs_v = mean_abs_v_calc / mean_absv
diff_mean_E_K = abs(mean_E_K_calc - mean_E_K)
rel_error_E_K = mean_E_K_calc / mean_E_K
diff_P = abs(P_calc - P)
rel_error_P = P_calc / P
delta_v , dt_ = boost()
print '---'
print 'Number particles, time steps', N,n
print 'Time', time
print 'Length box, temperature',L,T
print 'Fuel required', tot_fuel_loss
#Printing difference in calculated and analytic results
"""
print 'Momentum loss in time', dp_hole
print 'Escaped particles', esc_part
print 'Speed gain from momentum loss', dv
print 'Boxes required', nr_boxes
print 'Escaped part in 20 min', esc_part_in_20
print 'Mass,radius', M_0, R_0
"""
print '---'

print ("Difference in calculated and analytic mean kinetic energy: \n"
        "%g \n"
        "Relative difference: \n"
        "%g" %(diff_mean_E_K,rel_error_E_K))
print ("Difference in calculated and analytic mean absolute velocity: \n"
        "%g \n"
        "Relative difference: \n"
        "%g" %(diff_mean_abs_v,rel_error_abs_v))
print ("Difference in calculated and analytic pressure: \n"
        "%g \n"
        "Relative difference: \n"
        "%g" %(diff_P,rel_error_P))

mys.massNeededCheck(nr_boxes,v_esc,dp_hole/time,esc_part/time,1100)
#mysViewer.escapeVelMovie(delta_v, dt_)
