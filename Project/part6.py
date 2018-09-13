from AST1100ShortcutSystem import AST1100SolarSystem
import numpy as np
from funcpart3 import PlanetTemp, FluxSB, Luminosity
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as m

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
min_dist = 3861734.77276

#Functions:


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

def lowerorbitboost(vvec,v,r,factor,nr_boost,t):
    v_dir = vvec/v
    if nr_boost in np.arange(2,100,2):   #Correction of v to maintain circle
        new_v = new_v = np.sqrt(G*masstarget/r)
        new_vvec = np.sqrt(G*masstarget/r)*v_dir
        omega = new_v/r
        period = 2*pi/omega

        delta_v = new_vvec - vvec
        vvec = new_vvec
        print 'Cor Size delta_v',np.linalg.norm(delta_v)
        outfile.write('boost %.2f %.6e %.6e %.6e'%(t-0.1,delta_v[0],
                        delta_v[1],delta_v[2]))
        outfile.write('\n')
        return vvec,period
    else:
        delta_v = -factor*v_dir   #change magnitude of velocity by factor
        vvec = vvec + delta_v
        outfile.write('boost %.2f %.6e %.6e %.6e'%(t-0.1,delta_v[0],
                        delta_v[1],delta_v[2]))
        outfile.write('\n')
        if nr_boost in np.arange(3,100,8):
            outfile.write('orient %.2f'%(t))
            outfile.write('\n')
            print 'Orient at t:',t,rvec[:,i]
            dir_pict = -rvec[:,i]
            theta = (m.atan2(dir_pict[2],
                    m.sqrt(dir_pict[0]**2+dir_pict[1]**2))) + pi/2.
            phi =  m.atan2(dir_pict[1],dir_pict[0])
            outfile.write('picture %.2f %.3f %.3f %.3f %.3f %.3f'%(t+1,theta,phi,
                rvec[0,i],rvec[1,i],rvec[2,i]))
            outfile.write('\n')
        print 'Red Size delta_v',np.linalg.norm(delta_v)
        return vvec

def directionboost(vvec,v,rvec,dv,t,nr_dir):
    delta_dir = np.cross(rvec,vvec) #Direction of boost with magnitude
    delta_dir = delta_dir/np.linalg.norm(delta_dir)#Unit vector direction boost
    vvec_new = vvec + dv*delta_dir
    v_new = np.linalg.norm(vvec_new)
    vdir_new = vvec_new/v_new
    vvec_new = vdir_new*v
    delta_v = vvec_new - vvec
    vvec = vvec_new
    print 'Dir Size delta_v',np.linalg.norm(delta_v)
    outfile.write('boost %.2f %.6e %.6e %.6e'%(t-0.1,delta_v[0]
                    ,delta_v[1],delta_v[2]))
    outfile.write('\n')
    if nr_dir in np.arange(2,100,4):
        outfile.write('orient %.2f'%(t+0.1))
        outfile.write('\n')
        print 'Orient at t:',t,rvec
    return vvec

def plotting():
    plt.figure(3)
    plt.subplot(2,1,1)
    plt.plot(t,r,label='Absolute distance to senter')
    plt.plot(t,rvec[0,:],'--',label='Distance to star x-comp')
    plt.plot(t,rvec[2,:],label='Distance to star z-comp')
    plt.ylabel('Distance [m]')
    plt.legend(loc='best')
    plt.subplot(2,1,2)
    plt.plot(t,v,label='Absolute velocity')
    plt.ylim(0,np.max(v)*2)
    plt.legend(loc='best')
    plt.ylabel('Velocity [m/s]')
    plt.xlabel('Time [s]')

#Orbit simulations:
time = 1.65e6 #Time in orbit [s]
dt = 1.             #Time step
n = time/dt
n = int(n)          #Number of time steps
rvec = np.zeros((3,n+1))
vvec = np.zeros((3,n+1))
r = np.zeros(n+1)       #Absolute distance to senter
v = np.zeros(n+1)       #Absolute velocity
a = np.zeros(3)
t = np.zeros(n+1)       #Time in orbit
t[0] = 0
rvec[:,0] = [204390971.407 ,  0.0 ,  0 ]   #Initial relative posistion [m]
vvec = [0.00463126267347 ,  691.395286682 ,  0]   #Initial relative vel [m/s]
r[0] = np.linalg.norm(rvec[:,0])
v[0] = np.linalg.norm(vvec)
nr_boost = 1
N_l_0 = 22
nr_plane = 36
N_l_1 = N_l_0+nr_plane+4

index_time_dir_boosts = []      #Used to find
pos_dir_boost = np.zeros((3,nr_plane-1))


outfile = open('DeltaV.txt','w')
for i in range(n):
    a = -G*masstarget*rvec[:,i]/r[i]**3
    t[i+1] = t[i] + dt
    #Reduce radii of orbit:
    if i == 1:  #Initial boost
        vvec = lowerorbitboost(vvec,v[i],r[i],200,nr_boost,t[i])
        nr_boost += 1
    if i == 2:
        print 'Orient at t:',t[i],rvec[:,i]
    #Circular correction boosts:
    if rvec[0,i] > rvec[0,i-2] and rvec[1,i]<rvec[1,i-1]:
        if nr_boost in np.arange(2,N_l_0+1,2):
            vvec,period = lowerorbitboost(vvec,v[i],r[i],30,nr_boost,t[i])
            index_time_dir_boosts.append(i)
            nr_boost += 1
    #Reduce radii of orbit
    if rvec[0,i]<rvec[0,i-1]:
        if nr_boost in np.arange(3,N_l_0,2):
            vvec = lowerorbitboost(vvec,v[i],r[i],150,nr_boost,t[i])
            nr_boost += 1
    #Change plane of orbit
    if rvec[0,i]<rvec[0,i-1]:
        if nr_boost in np.arange(N_l_0,N_l_0 + nr_plane,1):
            if (nr_boost in np.arange(2,100,8) and
                    i in np.arange(index_time_dir_boosts[-1],n,1043)):
                outfile.write('orient %.2f'%(t[i]))
                outfile.write('\n')
                print 'Orient at t:',t[i],rvec[:,i]
                dir_pict = -rvec[:,i]
                theta = (m.atan2(dir_pict[2],
                        m.sqrt(dir_pict[0]**2+dir_pict[1]**2))) + pi/2.
                phi = m.atan2(dir_pict[1],dir_pict[0])
                outfile.write('picture %.2f %.3f %.3f %.3f %.3f %.3f'%(t[i]+1,theta,phi,
                    rvec[0,i],rvec[1,i],rvec[2,i]))
                outfile.write('\n')
            if i in np.arange(index_time_dir_boosts[-1],n,int(period/dt)):
                vvec = directionboost(vvec,v[i],rvec[:,i],202,t[i],nr_boost)
                nr_boost += 1
                pos_dir_boost[:,(nr_boost-(N_l_0+2))] = rvec[:,i]
                nr_boosts_untill_now = nr_boost

    #Circular correction boosts:
    if rvec[0,i] > rvec[0,i-1]:
        if nr_boost in np.arange(N_l_0+nr_plane,N_l_1+1,2):
            vvec,period = lowerorbitboost(vvec,v[i],r[i],30,nr_boost,t[i])
            index_time_dir_boosts.append(i)
            nr_boost += 1
            if nr_boost >62:
                outfile.write('orient %.2f'%(t[i]))
                outfile.write('\n')
                print 'Orient at t:',t[i],rvec[:,i]

    #Reduce radii of orbit
    if rvec[0,i]<rvec[0,i-1]:
        if nr_boost in np.arange(N_l_0+nr_plane+1,N_l_1,2):
            vvec = lowerorbitboost(vvec,v[i],r[i],150,nr_boost,t[i])
            nr_boost += 1

    vvec = vvec + a*dt
    rvec[:,i+1] = rvec[:,i] + vvec*dt
    r[i+1] = np.linalg.norm(rvec[:,i+1])
    v[i+1] = np.linalg.norm(vvec)
outfile.close()
print 'Distance to minimum distance', np.min(r)-min_dist
print 'Distance to center min', min_dist
print 'Distance to center',np.min(r),r[-1]


#slicing:
rvec = rvec[:,0::30]
r = r[0::30]
v = v[0::30]
t = t[0::30]

fig = plt.figure(2)
ax = fig.gca(projection='3d')
ax.plot(rvec[0,:],rvec[1,:],rvec[2,:])
ax.plot([rvec[0,0]],[rvec[1,0]],[rvec[2,0]],'yo')
#for k in range(nr_plane-1): #nr_boosts_untill_now-12):
#    ax.plot([pos_dir_boost[0,k]],[pos_dir_boost[1,k]],[pos_dir_boost[2,k]],'ro')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')


#Sphere for planet
u = np.linspace(0, 2 * np.pi, 100)
w = np.linspace(0, np.pi, 100)
x = radtarget * np.outer(np.cos(u), np.sin(w))
y = radtarget * np.outer(np.sin(u), np.sin(w))
z = radtarget * np.outer(np.ones(np.size(u)), np.cos(w))
ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b')
axlim = np.max(r)
ax.auto_scale_xyz([-3e8,3e8],[-3e8, 3e8],[-3e8, 3e8])
plotting()
print 'done'
plt.show()
