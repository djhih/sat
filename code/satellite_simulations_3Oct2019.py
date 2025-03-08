# modify from Spooky action at a global distance: analysis of space-based entanglement distribution for the quantum internet
# src link https://arxiv.org/src/1912.06678v3/anc/satellite_simulations_3Oct2019.py
 
import numpy as np
import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib import cm
import dill
import pickle
from scipy.linalg import expm
from numpy.linalg import norm
from scipy.optimize import brentq,brenth
from datetime import date
from joblib import Parallel,delayed
import time
from scipy.optimize import minimize,minimize_scalar,fminbound,root_scalar
import warnings


def get_constants(h):

    # Last modified: 13 August 2019
    '''
    h is the altitude of the satellites in meters
    '''

    G=6.67408e-11 # Gravitational constant in SI units
    R=6.371e6 # Radius of the Earth in meters
    M=5.972e24 # Mass of the earth in kg
    Omega=2*np.pi/86400.  # Rotation speed of the Earth in radians/s
    t=10e3 # Thickness of the atmosphere in meters

    #Earth_axis=[0,0,1]  # Axis of rotation of the earth: taken to be along the z-axis.

    #h=1000e3 # Altitude of the satellites in meters
    omega=np.sqrt(G*M)/(R+h)**(3./2.)  # Rotation speed of the satellites in radians/s

    #period=2*np.pi/omega # Period of each satellite in seconds.

    return G,R,M,Omega,t,h,omega


def eta_sg(L,lam=810e-9,w0=0.025,rg=0.75):

    # Last modified: 18 July 2019

    '''
    Satellite-to-ground transmittance. All distance quantities are
    in meters.
    '''

    LR=np.pi*w0**2/lam

    wg=w0*np.sqrt(1+L**2/LR**2)

    return 1-np.exp(-2*rg**2/wg**2)


def eta_Atm(L,h,etaZen=0.5):

    # Last modified: 18 July 2019

    '''
    Atmospheric transmittance. L is the link distance (distance from
    satellite to ground station) in meters, h is the altitude of the 
    satellites (in meters).
    '''

    warnings.filterwarnings("ignore")

    _,R,_,_,_,_,_=get_constants(h)

    secZen=1/((h/L)-(L**2-h**2)/(2*R*L))

    return etaZen**secZen


def eta_Tot(L,h):

    # Last modified: 18 July 2019

    '''
    Total transmittance, which is a product of the atmospheric transmittance
    and the beam spreading transmittance. L is the link distance (distance between
    the ground station and the satellite) in meters, and h is the altitude of the
    satellite (in meters).
    '''

    return eta_sg(L)*eta_Atm(L,h)


def link_distance(d,h):

    # Last modified: 18 July 2019

    '''
    Returns the link distance L between two ground stations separated by an 
    arc-distance of d (in meters) and a satellite at an altitude of
    h (in meters) located at the midpoint of the ground stations (so both
    ground stations are the same distance away from the satellite).
    '''

    _,R,_,_,_,_,_=get_constants(h)

    return np.sqrt(4*R*(R+h)*np.sin(d/(4*R))**2+h**2)


def opt_alt_arc(d):

    # Last modified: 18 July 2019

    '''
    Finds the optimal satellite altitude for two ground stations separated
    by an arc-distance of d (in meters). Optimality is in terms of the 
    transmissivity eta_Tot.
    '''

    def objfunc(h):

        return -eta_Tot(link_distance(d,h),h)

    
    opt=fminbound(objfunc,0,1500000)
    #opt=minimize_scalar(objfunc,bounds=(0,1500000),method='bounded')

    return opt


def generate_link_constraint(d):

    # Last modified: 18 July 2019

    '''
    Given a maximum allowed arc-distance separation d (in meters) between
    two ground stations, this function generates the maximum allowed distance
    between a ground station and a satellite based on the optimal altiude
    of the satellite.
    '''

    h_opt=opt_alt_arc(d)

    return link_distance(d,h_opt)


def opt_arcsep(h,eta_star,dB=False,max_d=2e6):

    # Last modified: 18 July 2019

    '''
    This function finds the largest arc-distance separation d (in meters) 
    between two ground stations that is possible given a fixed altitude h
    (in meters) of the satellite and a given threshold eta_star for the
    transmissivity. The threshold eta_star is the smallest tolerable 
    transmissivity (i.e., the transmissivity should be larger than eta_star).

    If dB=True, then eta_star is the loss, specified in dB. The loss and
    transmissivity are related via:

    loss=-10*log_10(transmissivity)
    '''


    def objfunc(d):

        return eta_star-eta_Tot(link_distance(d,h),h)


    if dB:
        eta_star=10**(-eta_star/10)

    r=root_scalar(objfunc,bracket=[0,max_d])

    return r.root


def generate_rotation(axis,angle):
    
    # Last modified: 13 March 2019
    
    '''
    Generates the rotation matrix about the given axis (unit vector in Cartesian coordinates) by the given angle.
    '''
    
    # Generators of the SO(3) rotations
    Lx=np.matrix([[0,0,0],[0,0,-1],[0,1,0]])
    Ly=np.matrix([[0,0,1],[0,0,0],[-1,0,0]])
    Lz=np.matrix([[0,-1,0],[1,0,0],[0,0,0]])
    
    n1=axis[0]
    n2=axis[1]
    n3=axis[2]
    
    L=n1*Lx+n2*Ly+n3*Lz
    
    return np.matrix(expm(angle*L))


def Rx(t):
    
    # Last modified: 13 March 2019
    
    '''
    Rotation matrix about the x-axis by angle t.
    '''
    
    return generate_rotation([1,0,0],t)
    #return np.matrix([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
    
    
def Ry(t):
    
    # Last modified: 13 March 2019
    
    '''
    Rotation matrix about the y-axis by angle t.
    '''
    
    return generate_rotation([0,1,0],t)
    #return np.matrix([[np.cos(t),0,np.sin(t)],[0,1,0],[-np.sin(t),0,np.cos(t)]])
    
    
def Rz(t):
    
    # Last modified: 13 March 2019
    
    '''
    Rotation matrix about the z-axis by angle t.
    '''
    
    return generate_rotation([0,0,1],t)
    #return np.matrix([[np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0],[0,0,1]])


def generate_point(R,theta,phi):
    
    '''
    Generates (in Cartesian coordinates) the point in 3D specified by
    the spherical coordinates R, theta, phi
    
    Here, theta is the azimuthal angle (positive going counter-clockwise from the positive x-axis), and
    phi is the polar angle (positive going down from the z-axis)
    '''
    
    return [R*np.sin(phi)*np.cos(theta),R*np.sin(phi)*np.sin(theta),R*np.cos(phi)]
    
    
def ground_distance(g1,g2,R=6.371e6,spherical=False):

    # Last modified: 15 August 2019

    '''
    Given two points g1 and g2 on the Earth (in Cartesian coordinates by default),
    this function finds the straight line ground distance between the points based
    on the angular separation of the points.
    '''

    if spherical:  # If points are given in spherical coordinates, convert them to Cartesian first
        g1=generate_point(g1[0],g1[1],g1[2])
        g2=generate_point(g2[0],g2[1],g2[2])

    g1=np.array(g1)
    g2=np.array(g2)

    d=norm(g1-g2)  # Euclidean distance between g1 and g2

    
    D=R*np.arcsin((d/(2*R**2))*np.sqrt(4*R**2-d**2))

    return D



def evolve_vector(init_point,axis,omega,T):
    
    '''
    Performs rotation of a given initial point (init_point, in Cartesian coordinates) by the given axis 
    (unit vector in Cartesian coordinates) at a given rotation speed (omega) for a time T (starting at t=0).
    '''
    
    
    init_point=np.matrix(init_point).H
    
    R=generate_rotation(axis,omega*T)
    point_t=R*init_point

    point_t_x=np.array(point_t[0])[0][0]
    point_t_y=np.array(point_t[1])[0][0]
    point_t_z=np.array(point_t[2])[0][0]
    
    return [point_t_x,point_t_y,point_t_z]


def generate_network(num_rings,num_sats,R,h,offset_angle=0,multiple_alt=False):
    
    # Last modified: 3 October 2019
    
    '''
    Generates the satellite and ground station architecture consisting of num_rings equally-spaced 
    longitudinal rings of satellites and ground stations, and num_sats equally-spaced 
    satellites and ground stations in each ring.

    R is the radius of the earth (in meters), and h is the altitude of the
    satellites (in meters)

    If multiple_alt=True, then h should be a list of altitudes.
    '''
    
    S={}
    G={}

    count1=0

    for m in range(num_rings):
        count2=0
        count1+=1
        S[count1]={}
        G[count1]={}
        
        for k in range(1,num_sats+1):
            
            count2+=1

            if multiple_alt==True:
                S[count1][count2]={}
                for alt in h:  # h should be a list of altitudes (in meters)
                    S[count1][count2][alt/1000]=generate_point(R+alt,m*np.pi/num_rings,k*2*np.pi/num_sats+offset_angle)
            else:
                S[count1][count2]=generate_point(R+h,m*np.pi/num_rings,k*2*np.pi/num_sats+offset_angle)
            
            
            
            G[count1][count2]=generate_point(R,m*np.pi/num_rings+np.pi/(2*num_rings),k*2*np.pi/num_sats+offset_angle)
        
        S[count1]['axis']=[-np.sin(m*np.pi/num_rings),np.cos(m*np.pi/num_rings),0]
            
    
    return S,G


def plot_initial_configuration(G,S,R):

    # Last modified: 28 October 2019

    '''
    Plots the initial configuration of the ground stations and the satellites.

    R is the radius of the earth (in meters).
    '''

    fig=plt.figure(figsize=(8,8))


    ax=fig.add_subplot(111,projection='3d')

    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
    x = (R)*np.cos(u)*np.sin(v)
    y = (R)*np.sin(u)*np.sin(v)
    z = (R)*np.cos(v)
    graph_earth=ax.plot_wireframe(x, y, z, color="gray")

    for i in list(G.keys()):#range(len(G)):
        for j in list(G[i].keys()):#range(len(G[i])):
            ax.scatter(G[i][j][0],G[i][j][1],G[i][j][2],c='g')        

    if S!=None:
        for i in list(S.keys()):#range(len(S)):
            for j in list(S[i].keys()):#range(len(S[i])-1):
                if j=='axis':
                    continue
                else:
                    ax.scatter(S[i][j][0],S[i][j][1],S[i][j][2],c='r')  

    ax.set_xlabel('$x$ (km)')
    ax.set_ylabel('$y$ (km)')
    ax.set_zlabel('$z$ (km)')

    ax.set_xlim(-8000000,8000000)
    ax.set_ylim(-8000000,8000000)
    ax.set_zlim(-8000000,8000000)

    ax.set_xticklabels(['$-8000$','$-6000$','$-4000$','$-2000$','$0$','$2000$','$4000$','$6000$','$8000$'])
    ax.set_yticklabels(['$-8000$','$-6000$','$-4000$','$-2000$','$0$','$2000$','$4000$','$6000$','$8000$'])
    ax.set_zticklabels(['$-8000$','$-6000$','$-4000$','$-2000$','$0$','$2000$','$4000$','$6000$','$8000$'])
     
    
    plt.show()


def generate_ground_station_pairs(G):

    # Last modified: 23 March 2019

    '''
    Pairs the ground stations in G (as generated from the function generate_network)
    according to a nearest-neighbour square-grid topology.
    '''

    G_pairs=[]
    G_pair_labels=[]

    num_rings=len(G)
    #num_ground_ring=len(G[1])

    for m in list(G.keys()):#range(num_rings):
        num_ground_ring=len(G[m])
        for k in list(G[m].keys()):#range(num_ground_ring):
            if m==num_rings:
                gs1=list(G.keys())[0]
            else:
                gs1=m+1
            
            if k==num_ground_ring:
                gs2=list(G[m].keys())[0]
            else:
                gs2=k+1

            G_pairs.append((G[m][k],G[m][gs2]))
            G_pairs.append((G[m][k],G[gs1][k]))

            G_pair_labels.append(([m,k],[m,gs2]))
            G_pair_labels.append(([m,k],[gs1,k]))


    return G_pair_labels,G_pairs



def plot_ground_pairs(points,R):

    # Last modified: 23 March 2019

    '''
    Plots all of the given points (each point a list of coordinates).
    '''

    fig=plt.figure(figsize=(8,8))


    ax=fig.add_subplot(111,projection='3d')

    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
    x = (R)*np.cos(u)*np.sin(v)
    y = (R)*np.sin(u)*np.sin(v)
    z = (R)*np.cos(v)
    graph_earth=ax.plot_wireframe(x, y, z, color="b")


    for point in points:
        ax.scatter(point[0],point[1],point[2],c='g')


    #for i in range(len(pairs)):
    #    ax.scatter(pairs[i][0][0],pairs[i][0][1],pairs[i][0][2],c='g')
    #    ax.scatter(pairs[i][1][0],pairs[i][1][1],pairs[i][1][2],c='g')

    #for i in range(len(S)):
    #    for j in range(len(S[i])-1):
    #        ax.scatter(S[i][j][0],S[i][j][1],S[i][j][2],c='r')  

    ax.set_xlabel('$x$ (km)')
    ax.set_ylabel('$y$ (km)')
    ax.set_zlabel('$z$ (km)')

    ax.set_xlim(-8000000,8000000)
    ax.set_ylim(-8000000,8000000)
    ax.set_zlim(-8000000,8000000)

    ax.set_xticklabels(['$-8000$','$-6000$','$-4000$','$-2000$','$0$','$2000$','$4000$','$6000$','$8000$'])
    ax.set_yticklabels(['$-8000$','$-6000$','$-4000$','$-2000$','$0$','$2000$','$4000$','$6000$','$8000$'])
    ax.set_zticklabels(['$-8000$','$-6000$','$-4000$','$-2000$','$0$','$2000$','$4000$','$6000$','$8000$'])
     
    
    plt.show()



def evolve_network(S,G,G_pair_labels,T,eta_star,h,Omega,t_start=0,display=False,parallel=False,multiple_alt=False):
    
    # Last modified: 3 October 2019
    
    '''
    Simulates the network with satellites in S and ground stations in G for a total time of T.
    t_start is the starting time (default is t_start=0).
    
    eta_star is the threshold transmissivity between satellites and ground stations (i.e., the total transmissivity
    between both neighboring stations and the satellite has to be greater than or equal to eta_star).

    h is the altitude of the satellites. If multiple_alt=True, then h should be a list of
    altitudes of the satellites.
    '''
    
    Earth_axis=[0,0,1]
    
    
    t_array=np.linspace(t_start,t_start+T,T+1)
    
    Data={}
        
    def time_evolve(t):

        if display:
            print('Time',t,'of',T)
        
        Data[t]={}
        Data[t]['ground station pairs']={}

        # Iterate over the ground station pairs
        for i in range(len(G_pair_labels)):

            label=str(G_pair_labels[i])

            Data[t]['ground station pairs'][label]={}

            g1_label=G_pair_labels[i][0]
            g2_label=G_pair_labels[i][1]

            g1_pos=G[g1_label[0]][g1_label[1]]
            g2_pos=G[g2_label[0]][g2_label[1]]

            g1_pos_evolve=evolve_vector(g1_pos,Earth_axis,Omega,t)
            g2_pos_evolve=evolve_vector(g2_pos,Earth_axis,Omega,t)   

            Data[t]['ground station pairs'][label]['position']=(g1_pos_evolve,g2_pos_evolve)

            # List of satellites that are in range of the ground station pair. Only those that are in range are stored.
            if multiple_alt==True:
                Data[t]['ground station pairs'][label]['in range']={}
                for alt in h:
                    Data[t]['ground station pairs'][label]['in range'][alt/1000]=[]
            else:
                Data[t]['ground station pairs'][label]['in range']=[]
            

            # Iterate over the rings of satellites
            for k in list(S.keys()):
                
                ring_axis=S[k]['axis']
                
                # Iterate over the satellites in the ring
                for l in list(S[k].keys()):
                    if l=='axis':
                        continue
                    
                    if multiple_alt==True:
                        for alt in h:  # h should be a list of altitudes (in meters)
                            
                            _,_,_,_,_,_,omega=get_constants(alt)
                            
                            # Evolve the satellite
                            sat_pos=evolve_vector(S[k][l][alt/1000],ring_axis,omega,t)
                            
                            g1_diff=np.array(sat_pos)-np.array(g1_pos_evolve)
                            g2_diff=np.array(sat_pos)-np.array(g2_pos_evolve)

                            d1=norm(g1_diff)  # distance between the first ground station in the pair and the satellite
                            d2=norm(g2_diff)  # distance between the second ground station in the pair and the satellite

                            eta_tot=eta_Tot(d1,alt)*eta_Tot(d2,alt)  # The total transmissivity.

                            # Computes the cosine of the angle between the position vectors of the ground stations and the satellite relative to the ground stations.
                            # The "horizon condition" is that cos(theta)>=0, i.e., theta<=pi/2. We call this the horizon condition
                            # because the satellite should be in the horizon of both ground stations (i.e., physically visible)
                            cos_angle1=np.dot(g1_diff,g1_pos_evolve)/(norm(g1_diff)*norm(g1_pos_evolve))
                            cos_angle2=np.dot(g2_diff,g2_pos_evolve)/(norm(g2_diff)*norm(g2_pos_evolve))

                            if eta_tot>=eta_star and cos_angle1>=0 and cos_angle2>=0:
                                Data[t]['ground station pairs'][label]['in range'][alt/1000].append([(k,l),d1,d2,eta_tot,-10*np.log10(eta_tot)])

                    else:
                        _,_,_,_,_,_,omega=get_constants(h)

                        # Evolve the satellite
                        sat_pos=evolve_vector(S[k][l],ring_axis,omega,t)
                        
                        g1_diff=np.array(sat_pos)-np.array(g1_pos_evolve)
                        g2_diff=np.array(sat_pos)-np.array(g2_pos_evolve)

                        d1=norm(g1_diff)
                        d2=norm(g2_diff)

                        eta_tot=eta_Tot(d1,h)*eta_Tot(d2,h)

                        # Computes the cosine of the angle between the position vectors of the ground stations and the satellite relative to the ground stations.
                        # The "horizon condition" is that cos(theta)>=0, i.e., theta<=pi/2. We call this the horizon condition
                        # because the satellite should be in the horizon of both ground stations (i.e., physically visible)
                        cos_angle1=np.dot(g1_diff,g1_pos_evolve)/(norm(g1_diff)*norm(g1_pos_evolve))
                        cos_angle2=np.dot(g2_diff,g2_pos_evolve)/(norm(g2_diff)*norm(g2_pos_evolve))

                        if eta_tot>=eta_star and cos_angle1>=0 and cos_angle2>=0:
                            Data[t]['ground station pairs'][label]['in range'].append([(k,l),d1,d2,eta_tot,-10*np.log10(eta_tot)])

        return Data[t]

    if parallel:
        dat=Parallel(n_jobs=4)(delayed(time_evolve)(t) for t in t_array)
        return dat
    else:
        for t in t_array:
            time_evolve(t)      
                        
        return Data


def range_times(data,pair,h=None):
    
    # Last modified: 3 October 2019

    '''
    Grabs the information about when the given pair of ground stations is in
    the range of all satellites

    Pair should be specified as a string.

    h=None means that data only contains satellites with one altitude. If data
    contains satellites with multiple altitudes, then h should specify (in meters)
    the desired altitude.
    '''
    
    in_range=[] # Tells us whether the ground station pair is in the range of some satellite as a function of time.
    #count=0  # Gives the total number of time steps for which the given pair of ground stations is in the range of some satellite.
    times=[] # Gives the exact times that the ground station pair is in the range of some satellite.
    time_gaps=[]
    num_sats=[] # Number of satellites in the range of the ground station pair as a function of time.
    
    for t in range(len(data)):
        
        if h==None:
            range_list=data[t]['ground station pairs'][pair]['in range']
        else:
            range_list=data[t]['ground station pairs'][pair]['in range'][h/1000]

        if not range_list:  # If range_list is empty, then there is no satellite in range at the particular time
            num_sats.append(0)
            time_gaps.append(t)
            in_range.append(0)
        else:
            num_sats.append(len(range_list))
            times.append(t)
            in_range.append(1)
    
    return np.array(in_range),np.array(times),np.array(time_gaps),np.array(num_sats)

def haversine_distance(coord1, coord2, R=6371000):
    """
    計算地球表面兩點之間的最短大圓距離（單位：公尺）。
    
    參數:
        coord1: (lat1, lon1) - 第一個點的經緯度（單位：度）
        coord2: (lat2, lon2) - 第二個點的經緯度（單位：度）
        R: 地球半徑（預設值：6371000 公尺）
    
    回傳:
        兩點之間的距離（公尺）
    """
    # 將經緯度從度轉換成弧度
    lat1, lon1 = np.radians(coord1)
    lat2, lon2 = np.radians(coord2)
    
    # 計算經緯度差異
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    # haversine公式
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    distance = R * c
    return distance

def ecef_to_geodetic(x, y, z):
    """
    將 ECEF 座標轉換為地理座標（緯度、經度，單位：弧度）
    假設地球為完美球體。
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.arcsin(z / r)           # 或用 np.arctan2(z, np.sqrt(x**2+y**2))
    lon = np.arctan2(y, x)
    return lat, lon

def ecef_to_enu(dx, dy, dz, lat, lon):
    """
    將 ECEF 差向量 (dx,dy,dz) 轉換成以當地地面站為原點的 ENU 座標系。
    參數 lat, lon 需以弧度表示。
    回傳值順序為 [East, North, Up]。
    """
    # 轉換矩陣
    t = np.array([
        [-np.sin(lon),              np.cos(lon),               0],
        [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon),  np.cos(lat)],
        [ np.cos(lat)*np.cos(lon),  np.cos(lat)*np.sin(lon),  np.sin(lat)]
    ])
    enu = t @ np.array([dx, dy, dz])
    return enu

def geodetic_to_ecef(lat_deg, lon_deg, alt):
    """
    將經緯度（單位：度）與高度（公尺）轉換為 ECEF 座標 (x, y, z)（單位：公尺）。
    
    使用 WGS84 橢球模型參數：
        a = 6378137.0  (地球赤道半徑, 公尺)
        f = 1 / 298.257223563  (扁率)
    
    參數:
        lat_deg : 緯度（度）
        lon_deg : 經度（度）
        alt     : 高度（公尺）
    
    回傳:
        x, y, z : ECEF 座標（公尺）
    """
    # WGS84 參數
    a = 6378137.0                # 赤道半徑 (公尺)
    f = 1 / 298.257223563        # 扁率
    e2 = 2*f - f**2              # 第一偏心率平方

    # 轉換成弧度
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    
    # 計算卯酉圈曲率半徑
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
    
    # ECEF 座標公式
    x = (N + alt) * np.cos(lat) * np.cos(lon)
    y = (N + alt) * np.cos(lat) * np.sin(lon)
    z = ( (1 - e2) * N + alt ) * np.sin(lat)
    
    return x, y, z

def ground_to_satellite_enu(ground_ecef, sat_ecef):
    """
    計算從地面站到衛星的相對座標（以 ENU 座標系表示）。
    
    參數:
      ground_ecef: 地面站在 ECEF 座標系下的位置，例如由 generate_point() 得到的 [x, y, z]
      sat_ecef: 衛星在 ECEF 座標系下的位置（同上，注意衛星位置是 (R+h)）
      
    回傳:
      ENU 座標 [E, N, U]，表示從地面站出發，東、北、上方向的分量。
    """
    # 將輸入轉換為 numpy array
    ground_ecef = np.array(ground_ecef)
    sat_ecef = np.array(sat_ecef)
    
    # 差向量（從地面站指向衛星）
    diff = sat_ecef - ground_ecef
    
    # 由地面站的 ECEF 座標求出其緯度與經度（假設地球為完美球體）
    lat, lon = ecef_to_geodetic(ground_ecef[0], ground_ecef[1], ground_ecef[2])
    
    # 將 ECEF 差向量轉換成 ENU 座標系
    enu = ecef_to_enu(diff[0], diff[1], diff[2], lat, lon)
    return enu

def compute_elevation_angle(ground_ecef, sat_ecef, degrees=False):
    """
    計算從地面站到衛星的仰角
    
    參數:
      ground_ecef : 地面站的 ECEF 座標 [x, y, z] (公尺)
      sat_ecef    : 衛星的 ECEF 座標 [x, y, z] (公尺)
      degrees     : 是否將結果轉換成角度（True 表示返回度數，False 表示返回弧度）
      
    回傳:
      仰角（elevation angle），當衛星在地平線之上時，角度為正，否則為負
    """
    # 轉換成 numpy 陣列
    ground = np.array(ground_ecef)
    sat = np.array(sat_ecef)
    
    # 計算從地面站指向衛星的向量
    diff = sat - ground
    norm_diff = np.linalg.norm(diff)
    
    # 地面站的 local up 向量（假設地球為球形，可取 ground 向量的方向）
    up = ground / np.linalg.norm(ground)
    
    # 計算 diff 在 up 方向的分量
    dot_up = np.dot(diff, up)
    
    # 仰角：diff 與 local horizontal 平面的夾角，其正弦值為 diff 在 up 方向的分量除以 diff 的總長度
    elev_rad = np.arcsin(dot_up / norm_diff)
    
    if degrees:
        return np.degrees(elev_rad)
    else:
        return elev_rad

def ground_to_atmosphere_distance(ground_ecef, sat_ecef, R=6371000, atmosphere_thickness=10000):
    """
    計算從地面站到大氣外層交點的距離
    
    假設大氣外層為一個球面，其半徑為 (R + atmosphere_thickness) (單位：公尺)。
    參數:
        ground_ecef         : 地面站的 ECEF 座標 [x, y, z] (公尺)
        sat_ecef            : 衛星的 ECEF 座標 [x, y, z] (公尺)
        R                   : 地球半徑 (預設值：6371000 公尺)
        atmosphere_thickness: 大氣厚度 (預設值：10000 公尺)
        
    回傳:
        distance: 從地面站沿著指向衛星的方向，到達大氣外層的直線距離 (公尺)
        
    計算原理:
        設 p0 為地面站位置, v = sat_ecef - ground_ecef 為連線方向向量。
        求解參數 t，使得 ||p0 + t*v|| = R + atmosphere_thickness，
        最終距離為 t * ||v||。
    """
    # 將輸入轉換成 numpy 陣列
    p0 = np.array(ground_ecef, dtype=float)
    v  = np.array(sat_ecef, dtype=float) - p0
    
    # Atmosphere 外層的半徑
    R_atm = R + atmosphere_thickness
    
    # 求解 t，使得 ||p0 + t*v||^2 = R_atm^2
    # 令 a = ||v||^2, b = 2*(p0·v), c = ||p0||^2 - R_atm^2
    a = np.dot(v, v)
    b = 2 * np.dot(p0, v)
    c = np.dot(p0, p0) - R_atm**2
    
    # 求解二次方程式: a*t^2 + b*t + c = 0
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        raise ValueError("無交點：檢查輸入的地面站或大氣厚度參數是否正確。")
    
    # 由於地面站位於球內，t 有一正一負的根，我們取正根
    t1 = (-b + np.sqrt(discriminant)) / (2 * a)
    t2 = (-b - np.sqrt(discriminant)) / (2 * a)
    t = t1 if t1 > 0 else t2  # 選取正的 t
    
    # 從地面站到交點的距離
    distance = t * np.linalg.norm(v)
    return distance

num_rings = 6       # 地面站環數
num_sats_per_ring = 10  # 每個環的衛星數
R = 6.371e6         # 地球半徑 (m)
h = 500e3           # 衛星高度 (m)

# 生成衛星與地面站網絡
S, G = generate_network(num_rings, num_sats_per_ring, R, h)

# 繪製網絡
plot_initial_configuration(G, S, R)



