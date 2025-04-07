#!/usr/bin/env python
import numpy as np
import math
import random
from numpy.linalg import norm
from scipy.linalg import expm
from scipy.optimize import fminbound, root_scalar
import warnings

# -------------------------
# eta 計算
# eta_sg 計算的是衛星到地面站的穿透率
# eta_atm 計算的是大氣穿透
# eta_tot 是總穿透率也就是 eta_sg * eta_atm
# -------------------------
def eta_sg(L,lam=810e-9,w0=0.025,rg=0.75): # L: dis in meter
    '''
    L : 衛星到地面的距離 m。
    Satellite-to-ground transmittance. All distance quantities are
    in meters.
    '''

    LR=np.pi*w0**2/lam

    wg=w0*np.sqrt(1+L**2/LR**2)

    return 1-np.exp(-2*rg**2/wg**2)
def eta_Atm(L,h,etaZen=0.5): # L: dis, h: height in meter

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
def eta_Tot(L,h): # L: dis, h: height

    # Last modified: 18 July 2019

    '''
    Total transmittance, which is a product of the atmospheric transmittance
    and the beam spreading transmittance. L is the link distance (distance between
    the ground station and the satellite) in meters, and h is the altitude of the
    satellite (in meters).
    '''

    return eta_sg(L)*eta_Atm(L,h)

# -------------------------
# 產生地面站與衛星之間建立量子糾纏的保真度
# 產生地面戰與衛星之間的 generation rate
# -------------------------
def fidelity_Fij(theta_k1, theta_e, n_sg, n_ij, F0):
    """
    根據以下條件計算 F_{i,j}:
    
        F_{i,j} = 1/4 * (1 + (4*F0 - 1) / (1 + n_sg / n_ij)^2),
        若 theta_k1 >= theta_e，且 (k1, k2) ∈ j；
        否則 = 0。
    
    參數說明:
      theta_k1: 衛星或地面站的兩個仰角 (或其他角度)，用於判斷是否 >= theta_e
      theta_e: 判斷的門檻角度
      n_sg: 公式中的 n_{sg}，通常為某個固定常數或可變參數
      n_ij: 公式中的 n_{i,j}，視情況可來自鏈路特性或隨機生成
      F0: 公式中的 F0，通常是某個固定常數，用於計算保真度修正
    
    回傳:
      Fij (float): 按照上式計算得到的數值
    """
    # 條件判斷: 若兩個角度都 >= theta_e，才計算公式
    if theta_k1 >= theta_e:
        # 避免除以零，假設 n_ij>0
        return 0.25 * (
            1.0 + (4.0 * F0 - 1.0) / ((1.0 + n_sg / n_ij)**2)
        )
    else:
        return 0.0
def gen_rate(d):
    ratio = ((float)(d) / (float)(50000)) 
    return 20 * (1/(ratio**2)) 
    ## generate 我先基質先設 20 

# -------------------------
# 計算給定距離 d 下的最好透射率的高度
# generate_link_constraint(d)
# 先計算最佳衛星高度 h_opt（使用 opt_alt_arc(d)）。
# 再計算並返回最佳 h_opt 下的鏈路距離使用 link_distance(d, h_opt)
# 要注意這邊的高度是算出來的不是用我們給定的資料
# -------------------------
def link_distance(d,h):
    '''
    計算其中一個地面站到衛星的距離
    假設當前衛星於兩個地面站中間，且衛星於兩地面站中點之高空
    傳入 D 代表兩地面站距離, h 代表衛星高度

    Returns the link distance L between two ground stations separated by an 
    arc-distance of d (in meters) and a satellite at an altitude of
    h (in meters) located at the midpoint of the ground stations (so both
    ground stations are the same distance away from the satellite).
    '''

    _,R,_,_,_,_,_=get_constants(h)

    return np.sqrt(4*R*(R+h)*np.sin(d/(4*R))**2+h**2)
def opt_alt_arc(d):
    '''
     計算最佳衛星高度，使得兩個地面站之間的光傳輸效率 eta_Tot 最大。
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
    '''
    Given a maximum allowed arc-distance separation d (in meters) between
    two ground stations, this function generates the maximum allowed distance
    between a ground station and a satellite based on the optimal altiude
    of the satellite.
    '''

    h_opt=opt_alt_arc(d)

    return link_distance(d,h_opt)

# -------------------------
# 基本物理與座標轉換函式
# -------------------------
def get_constants(h):
    # 輸入：衛星高度 (m)
    G_const = 6.67408e-11   # 重力常數 (SI)
    R_earth = 6.371e6       # 地球半徑 (m)
    M = 5.972e24            # 地球質量 (kg)
    Omega = 2*np.pi/86400.  # 地球自轉速率 (rad/s)
    t = 10e3                # 大氣厚度 (m)
    omega = np.sqrt(G_const * M) / (R_earth + h)**(3./2.)  # Rotation speed of the satellites in radians/s
    #h=1000e3 # Altitude of the satellites in meters
    #period=2*np.pi/omega # Period of each satellite in seconds.
    return G_const, R_earth, M, Omega, t, h, omega

def generate_point(R, theta, phi): 
    # 將球座標 (R, theta, phi) 轉換成 Cartesian 座標
    return [R * np.sin(phi) * np.cos(theta),
            R * np.sin(phi) * np.sin(theta),
            R * np.cos(phi)]

# -------------------------
# 產生衛星與地面站網路 (環狀結構)
# -------------------------
def generate_network(num_rings, num_sats, R, h, offset_angle=0, multiple_alt=False):
    S = {}
    G = {}
    count1 = 0
    for m in range(num_rings):
        count2 = 0
        count1 += 1
        S[count1] = {}
        G[count1] = {}
        for k in range(1, num_sats+1):
            count2 += 1
            if multiple_alt:
                S[count1][count2] = {}
                for alt in h:  # h 為一個高度列表
                    S[count1][count2][alt/1000] = generate_point(R+alt, m*np.pi/num_rings, k*2*np.pi/num_sats+offset_angle)
            else:
                S[count1][count2] = generate_point(R+h, m*np.pi/num_rings, k*2*np.pi/num_sats+offset_angle)
            G[count1][count2] = generate_point(R, m*np.pi/num_rings+np.pi/(2*num_rings), k*2*np.pi/num_sats+offset_angle)
        # 儲存該環的旋轉軸（後續可能用於模擬動態）
        S[count1]['axis'] = [-np.sin(m*np.pi/num_rings), np.cos(m*np.pi/num_rings), 0]
    return S, G

# -------------------------
# 幾何與傳輸計算函式
# -------------------------
def compute_elevation_angle(ground_ecef, sat_ecef, degrees=False):
    ground = np.array(ground_ecef)
    sat = np.array(sat_ecef)
    diff = sat - ground
    norm_diff = np.linalg.norm(diff)
    up = ground / np.linalg.norm(ground)
    dot_up = np.dot(diff, up)
    elev_rad = np.arcsin(dot_up / norm_diff)
    if degrees:
        return np.degrees(elev_rad)
    else:
        return elev_rad

def ground_to_atmosphere_distance(ground_ecef, sat_ecef, R=6371000, atmosphere_thickness=10000):
    p0 = np.array(ground_ecef, dtype=float)
    v = np.array(sat_ecef, dtype=float) - p0
    R_atm = R + atmosphere_thickness
    a = np.dot(v, v)
    b = 2 * np.dot(p0, v)
    c = np.dot(p0, p0) - R_atm**2
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        raise ValueError("無交點，請檢查輸入參數")
    t1 = (-b + np.sqrt(discriminant)) / (2 * a)
    t2 = (-b - np.sqrt(discriminant)) / (2 * a)
    t = t1 if t1 > 0 else t2
    distance = t * np.linalg.norm(v)
    return distance

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
# -------------------------
# 產生地面站對 (Requirement) 配對
# -------------------------
def generate_ground_station_pairs(G):
    # 利用環狀網路的結構產生配對
    G_pair_labels = []
    num_rings = len(G)
    for m in list(G.keys()):
        num_ground_ring = len(G[m])
        for k in list(G[m].keys()):
            if m == num_rings:
                gs1 = (m, k)
                gs2 = (1, (k % num_ground_ring) + 1)
            else:
                gs1 = (m, k)
                gs2 = (m+1, k)
            G_pair_labels.append((gs1, gs2))
    return G_pair_labels

# -------------------------
# 主程式：產生 dataset 並輸出成文字檔
# -------------------------
def main():
    num_rings = 3            # 環數
    num_sats_per_ring = 3    # 每環衛星數
    R_earth = 6.371e6         # 地球半徑 (m)
    h = 500e3                 # 衛星高度 (m)
    
    # 產生環狀網路 (S 與 G 為字典，key 為環編號)
    S_dict, G_dict = generate_network(num_rings, num_sats_per_ring, R_earth, h)
    
    # 將地面站展平成一維列表，同時建立 (ring, index) 到唯一 id 的對應
    ground_stations = []
    gs_mapping = {}  # (ring, index) -> id
    gs_id = 0
    for m in sorted(G_dict.keys()):
        for k in sorted(G_dict[m].keys()):
            ground_stations.append(G_dict[m][k])
            gs_mapping[(m, k)] = gs_id
            gs_id += 1
    
    # 將衛星展平成一維列表，同時建立 (ring, index) 到唯一 id 的對應
    satellites = []
    sat_mapping = {}  # (ring, index) -> id
    sat_id = 0
    for m in sorted(S_dict.keys()):
        for k in sorted([key for key in S_dict[m].keys() if isinstance(key, int)]):
            if k == 'axis':
                continue
            satellites.append(S_dict[m][k])
            sat_mapping[(m, k)] = sat_id
            sat_id += 1

    S_total = len(satellites)
    G_total = len(ground_stations)
    
    # 產生需求對 (Requirement)
    req_pair_labels = generate_ground_station_pairs(G_dict)
    requirements = []
    for pair in req_pair_labels:
        gs1_tuple, gs2_tuple = pair
        if gs1_tuple in gs_mapping and gs2_tuple in gs_mapping:
            requirements.append((gs_mapping[gs1_tuple], gs_mapping[gs2_tuple]))
    R_total = len(requirements)
    
    # -------------------------
    # 為每顆衛星產生參數：
    # 選擇一個隨機地面站作為 tmp_gs，並計算：
    #   tmp_dis：歐式距離
    #   tmp_ang：仰角 (degrees)
    #   tmp_h_atm：地面站到大氣外層距離 (m)
    #   tmp_fid：保真度 (隨機取 [0.5, 1.0])
    #   tmp_gen_rate：生成速率 (隨機取 [5, 15])
    #   base_dis : 作為計算 gen_rate 的最短距離，這樣就可以讓 gen_rate 可以用 d^2 反比去處理
    # -------------------------
    sat_parameters = []  # 衛星資料
    gs_sat_data = {}  # (gs, sat) -> 參數
    min_dis = 0.0

    for i, sat_pos in enumerate(satellites):
        sat_parameters.append((sat_pos[0], sat_pos[1], sat_pos[2]))

        for j, gs_pos in enumerate(ground_stations):
            tmp_dis = np.linalg.norm(np.array(sat_pos) - np.array(gs_pos))  # 用 x, y, z 公式直接算的
            tmp_ang = compute_elevation_angle(gs_pos, sat_pos, degrees=True)  # 仰角
            tmp_h_atm = ground_to_atmosphere_distance(gs_pos, sat_pos, R_earth, atmosphere_thickness=10000)  # 大氣層距離
            fid0 = 0.8 # 忘記怎麼算了
            height = np.linalg.norm(np.array(sat_pos)) - R_earth
            tmp_fid = fidelity_Fij(tmp_ang, 20, eta_Tot(tmp_dis, height), 10, fid0) 
            tmp_gen_rate = gen_rate(tmp_dis)  
            min_dis = min(min_dis, tmp_dis)
            gs_sat_data[(j, i)] = (int(tmp_dis), int(tmp_ang), int(tmp_h_atm), tmp_fid, tmp_gen_rate)

    output_filename = "dataset.txt"
    with open(output_filename, "w") as f:
        f.write("S G R\n")
        f.write(f"{S_total} {G_total} {R_total}\n")

        # sat
        f.write("\nSatellite location\n")
        for sat in sat_parameters:
            f.write(f"{sat[0]} {sat[1]} {sat[2]}\n")

        # gs
        f.write("\nGS location\n")
        for gs in ground_stations:
            f.write(f"{gs[0]} {gs[1]} {gs[2]}\n")

        # gs pair
        f.write("\n Request\n")
        for (gs1, gs2) in requirements:
            f.write(f"{gs1} {gs2}\n")

        # gs-sat data
        f.write("\nData dis, ang, h_atm, fid, gen_rate\n")
        for i in range(S_total):
            for j in range(G_total):
                tmp_dis, tmp_ang, tmp_h_atm, tmp_fid, tmp_gen_rate = gs_sat_data[(j, i)]
                f.write(f"{tmp_dis} {tmp_ang} {tmp_h_atm} {tmp_fid} {tmp_gen_rate}\n")

    print(f"測資已輸出至 {output_filename}")
    
if __name__ == '__main__':
    main()
