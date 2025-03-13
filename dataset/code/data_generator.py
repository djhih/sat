#!/usr/bin/env python
import numpy as np
import math
import random
from numpy.linalg import norm
from scipy.linalg import expm
from scipy.optimize import fminbound, root_scalar
import warnings

# -------------------------
# 基本物理與座標轉換函式
# -------------------------
def get_constants(h):
    # 輸入：衛星高度 (m)
    G_const = 6.67408e-11  # 重力常數 (SI)
    R_earth = 6.371e6      # 地球半徑 (m)
    M = 5.972e24           # 地球質量 (kg)
    Omega = 2*np.pi/86400.  # 地球自轉速率 (rad/s)
    t = 10e3               # 大氣厚度 (m)
    omega = np.sqrt(G_const * M) / (R_earth + h)**(3./2.)
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
    # 參數設定
    num_rings = 3            # 環數
    num_sats_per_ring = 5    # 每環衛星數
    R_earth = 6.371e6        # 地球半徑 (m)
    h = 500e3                # 衛星高度 (m)
    
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
    # -------------------------
    sat_parameters = []  # 每個元素為 (x, y, z, tmp_gs, tmp_dis, tmp_ang, tmp_h_atm, tmp_fid, tmp_gen_rate)
    for i, sat_pos in enumerate(satellites):
        tmp_gs = random.randint(0, G_total - 1)
        gs_pos = ground_stations[tmp_gs]
        tmp_dis = np.linalg.norm(np.array(sat_pos) - np.array(gs_pos))
        tmp_ang = compute_elevation_angle(gs_pos, sat_pos, degrees=True)
        tmp_h_atm = ground_to_atmosphere_distance(gs_pos, sat_pos, R_earth, atmosphere_thickness=10000)
        tmp_fid = random.uniform(0.5, 1.0)
        tmp_gen_rate = random.uniform(5, 15)
        sat_parameters.append((sat_pos[0], sat_pos[1], sat_pos[2],
                               tmp_gs, tmp_dis, tmp_ang, tmp_h_atm, tmp_fid, tmp_gen_rate))
    
    # -------------------------
    # 輸出 dataset 至檔案 (例如 dataset.txt)
    # 格式：
    # 第一行：S G R
    # 接下來 S 行：每行 9 個數值 (衛星位置與參數)
    # 接下來 G 行：每行 3 個數值 (地面站位置)
    # 接下來 R 行：每行 2 個數值 (需求對：gs1 gs2)
    # -------------------------
    output_filename = "dataset.txt"
    with open(output_filename, "w") as f:
        # 輸出第一行：S G R
        f.write(f"{S_total} {G_total} {R_total}\n")
        # 輸出每顆衛星資料
        for params in sat_parameters:
            f.write(" ".join(map(str, params)) + "\n")
        # 輸出每個地面站資料
        for gs in ground_stations:
            f.write(" ".join(map(str, gs)) + "\n")
        # 輸出每筆需求 (gs1 gs2)
        for (gs1, gs2) in requirements:
            f.write(f"{gs1} {gs2}\n")
    
    print(f"Dataset generated in {output_filename}")
    
if __name__ == '__main__':
    main()
