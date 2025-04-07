#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from math import atan2, sqrt, sin, cos, degrees

def ecef_to_geodetic(x, y, z):
    """
    將 ECEF 座標 (x, y, z) 轉換為 WGS84 地理座標 (lat, lon, alt)
    輸入單位：公尺
    輸出：lat, lon（單位：度），alt（單位：公尺）
    """
    # WGS84 參數
    a = 6378137.0          # 長半軸 (m)
    f = 1 / 298.257223563  # 扁率
    b = a * (1 - f)        # 短半軸
    e_sq = 1 - (b**2)/(a**2)
    
    # 初步估算
    lon = atan2(y, x)
    p = sqrt(x*x + y*y)
    lat = atan2(z, p * (1 - e_sq))
    
    # 迭代改善 (一般 5 次迭代足夠)
    for _ in range(5):
        N = a / sqrt(1 - e_sq * sin(lat)**2)
        alt = p / cos(lat) - N
        lat = atan2(z, p * (1 - e_sq * N / (N + alt)))
    
    return degrees(lat), degrees(lon), alt

def farthest_point_sampling(points, N):
    """
    points: np.array, shape (total_points, 3)，表示所有衛星的座標
    N: int，要選取的衛星數量
    """
    # 隨機選取第一個點的索引
    selected_indices = [np.random.randint(len(points))]
    # 建立一個陣列保存每個點到已選點集合的最小距離，初始設定為無窮大
    distances = np.full(len(points), np.inf)
    
    for _ in range(1, N):
        last_point = points[selected_indices[-1]]
        dist_to_last = np.linalg.norm(points - last_point, axis=1)
        distances = np.minimum(distances, dist_to_last)
        next_index = np.argmax(distances)
        selected_indices.append(next_index)
    
    return selected_indices

def main():
    # 使用 argparse 取得命令列參數
    parser = argparse.ArgumentParser(
        description="從 NORAD TLE 產生的真實 ECEF 衛星資料中，選出位於北美上空且均勻分布的 N 顆衛星"
    )
    parser.add_argument("N", type=int, help="要選取的衛星數量")
    args = parser.parse_args()
    N = args.N

    # 讀取檔案 (檔案內資料已是 ECEF 座標，單位 m)
    with open("dataset/code/output/satellite_coordinates_all.txt", "r") as f:
        satellites = []
        for line in f:
            line = line.strip()
            if line:
                tokens = line.split()
                try:
                    x, y, z = float(tokens[0]), float(tokens[1]), float(tokens[2])
                    satellites.append([x, y, z])
                except (IndexError, ValueError) as e:
                    print(f"解析行失敗: '{line}', error: {e}")
                    continue
    print(f"Loaded {len(satellites)} satellites from file.")

    # 過濾出位於北美上空的衛星
    # 北美洲邊界定義 (簡單版)：
    # 緯度介於 15°N 到 75°N, 經度介於 -170°W 到 -50°W
    north_america_satellites = []
    for rec in satellites:
        x, y, z = rec
        lat, lon, alt = ecef_to_geodetic(x, y, z)
        if 15 <= lat <= 75 and -170 <= lon <= -50:
            north_america_satellites.append(rec)
    print(f"Filtered satellites over North America: {len(north_america_satellites)}")

    if len(north_america_satellites) == 0:
        print("沒有找到符合北美上空條件的衛星。")
        return

    # 檢查 N 是否超過北美區域內的衛星數量
    if N > len(north_america_satellites):
        print("N 超過北美上空的衛星數量，將以全部衛星數量作為輸出。")
        N = len(north_america_satellites)

    # 將過濾後的衛星轉換為 numpy 陣列 (形狀: (total_points, 3))
    coords = np.array(north_america_satellites)

    # 使用 Farthest Point Sampling 選取 N 個均勻分布的衛星索引
    selected_indices = farthest_point_sampling(coords, N)
    selected_records = [north_america_satellites[i] for i in selected_indices]

    # 輸出選取的 N 顆衛星的 ECEF 座標 (m) 到 txt 檔案
    output_filename = "dataset/code/output/satellite_coordinates_selected.txt"
    with open(output_filename, "w") as f:
        f.write(f"{N}\n")
        for rec in selected_records:
            f.write(f"{rec[0]} {rec[1]} {rec[2]}\n")
    print(f"Output {N} satellites over North America to {output_filename}")

if __name__ == "__main__":
    main()
