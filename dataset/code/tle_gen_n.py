#!/usr/bin/env python3
import argparse
from skyfield.api import load, wgs84
import numpy as np

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
        # 取出最後選的點
        last_point = points[selected_indices[-1]]
        # 更新每個點到已選集合中最接近點的距離
        dist_to_last = np.linalg.norm(points - last_point, axis=1)
        distances = np.minimum(distances, dist_to_last)
        # 選取距離最大的那個點
        next_index = np.argmax(distances)
        selected_indices.append(next_index)
    
    return selected_indices

def main():
    # 使用 argparse 取得命令列參數
    parser = argparse.ArgumentParser(
        description="從 NORAD TLE 資料中選取均勻分布的 N 顆衛星"
    )
    parser.add_argument("N", type=int, help="要選取的衛星數量")
    args = parser.parse_args()
    N = args.N

    if N > 11200:
        print("N 請設定更小的數值")
        return

    # read from file (already ecef with m)
    with open("dataset/code/output/satellite_coordinates_all.txt", "r") as f:
        satellites = []
        for line in f:
            line = line.strip()
            if line:
                tokens = line.split()
                try:
                    parsed = [
                        
                        float(tokens[0]),
                        float(tokens[1]),
                        float(tokens[2])
                    ]
                    satellites.append(parsed)
                except (IndexError, ValueError) as e:
                    print(f"解析行失敗: '{line}', error: {e}")
                    continue
        
    print(f"Loaded {len(satellites)} satellites")

    # 將所有衛星座標轉換為 numpy 陣列 (形狀: (total_points, 3))
    coords = np.array([[rec[0], rec[1], rec[2]] for rec in satellites])

    # 檢查 N 是否超過資料點數
    if N > len(coords):
        print("N 超過衛星總數，請設定更小的數值")
        N = len(coords)

    # 使用 Farthest Point Sampling 選取 N 個衛星的索引
    selected_indices = farthest_point_sampling(coords, N)
    selected_records = [satellites[i] for i in selected_indices]

    

    # 輸出選取的 N 顆均勻分布的衛星座標到另一個 txt 檔案
    with open("dataset/code/output/satellite_coordinates_selected.txt", "w") as f:
        f.write(f"{N}\n")
        for rec in selected_records:
            f.write(f"{rec[0]} {rec[1]} {rec[2]}\n")
    print(f"output {N} sat to txt done")

if __name__ == "__main__":
    main()
