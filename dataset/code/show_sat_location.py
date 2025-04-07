import numpy as np
import matplotlib.pyplot as plt

def read_satellite_positions(filename):
    sat_positions = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # 跳過空行
            parts = line.split()
            if len(parts) != 3:
                print(f"警告: 跳過格式錯誤的行: {line}")
                continue
            try:
                x, y, z = map(float, parts)
                sat_positions.append([x, y, z])
            except ValueError:
                print(f"警告: 無法解析行: {line}")
                continue
    return np.array(sat_positions)

def plot_earth_and_satellites(sat_positions):
    """
    畫出地球 (以半徑 6.371e6 公尺的球面表示) 以及衛星位置
    """
    # 地球半徑
    R_earth = 6.371e6

    # 建立地球的球面
    phi = np.linspace(0, np.pi, 50)       # 夾角：0 到 pi
    theta = np.linspace(0, 2*np.pi, 50)     # 方位角：0 到 2pi
    phi, theta = np.meshgrid(phi, theta)
    x_earth = R_earth * np.sin(phi) * np.cos(theta)
    y_earth = R_earth * np.sin(phi) * np.sin(theta)
    z_earth = R_earth * np.cos(phi)

    # 建立 3D 圖形
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 畫出地球（半透明藍色）
    ax.plot_surface(x_earth, y_earth, z_earth, color='b', alpha=0.3, rstride=4, cstride=4, linewidth=0)

    # 畫出衛星位置 (紅色點)
    ax.scatter(sat_positions[:, 0], sat_positions[:, 1], sat_positions[:, 2],
               color='r', s=50, label='Satellites')

    # 設定軸標籤與標題
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('Earth and Satellite Positions')
    ax.legend()

    # 設定座標軸等比例顯示
    ax.set_box_aspect([1,1,1])
    plt.show()

def main():
    filename = 'gs_loc.txt'
    sat_positions = read_satellite_positions(filename)
    if sat_positions.size == 0:
        print("沒有讀取到有效的衛星資料！")
        return
    plot_earth_and_satellites(sat_positions)

if __name__ == '__main__':
    main()
