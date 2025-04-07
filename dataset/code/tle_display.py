from skyfield.api import load
import numpy as np
import plotly.graph_objects as go

# 載入所有活躍衛星 TLE
stations_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle'
satellites = load.tle_file(stations_url)
print(f'Loaded {len(satellites)} satellites')

# 時間設定
ts = load.timescale()
t = ts.now()

# 地球半徑
R_earth = 6371  # km

# 計算衛星 ECI 座標
x_list, y_list, z_list = [], [], []

for sat in satellites:
    geocentric = sat.at(t)
    position = geocentric.position.km
    x_list.append(position[0])
    y_list.append(position[1])
    z_list.append(position[2])

# 產生衛星分布點
satellite_trace = go.Scatter3d(
    x=x_list, y=y_list, z=z_list,
    mode='markers',
    marker=dict(size=2, color='red'),
    name='Satellites'
)

# 畫地球（球面）
u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
x_sphere = R_earth * np.cos(u) * np.sin(v)
y_sphere = R_earth * np.sin(u) * np.sin(v)
z_sphere = R_earth * np.cos(v)

earth_trace = go.Surface(
    x=x_sphere, y=y_sphere, z=z_sphere,
    colorscale='Blues',
    opacity=0.7,
    showscale=False,
    name='Earth'
)

# 畫圖
fig = go.Figure(data=[earth_trace, satellite_trace])
fig.update_layout(
    title='3D Satellite Distribution (ECI Coordinates)',
    scene=dict(
        xaxis_title='X (km)',
        yaxis_title='Y (km)',
        zaxis_title='Z (km)',
        aspectmode='data'
    ),
    margin=dict(l=0, r=0, b=0, t=40)
)
fig.show()
