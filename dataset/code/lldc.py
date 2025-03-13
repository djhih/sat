import math

def haversine(lat1, lon1, lat2, lon2):
    """
    計算兩個經緯度之間的地表距離（單位：公里）
    :param lat1: 第一個點的緯度（度）
    :param lon1: 第一個點的經度（度）
    :param lat2: 第二個點的緯度（度）
    :param lon2: 第二個點的經度（度）
    :return: 兩點之間的距離（公里）
    """
    # 將角度轉換為弧度
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    
    # Haversine 公式計算
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    # 地球半徑（以公里為單位）
    R = 6371.0
    distance = R * c
    
    return distance

# 測試範例
lat1, lon1 = 25.0330, 121.5654  # 台北 101
lat2, lon2 = 22.3964, 114.1095  # 香港

distance = haversine(lat1, lon1, lat2, lon2)
print(f"兩點之間的距離: {distance:.2f} 公里")
