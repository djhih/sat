#!/usr/bin/env python
import math

def geodetic_to_ecef(lat_deg, lon_deg, alt=0, R=6371e3):
    """
    將地理座標 (lat, lon, alt) 轉換為 ECEF 坐標 (x, y, z)。
    假設地球為完美球體，R 為地球半徑 (公尺)，lat, lon 需以度為單位。
    """
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    x = (R + alt) * math.cos(lat) * math.cos(lon)
    y = (R + alt) * math.cos(lat) * math.sin(lon)
    z = (R + alt) * math.sin(lat)
    return x, y, z

def main():
    states = {
        "Alabama": ("Montgomery", 32.361538, -86.279118),
        "Alaska": ("Juneau", 58.301935, -134.419740),
        "Arizona": ("Phoenix", 33.448457, -112.073844),
        "Arkansas": ("Little Rock", 34.736009, -92.331122),
        "California": ("Sacramento", 38.555605, -121.468926),
        "Colorado": ("Denver", 39.7391667, -104.984167),
        "Connecticut": ("Hartford", 41.763711, -72.685093),
        "Delaware": ("Dover", 39.158168, -75.524368),
        "Florida": ("Tallahassee", 30.438118, -84.281296),
        "Georgia": ("Atlanta", 33.749027, -84.388229),
        "Hawaii": ("Honolulu", 21.307442, -157.857376),
        "Idaho": ("Boise", 43.613739, -116.237651),
        "Illinois": ("Springfield", 39.783250, -89.650373),
        "Indiana": ("Indianapolis", 39.790942, -86.147685),
        "Iowa": ("Des Moines", 41.590939, -93.620866),
        "Kansas": ("Topeka", 39.04, -95.69),
        "Kentucky": ("Frankfort", 38.197274, -84.86311),
        "Louisiana": ("Baton Rouge", 30.45809, -91.140229),
        "Maine": ("Augusta", 44.323535, -69.765261),
        "Maryland": ("Annapolis", 38.972945, -76.501157),
        "Massachusetts": ("Boston", 42.2352, -71.0275),
        "Michigan": ("Lansing", 42.7335, -84.5467),
        "Minnesota": ("St. Paul", 44.95, -93.094),
        "Mississippi": ("Jackson", 32.32, -90.207),
        "Missouri": ("Jefferson City", 38.572954, -92.189283),
        "Montana": ("Helena", 46.595805, -112.027031),
        "Nebraska": ("Lincoln", 40.809868, -96.675345),
        "Nevada": ("Carson City", 39.160949, -119.753877),
        "New Hampshire": ("Concord", 43.220093, -71.549127),
        "New Jersey": ("Trenton", 40.221741, -74.756138),
        "New Mexico": ("Santa Fe", 35.667231, -105.964575),
        "New York": ("Albany", 42.659829, -73.781339),
        "North Carolina": ("Raleigh", 35.771, -78.638),
        "North Dakota": ("Bismarck", 46.813343, -100.779004),
        "Ohio": ("Columbus", 39.962245, -83.000647),
        "Oklahoma": ("Oklahoma City", 35.482309, -97.534994),
        "Oregon": ("Salem", 44.931109, -123.029159),
        "Pennsylvania": ("Harrisburg", 40.269789, -76.875613),
        "Rhode Island": ("Providence", 41.82355, -71.422132),
        "South Carolina": ("Columbia", 34.000, -81.0348),
        "South Dakota": ("Pierre", 44.367966, -100.336378),
        "Tennessee": ("Nashville", 36.165, -86.784),
        "Texas": ("Austin", 30.266667, -97.75),
        "Utah": ("Salt Lake City", 40.7547, -111.892622),
        "Vermont": ("Montpelier", 44.26639, -72.57194),
        "Virginia": ("Richmond", 37.54, -77.46),
        "Washington": ("Olympia", 47.042418, -122.893077),
        "West Virginia": ("Charleston", 38.349497, -81.633294),
        "Wisconsin": ("Madison", 43.074722, -89.384444),
        "Wyoming": ("Cheyenne", 41.145548, -104.802042)
    }

    # print("State, Capital, Lat (deg), Lon (deg), ECEF X (m), ECEF Y (m), ECEF Z (m)")
    with open("dataset/code/output/gs_loc.txt", "w") as f:
        f.write(f"{len(states)}\n")
        for state, (capital, lat, lon) in states.items():
            x, y, z = geodetic_to_ecef(lat, lon, alt=0)
            f.write(f"{x:.2f} {y:.2f} {z:.2f}\n")

if __name__ == "__main__":
    main()
