#!/usr/bin/env python3
from skyfield.api import load, wgs84
import numpy as np

def main():
    stations_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle'
    satellites = load.tle_file(stations_url)
    print(f"Loaded {len(satellites)} satellites")

    ts = load.timescale()
    t = ts.now()

    records = []
    for sat in satellites:
        geocentric = sat.at(t)
        ecef = wgs84.subpoint(geocentric).itrs_xyz.km  
        # km → m
        x_m, y_m, z_m = [coord * 1000 for coord in ecef]
        if np.isnan(x_m) or np.isnan(y_m) or np.isnan(z_m):
            continue
        records.append({
            "satellite_name": sat.name,
            "x_m": x_m,
            "y_m": y_m,
            "z_m": z_m
        })

    with open("dataset/code/output/satellite_coordinates_all.txt", "w") as f:
        for rec in records:
            f.write(f"{rec['x_m']} {rec['y_m']} {rec['z_m']}\n")
    
    print("write all sat to txt done")

if __name__ == "__main__":
    main()
