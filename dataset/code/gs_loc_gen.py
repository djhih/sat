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
    govs = {
        "Afghanistan": ("Kabul", 34.516667, 69.183333),
        "Albania": ("Tirana", 41.3275, 19.8189),
        "Algeria": ("Algiers", 36.75, 3.04),
        "Andorra": ("Andorra la Vella", 42.50779, 1.52109),
        "Angola": ("Luanda", -8.838333, 13.234444),
        "Antigua and Barbuda": ("St. John's", 17.12741, -61.84677),
        "Argentina": ("Buenos Aires", -34.603722, -58.381592),
        "Armenia": ("Yerevan", 40.18111, 44.51361),
        "Australia": ("Canberra", -35.282, 149.128684),
        "Austria": ("Vienna", 48.208174, 16.373819),
        "Azerbaijan": ("Baku", 40.409264, 49.867092),
        "Bahamas": ("Nassau", 25.03428, -77.39628),
        "Bahrain": ("Manama", 26.228516, 50.586049),
        "Bangladesh": ("Dhaka", 23.728, 90.394),
        "Barbados": ("Bridgetown", 13.0975, -59.616667),
        "Belarus": ("Minsk", 53.9, 27.566667),
        "Belgium": ("Brussels", 50.850346, 4.351721),
        "Belize": ("Belmopan", 17.251, -88.766667),
        "Benin": ("Porto-Novo", 6.483333, 2.616667),
        "Bhutan": ("Thimphu", 27.4728, 89.639286),
        "Bolivia": ("La Paz", -16.5, -68.15),
        "Bosnia and Herzegovina": ("Sarajevo", 43.848642, 18.356443),
        "Botswana": ("Gaborone", -24.658056, 25.908611),
        "Brazil": ("Brasília", -15.793889, -47.882778),
        "Brunei": ("Bandar Seri Begawan", 4.9031, 114.9398),
        "Bulgaria": ("Sofia", 42.697708, 23.321868),
        "Burkina Faso": ("Ouagadougou", 12.3714, -1.51966),
        "Burundi": ("Gitega", -3.4264, 29.9306),
        "Cambodia": ("Phnom Penh", 11.55, 104.916667),
        "Cameroon": ("Yaoundé", 3.866667, 11.516667),
        "Canada": ("Ottawa", 45.4215, -75.6972),
        "Cape Verde": ("Praia", 14.933333, -23.516667),
        "Central African Republic": ("Bangui", 4.366667, 18.583333),
        "Chad": ("N'Djamena", 12.106389, 15.044444),
        "Chile": ("Santiago", -33.45, -70.666667),
        "China": ("Beijing", 39.9042, 116.4074),
        "Colombia": ("Bogotá", 4.711, -74.0721),
        "Comoros": ("Moroni", -11.717, 43.247),
        "Congo (Republic of the)": ("Brazzaville", -4.266667, 15.283333),
        "Congo (Democratic Republic of the)": ("Kinshasa", -4.325, 15.322222),
        "Costa Rica": ("San José", 9.9281, -84.0907),
        "Croatia": ("Zagreb", 45.815, 15.9785),
        "Cuba": ("Havana", 23.1136, -82.3666),
        "Cyprus": ("Nicosia", 35.1856, 33.3823),
        "Czech Republic": ("Prague", 50.087465, 14.421254),
        "Denmark": ("Copenhagen", 55.6761, 12.5683),
        "Djibouti": ("Djibouti", 11.588, 43.145),
        "Dominica": ("Roseau", 15.3, -61.4),
        "Dominican Republic": ("Santo Domingo", 18.4861, -69.9312),
        "East Timor": ("Dili", -8.556944, 125.560278),
        "Ecuador": ("Quito", -0.22985, -78.52495),
        "Egypt": ("Cairo", 30.0444, 31.2357),
        "El Salvador": ("San Salvador", 13.69294, -89.218191),
        "Equatorial Guinea": ("Malabo", 3.75, 8.783333),
        "Eritrea": ("Asmara", 15.322876, 38.925052),
        "Estonia": ("Tallinn", 59.437, 24.7536),
        "Eswatini": ("Mbabane", -26.316667, 31.133333),
        "Ethiopia": ("Addis Ababa", 8.9806, 38.7578),
        "Fiji": ("Suva", -18.1416, 178.4419),
        "Finland": ("Helsinki", 60.1699, 24.9384),
        "France": ("Paris", 48.8566, 2.3522),
        "Gabon": ("Libreville", 0.4162, 9.4673),
        "Gambia": ("Banjul", 13.4549, -16.579),
        "Georgia": ("Tbilisi", 41.7151, 44.8271),
        "Germany": ("Berlin", 52.52, 13.4050),
        "Ghana": ("Accra", 5.6037, -0.1870),
        "Greece": ("Athens", 37.9838, 23.7275),
        "Grenada": ("St. George's", 12.0561, -61.7488),
        "Guatemala": ("Guatemala City", 14.6349, -90.5069),
        "Guinea": ("Conakry", 9.509167, -13.712222),
        "Guinea-Bissau": ("Bissau", 11.8633, -15.5977),
        "Guyana": ("Georgetown", 6.8013, -58.1551),
        "Haiti": ("Port-au-Prince", 18.5392, -72.335),
        "Honduras": ("Tegucigalpa", 14.0723, -87.1921),
        "Hungary": ("Budapest", 47.4979, 19.0402),
        "Iceland": ("Reykjavik", 64.1466, -21.9426),
        "India": ("New Delhi", 28.6139, 77.2090),
        "Indonesia": ("Jakarta", -6.2088, 106.8456),
        "Iran": ("Tehran", 35.6892, 51.3890),
        "Iraq": ("Baghdad", 33.3128, 44.3615),
        "Ireland": ("Dublin", 53.3498, -6.2603),
        "Israel": ("Jerusalem", 31.7683, 35.2137),
        "Italy": ("Rome", 41.9029, 12.4964),
        "Ivory Coast": ("Yamoussoukro", 6.816667, -5.266667),
        "Jamaica": ("Kingston", 17.9714, -76.7936),
        "Japan": ("Tokyo", 35.6895, 139.6917),
        "Jordan": ("Amman", 31.9539, 35.9106),
        "Kazakhstan": ("Nur-Sultan", 51.1605, 71.4704),
        "Kenya": ("Nairobi", -1.2921, 36.8219),
        "Kiribati": ("Tarawa", 1.4518, 173.02),
        "North Korea": ("Pyongyang", 39.0194, 125.738),
        "South Korea": ("Seoul", 37.5665, 126.9780),
        "Kosovo": ("Pristina", 42.6629, 21.1655),
        "Kuwait": ("Kuwait City", 29.3759, 47.9774),
        "Kyrgyzstan": ("Bishkek", 42.8746, 74.6122),
        "Laos": ("Vientiane", 17.966667, 102.6),
        "Latvia": ("Riga", 56.946, 24.1059),
        "Lebanon": ("Beirut", 33.8938, 35.5018),
        "Lesotho": ("Maseru", -29.316667, 27.483333),
        "Liberia": ("Monrovia", 6.3, -10.8),
        "Libya": ("Tripoli", 32.8872, 13.1913),
        "Liechtenstein": ("Vaduz", 47.141, 9.5215),
        "Lithuania": ("Vilnius", 54.6894, 25.2799),
        "Luxembourg": ("Luxembourg", 49.6116, 6.1319),
        "Madagascar": ("Antananarivo", -18.8792, 47.5079),
        "Malawi": ("Lilongwe", -13.9626, 33.7741),
        "Malaysia": ("Kuala Lumpur", 3.139, 101.6869),
        "Maldives": ("Malé", 4.1755, 73.5093),
        "Mali": ("Bamako", 12.6392, -8.00289),
        "Malta": ("Valletta", 35.8997, 14.5146),
        "Marshall Islands": ("Majuro", 7.0897, 171.38),
        "Mauritania": ("Nouakchott", 18.0858, -15.9785),
        "Mauritius": ("Port Louis", -20.1609, 57.5012),
        "Mexico": ("Mexico City", 19.4326, -99.1332),
        "Micronesia": ("Palikir", 6.9147, 158.1610),
        "Moldova": ("Chisinau", 47.0105, 28.8638),
        "Monaco": ("Monaco", 43.7384, 7.4246),
        "Mongolia": ("Ulaanbaatar", 47.8864, 106.9057),
        "Montenegro": ("Podgorica", 42.4304, 19.2594),
        "Morocco": ("Rabat", 34.020882, -6.84165),
        "Mozambique": ("Maputo", -25.9692, 32.5732),
        "Myanmar": ("Naypyidaw", 19.7633, 96.0785),
        "Namibia": ("Windhoek", -22.5609, 17.0658),
        "Nauru": ("Yaren", -0.5477, 166.9209),
        "Nepal": ("Kathmandu", 27.7172, 85.324),
        "Netherlands": ("Amsterdam", 52.3676, 4.9041),
        "New Zealand": ("Wellington", -41.2865, 174.7762),
        "Nicaragua": ("Managua", 12.136389, -86.251389),
        "Niger": ("Niamey", 13.5125, 2.112222),
        "Nigeria": ("Abuja", 9.05785, 7.49508),
        "North Macedonia": ("Skopje", 41.99646, 21.43141),
        "Norway": ("Oslo", 59.9139, 10.7522),
        "Oman": ("Muscat", 23.5859, 58.4059),
        "Pakistan": ("Islamabad", 33.6844, 73.0479),
        "Palau": ("Ngerulmud", 7.5, 134.624),
        "Panama": ("Panama City", 8.983333, -79.516667),
        "Papua New Guinea": ("Port Moresby", -9.4438, 147.1803),
        "Paraguay": ("Asunción", -25.2637, -57.5759),
        "Peru": ("Lima", -12.0464, -77.0428),
        "Philippines": ("Manila", 14.5995, 120.9842),
        "Poland": ("Warsaw", 52.2297, 21.0122),
        "Portugal": ("Lisbon", 38.7223, -9.1393),
        "Qatar": ("Doha", 25.285447, 51.441883),
        "Romania": ("Bucharest", 44.4268, 26.1025),
        "Russia": ("Moscow", 55.7558, 37.6173),
        "Rwanda": ("Kigali", -1.9403, 29.8739),
        "Saint Kitts and Nevis": ("Basseterre", 17.3, -62.716667),
        "Saint Lucia": ("Castries", 14.0101, -61.0090),
        "Saint Vincent and the Grenadines": ("Kingstown", 13.16, -61.224),
        "Samoa": ("Apia", -13.8333, -171.7333),
        "San Marino": ("San Marino", 43.9424, 12.4578),
        "Sao Tome and Principe": ("São Tomé", 0.3365, 6.7273),
        "Saudi Arabia": ("Riyadh", 24.7136, 46.6753),
        "Senegal": ("Dakar", 14.7167, -17.4677),
        "Serbia": ("Belgrade", 44.7866, 20.4489),
        "Seychelles": ("Victoria", -4.6191, 55.4513),
        "Sierra Leone": ("Freetown", 8.4844, -13.2344),
        "Singapore": ("Singapore", 1.3521, 103.8198),
        "Slovakia": ("Bratislava", 48.1486, 17.1077),
        "Slovenia": ("Ljubljana", 46.0569, 14.5058),
        "Solomon Islands": ("Honiara", -9.433333, 159.95),
        "Somalia": ("Mogadishu", 2.03711, 45.34375),
        "South Africa": ("Pretoria", -25.7479, 28.2293),
        "South Sudan": ("Juba", 4.859363, 31.57125),
        "Spain": ("Madrid", 40.4168, -3.7038),
        "Sri Lanka": ("Sri Jayawardenepura Kotte", 6.9028, 79.8619),
        "Sudan": ("Khartoum", 15.5007, 32.5599),
        "Suriname": ("Paramaribo", 5.866667, -55.166667),
        "Sweden": ("Stockholm", 59.3293, 18.0686),
        "Switzerland": ("Bern", 46.9480, 7.4474),
        "Syria": ("Damascus", 33.5138, 36.2765),
        "Taiwan": ("Taipei", 25.0330, 121.5654),
        "Tajikistan": ("Dushanbe", 38.5598, 68.7870),
        "Tanzania": ("Dodoma", -6.1630, 35.7516),
        "Thailand": ("Bangkok", 13.7563, 100.5018),
        "Togo": ("Lomé", 6.1725, 1.2314),
        "Tonga": ("Nuku'alofa", -21.1394, -175.2018),
        "Trinidad and Tobago": ("Port of Spain", 10.666667, -61.516667),
        "Tunisia": ("Tunis", 36.8065, 10.1815),
        "Turkey": ("Ankara", 39.9208, 32.8541),
        "Turkmenistan": ("Ashgabat", 37.9601, 58.3261),
        "Tuvalu": ("Funafuti", -8.5211, 179.1962),
        "Uganda": ("Kampala", 0.316667, 32.5825),
        "Ukraine": ("Kyiv", 50.4501, 30.5234),
        "United Arab Emirates": ("Abu Dhabi", 24.4539, 54.3773),
        "United Kingdom": ("London", 51.5074, -0.1278),
        "United States": ("Washington, D.C.", 38.9072, -77.0369),
        "Uruguay": ("Montevideo", -34.9011, -56.1645),
        "Uzbekistan": ("Tashkent", 41.2995, 69.2401),
        "Vanuatu": ("Port Vila", -17.733333, 168.316667),
        "Vatican City": ("Vatican City", 41.9029, 12.4534),
        "Venezuela": ("Caracas", 10.4806, -66.9036),
        "Vietnam": ("Hanoi", 21.0285, 105.8542),
        "Yemen": ("Sana'a", 15.3547, 44.2064),
        "Zambia": ("Lusaka", -15.3875, 28.3228),
        "Zimbabwe": ("Harare", -17.8252, 31.0335)
    }

    # print("State, Capital, Lat (deg), Lon (deg), ECEF X (m), ECEF Y (m), ECEF Z (m)")
    with open("dataset/code/output/gs_loc.txt", "w") as f:
        f.write(f"{len(govs)}\n")
        for _, (capital, lat, lon) in govs.items():
            x, y, z = geodetic_to_ecef(lat, lon, alt=0)
            f.write(f"{x:.2f} {y:.2f} {z:.2f}\n")

if __name__ == "__main__":
    main()