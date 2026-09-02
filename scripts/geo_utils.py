#!/usr/bin/env python3
"""Geography helpers for the metadata master: parse INSDC `lat_lon` strings to decimal, normalise
country tokens, and geocode a country name to a centroid. Country centroids are computed **empirically
from our own data** (median of the real coordinates we already hold per country) in
`build_metadata_master.py`; this module supplies the parser, the token normaliser, and a small
hand-curated fallback centroid table for the ~tail countries that carry no explicit coordinate in our
data (so a place name still yields coordinates). No network, no external library.
"""
import re

# tokens that are not a real place -> no geography
NULL_TOKENS = {"", "missing", "not applicable", "not determined", "not available", "not collected",
               "na", "n/a", "unknown", "none", "unspecified", "not provided"}

# raw geo_loc_name country token -> canonical country (historical names, leaked Chinese city names, case)
ALIASES = {
    "ussr": "Russia", "czechoslovakia": "Czech Republic", "serbia and montenegro": "Serbia",
    "korea": "South Korea",
    "xi'an": "China", "xianyang": "China", "hanzhong": "China", "lintong": "China",
    "weinan": "China", "baoji": "China",
}

# fallback centroids (lat, lng) for countries/territories that may lack an explicit coordinate in our
# data. Only used when the empirical centroid is absent. Coarse country-level points.
FALLBACK_CENTROIDS = {
    "Falkland Islands (Islas Malvinas)": (-51.8, -59.5), "Uruguay": (-32.5, -55.8),
    "State of Palestine": (31.9, 35.2), "Serbia": (44.0, 21.0), "Somalia": (5.15, 46.2),
    "Afghanistan": (33.9, 67.7), "Mauritania": (20.0, -10.9), "Kuwait": (29.3, 47.5),
    "Morocco": (31.8, -7.1), "Trinidad and Tobago": (10.7, -61.2), "Zambia": (-13.1, 27.85),
    "Armenia": (40.1, 45.0), "Sudan": (12.9, 30.2), "Syria": (34.8, 38.99), "Kyrgyzstan": (41.2, 74.8),
    "Comoros": (-11.6, 43.3), "Madagascar": (-18.8, 46.9), "Rwanda": (-1.9, 29.9),
    "Solomon Islands": (-9.6, 160.2), "Ethiopia": (9.15, 40.5), "Svalbard": (77.9, 20.9),
    "Republic of the Congo": (-0.7, 15.0), "Haiti": (19.1, -72.3), "Azerbaijan": (40.1, 47.6),
    "Burundi": (-3.4, 29.9), "Cameroon": (5.7, 12.35), "Curacao": (12.2, -69.0),
    "Zimbabwe": (-19.0, 29.15), "Libya": (26.3, 17.2), "Vanuatu": (-16.5, 168.0),
    "Bahamas": (24.5, -76.5), "Micronesia": (6.9, 158.2), "Bahrain": (26.0, 50.55),
    "Samoa": (-13.76, -172.1), "Seychelles": (-4.7, 55.5), "Kiribati": (1.87, -157.4),
    "Papua New Guinea": (-6.3, 144.0), "Cote d'Ivoire": (7.5, -5.5), "Albania": (41.15, 20.2),
    "El Salvador": (13.8, -88.9), "Estonia": (58.6, 25.5), "Greenland": (71.7, -42.6),
    "Dominican Republic": (18.7, -70.2), "Niger": (17.6, 8.1), "South Sudan": (7.9, 29.7),
    "Central African Republic": (6.6, 20.9), "Faroe Islands": (62.0, -6.9), "Barbados": (13.2, -59.5),
    "Grenada": (12.1, -61.7), "Turkmenistan": (38.97, 59.6), "Guatemala": (15.8, -90.2),
    "Gabon": (-0.8, 11.6), "Eswatini": (-26.5, 31.5), "Palau": (7.5, 134.6),
    "Democratic Republic of the Congo": (-4.0, 21.8), "Oman": (21.5, 55.9), "Venezuela": (6.4, -66.6),
    "Mongolia": (46.9, 103.8), "Mayotte": (-12.8, 45.2), "Reunion": (-21.1, 55.5),
    "Sierra Leone": (8.5, -11.8), "Guinea-Bissau": (11.8, -15.2), "Botswana": (-22.3, 24.7),
    "Chad": (15.5, 18.7), "Malawi": (-13.25, 34.3), "Mauritius": (-20.3, 57.6),
    "Puerto Rico": (18.2, -66.5), "Honduras": (15.2, -86.2), "Benin": (9.3, 2.3),
    "Nicaragua": (12.9, -85.2), "Gambia": (13.4, -15.3),
}


def parse_latlon(s):
    """INSDC lat_lon string -> (lat, lng) decimal, or None. Handles 'D N D E' and 'lat, lng'."""
    s = (s or "").strip()
    if not s:
        return None
    m = re.match(r"^\s*([-\d.]+)\s*([NS])\s*[, ]\s*([-\d.]+)\s*([EW])\s*$", s, re.I)
    if m:
        lat, lng = float(m.group(1)), float(m.group(3))
        if m.group(2).upper() == "S":
            lat = -lat
        if m.group(4).upper() == "W":
            lng = -lng
        if -90 <= lat <= 90 and -180 <= lng <= 180:
            return (lat, lng)
    m = re.match(r"^\s*([-\d.]+)\s*,\s*([-\d.]+)\s*$", s)
    if m:
        try:
            lat, lng = float(m.group(1)), float(m.group(2))
            if -90 <= lat <= 90 and -180 <= lng <= 180:
                return (lat, lng)
        except ValueError:
            pass
    return None


def split_geo_name(geo):
    """'Country: Region, City' -> (country_canonical, admin_remainder) or (None, '') if null token."""
    geo = (geo or "").strip()
    if not geo:
        return None, ""
    parts = geo.split(":", 1)
    country = parts[0].strip()
    admin = parts[1].strip() if len(parts) > 1 else ""
    if country.lower() in NULL_TOKENS:
        return None, ""
    country = ALIASES.get(country.lower(), country)
    return country, admin
