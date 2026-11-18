"""Build a real Hyades cluster example from GCNS membership and Gaia DR3."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parents[1]
MEMBERS_ENDPOINT = "https://dc.g-vo.org/tap/sync"
GAIA_ENDPOINT = "https://gea.esac.esa.int/tap-server/tap/sync"
MEMBERS_QUERY = "SELECT * FROM gcns.hyacob"


def fetch(path, endpoint, query, download, post=False):
    if path.exists():
        return
    if not download:
        raise SystemExit(f"Missing {path}; use --download")
    params = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query})
    request = Request(endpoint, data=params.encode()) if post else endpoint + "?" + params
    with urlopen(request, timeout=120) as response:
        data = response.read()
    if not data.startswith(b"source_id,"):
        raise ValueError("Unexpected archive response")
    path.write_bytes(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT / ".star-data-cache")
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    member_path, gaia_path = cache / "hyades-members.csv", cache / "hyades-gaia.csv"
    fetch(member_path, MEMBERS_ENDPOINT, MEMBERS_QUERY, args.download)
    members = {r["source_id"]: r for r in csv.DictReader(member_path.open()) if r["name"] == "Hyades"}
    if len(members) != 920:
        raise ValueError("Unexpected GCNS Hyades membership count")
    query = """SELECT g.source_id, g.ra, g.dec, g.parallax, g.parallax_error,
g.phot_g_mean_mag, g.bp_rp, g.ruwe, ap.teff_gspphot,
ap.lum_flame, ap.lum_flame_lower, ap.lum_flame_upper
FROM gaiadr3.gaia_source AS g
LEFT OUTER JOIN gaiadr3.astrophysical_parameters AS ap ON g.source_id=ap.source_id
WHERE g.source_id IN (""" + ",".join(members) + ")"
    fetch(gaia_path, GAIA_ENDPOINT, query, args.download, post=True)
    chosen = []
    for row in csv.DictReader(gaia_path.open()):
        member = members[row["source_id"]]
        values = {key: float(value) if value.strip() else None for key, value in row.items() if key != "source_id"}
        if float(member["dist_cluster"]) > 20:
            continue
        if any(values[key] is None or not math.isfinite(values[key]) for key in ["ra", "dec", "parallax", "parallax_error", "phot_g_mean_mag", "bp_rp"]):
            continue
        if not (values["parallax"] > 0 and values["parallax_error"] > 0 and values["parallax"] / values["parallax_error"] > 50 and -.5 <= values["bp_rp"] <= 5.5):
            continue
        chosen.append((row["source_id"], values, member))
    chosen.sort(key=lambda row: int(row[0]))
    if len(chosen) < 100:
        raise ValueError("Unexpectedly small cluster sample")
    positions = SkyCoord(ra=[r[1]["ra"] for r in chosen] * u.deg,
                         dec=[r[1]["dec"] for r in chosen] * u.deg,
                         distance=[1000/r[1]["parallax"] for r in chosen] * u.pc,
                         frame="icrs").galactic.cartesian.xyz.to_value(u.pc).T
    center = positions.mean(axis=0)
    radius = math.ceil(float(np.abs(positions - center).max()) / 5) * 5
    stars = []
    for (source_id, v, member), xyz in zip(chosen, positions):
        lower, upper, luminosity = (v[key] for key in ["lum_flame_lower", "lum_flame_upper", "lum_flame"])
        valid_interval = all(x is not None and math.isfinite(x) for x in [lower, upper, luminosity]) and 0 < lower <= luminosity <= upper
        magnitude = v["phot_g_mean_mag"] + 5 * math.log10(v["parallax"]) - 10
        stars.append([source_id, v["ra"], v["dec"], v["parallax"], v["parallax_error"],
                      *[round(float(x), 6) for x in xyz], v["bp_rp"], luminosity,
                      v["teff_gspphot"], lower, upper, valid_interval, magnitude,
                      v["phot_g_mean_mag"], float(member["dist_cluster"]), float(member["c_cluster"])])
    asset = {
        "version": 1, "title": "Hyades — inner 20 pc", "epoch": "J2016.0",
        "source": "Gaia Collaboration (2021), GCNS candidate memberships; ESA/Gaia/DPAC Gaia DR3",
        "membershipUrl": "https://dc.zah.uni-heidelberg.de/tableinfo/gcns.hyacob",
        "membershipLicense": "CC BY-NC 4.0, Richard Smart (2023)",
        "membershipQuery": MEMBERS_QUERY, "gaiaQuery": query,
        "publishedMembers": len(members), "publishedInnerMembers": sum(float(r["dist_cluster"]) <= 20 for r in members.values()),
        "luminosityCount": sum(r[9] is not None and r[9] > 0 for r in stars),
        "columns": ["sourceId", "ra", "dec", "parallax", "parallaxError", "xPc", "yPc", "zPc", "bpRp", "luminosity", "temperature", "luminosityLower", "luminosityUpper", "luminosityIntervalValid", "absoluteG", "apparentG", "publishedClusterDistance", "membershipStatistic"],
        "centerPc": center.tolist(), "radiusPc": radius,
        "coordinateSystem": "Solar-barycentric Galactic Cartesian, displayed relative to the sample mean",
        "colorRange": [-.5, 5], "magnitudeRange": [-1, 18], "stars": stars,
    }
    serialized = json.dumps(asset, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    (ROOT / "assets/hyades-catalog.json").write_text(serialized + "\n")
    html_path = ROOT / "hrd.html"
    html = html_path.read_text(encoding="utf-8")
    tag = '<script type="application/json" id="hyades-catalog">' + serialized + '</script>'
    pattern = r'<script type="application/json" id="hyades-catalog">.*?</script>'
    html = re.sub(pattern, lambda _: tag, html, count=1) if re.search(pattern, html) else html.replace('  <script>\n', '  ' + tag + '\n  <script>\n', 1)
    html_path.write_text(html, encoding="utf-8")
    manifest_path = ROOT / "assets/source-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for path, endpoint in [(member_path, MEMBERS_ENDPOINT), (gaia_path, GAIA_ENDPOINT)]:
        manifest[path.name] = {"endpoint": endpoint, "queryIn": "assets/hyades-catalog.json", "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Built {len(stars)} Hyades candidates; {asset['luminosityCount']} have luminosity estimates; coordinate half-width {radius} pc")


if __name__ == "__main__":
    main()
