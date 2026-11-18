"""Build the 5 kpc Gaia DR3 / Bailer-Jones geometric-distance teaching sample."""
import argparse
import csv
import hashlib
import heapq
import json
import math
from pathlib import Path
import re

import astropy.units as u
from astropy.coordinates import SkyCoord

from build_expanded_gaia_assets import corrected_excess, region
from gaia_archive import tiled_csv

ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = "https://gea.esac.esa.int/tap-server/tap/async"
TARGET = 24000
FIELD_LIMIT = 1000000
GAIA_SOURCE_COUNT = 1811709771
SELECTIONS = {
    "field": ("random_index<1000000", "", 500000),
    "bright": ("phot_g_mean_mag<10", "", 750000),
    "white-dwarfs": ("bp_rp<1.5", "AND g.phot_g_mean_mag-5*LOG10(d.r_med_geo)+5>9+3*g.bp_rp", 500000),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT / ".star-data-cache")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--refresh", action="store_true", help="Replace combined CSV caches with completed tiled results")
    args = parser.parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    template = (ROOT / "scripts/gaia-bailer-jones.adql").read_text()
    queries = {name: template.format(selection=selection, distance_selection=extra, limit=limit).strip()
               for name, (selection, extra, limit) in SELECTIONS.items()}

    def obtain(name):
        path = cache / f"gaia-bj-{name}.csv"
        if not path.exists() or args.refresh:
            if not args.download:
                raise ValueError(f"Missing {path}; run with --download")
            print(f"Downloading {name}", flush=True)
            selection, extra, limit = SELECTIONS[name]
            raw = tiled_csv(template, selection, extra, limit, cache, name)
            path.write_bytes(raw)
        return path

    paths = {name: obtain(name) for name in SELECTIONS}
    candidates, stats = {}, {}
    for name, path in paths.items():
        downloaded, accepted = 0, 0
        with path.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                downloaded += 1
                if any(not value.strip() for value in row.values()):
                    continue
                v = {key: float(value) for key, value in row.items() if key != "source_id"}
                if not all(math.isfinite(value) for value in v.values()):
                    continue
                lo, distance, hi = v["r_lo_geo"], v["r_med_geo"], v["r_hi_geo"]
                if not (0 < lo <= distance <= hi and distance <= 5000 and (hi-lo)/(2*distance) <= .2):
                    raise ValueError("Archive row violates distance-posterior selection")
                if not (v["parallax_error"] > 0 and v["ruwe"] < 1.4
                        and v["phot_bp_mean_flux_over_error"] > 10 and v["phot_rp_mean_flux_over_error"] > 10
                        and -1 <= v["bp_rp"] <= 6):
                    raise ValueError("Archive row violates quality selection")
                magnitude = v["phot_g_mean_mag"] - 5*math.log10(distance) + 5
                if name == "bright" and not v["phot_g_mean_mag"] < 10:
                    raise ValueError("Cached row does not match bright pool")
                if name == "white-dwarfs" and not (v["bp_rp"] < 1.5 and magnitude > 9+3*v["bp_rp"]):
                    raise ValueError("Cached row does not match white-dwarf pool")
                cstar, sigma = corrected_excess(v["bp_rp"], v["phot_bp_rp_excess_factor"], v["phot_g_mean_mag"])
                if abs(cstar) > 3*sigma or not -10 <= magnitude <= 18:
                    continue
                v.update(absoluteG=magnitude, cstar=cstar)
                candidates[row["source_id"]] = v
                accepted += 1
        stats[name] = {"downloaded": downloaded, "photometryAndRangeAccepted": accepted}

    def rank(item):
        source_id, v = item
        color, magnitude = v["bp_rp"], v["absoluteG"]
        inclusion = 1 if v["phot_g_mean_mag"] < 10 or region(color, magnitude) == "white-dwarf-region" else FIELD_LIMIT/GAIA_SOURCE_COUNT
        wd_boost = 200/(1+math.exp((8+2.5*color-magnitude)/.8))/(1+math.exp((color-1.8)/.3))
        importance = 1 + 3*10**(-.9*magnitude) + wd_boost
        mantissa = int.from_bytes(hashlib.sha256(source_id.encode("ascii")).digest()[:8], "big") >> 12
        return -math.log((mantissa+.5)/2**52)*inclusion/importance, int(source_id)

    chosen = sorted(heapq.nsmallest(TARGET, candidates.items(), key=rank), key=lambda item: int(item[0]))
    coordinates = SkyCoord(ra=[v["ra"] for _, v in chosen]*u.deg, dec=[v["dec"] for _, v in chosen]*u.deg,
                           distance=[v["r_med_geo"] for _, v in chosen]*u.pc, frame="icrs").galactic
    positions = coordinates.cartesian.xyz.to_value(u.pc).T
    stars, counts = [], {}
    for (source_id, v), xyz in zip(chosen, positions):
        label = region(v["bp_rp"], v["absoluteG"])
        counts[label] = counts.get(label, 0)+1
        stars.append([source_id, v["ra"], v["dec"], v["parallax"], v["parallax_error"],
                      *[round(float(value), 6) for value in xyz], v["bp_rp"],
                      None, None, None, None, False, v["absoluteG"], v["phot_g_mean_mag"], v["b"], v["cstar"],
                      v["r_med_geo"], v["r_lo_geo"], v["r_hi_geo"]])
    if len(stars) != TARGET or sum(row[18] > 2000 for row in stars) < 3000 or counts.get("white-dwarf-region", 0) < 1000:
        raise ValueError("Unexpectedly poor coverage of the deeper sample")
    data = {
        "version": 1, "title": "Gaia / Bailer-Jones geometric distances — out to 5,000 pc", "epoch": "J2016.0",
        "source": "ESA/Gaia/DPAC DR3 astrometry and photometry; Bailer-Jones et al. (2021) EDR3 geometric distances",
        "distanceUrl": "https://doi.org/10.3847/1538-3881/abd806",
        "distanceCatalogUrl": "https://dc.g-vo.org/tableinfo/gedr3dist.main",
        "distanceLicense": "CC BY 4.0; https://creativecommons.org/licenses/by/4.0/",
        "endpoint": ENDPOINT, "queries": queries, "selectionStats": stats, "regionCounts": counts,
        "queryExecution": "Completed asynchronous jobs in disjoint source-ID tiles covering [0, 12*2**59); timed-out tiles are bisected and retried, never accepted as partial CSV; each tile must stay below its row ceiling",
        "selection": "0 < geometric posterior median <= 5000 pc; positive ordered 16th/84th percentiles; (r_hi_geo-r_lo_geo)/(2*r_med_geo) <= 0.2; RUWE < 1.4, visibility periods > 8, G flux SNR > 50, BP/RP flux SNR > 10; corrected photometric excess within 3 sigma; no parallax-SNR or positivity cut",
        "distanceMethod": "Bailer-Jones et al. (2021) EDR3 geometric posterior median r_med_geo in pc; 16th/84th percentiles r_lo_geo/r_hi_geo; position-dependent Galactic prior and corrected EDR3 parallaxes; not photogeometric distances and not inverse parallax",
        "distanceRangePc": [0, 5000], "radiusPc": 5000,
        "coordinateSystem": "Solar-barycentric Galactic XYZ, J2016.0; Astropy ICRS-to-Galactic transform at the published geometric median distance",
        "magnitudeMethod": "M_G = G - 5*log10(r_med_geo_pc) + 5; no extinction correction; distance-only interval uses the same distance percentiles",
        "extinction": "Observed BP-RP and absolute G without extinction or reddening correction",
        "photometryQuality": "abs(corrected BP/RP excess C*) <= 3 sigma, Riello et al. 2021 Table 2 and Eq. 18; sigma evaluated at G>=4",
        "photometryQualityUrl": "https://doi.org/10.1051/0004-6361/202039587",
        "sampling": {
            "candidateUnique": len(candidates), "target": TARGET, "fieldRandomIndexLimit": FIELD_LIMIT,
            "gaiaSourceCount": GAIA_SOURCE_COUNT,
            "method": "Keep the smallest -log(U) * inclusion / importance; U is the 52-bit midpoint of SHA-256(source_id); inclusion is 1 for G<10 or the full white-dwarf query, otherwise 1000000/1811709771",
            "importance": "1 + 3*10**(-0.9*M_G) + 200/(1+exp((8+2.5*BP_RP-M_G)/0.8))/(1+exp((BP_RP-1.8)/0.3))",
            "interpretation": "Smooth importance sampling without replacement, with bright and white-dwarf enrichment; not representative population frequencies or a complete density census",
        },
        "completeness": "Cumulative volume including the inner 2 kpc, without an extra apparent-magnitude cutoff on the combined sample; detection, photometry and distance-posterior quality limits remain; rare HR regions are emphasized; region labels are not confirmed stellar classifications",
        "columns": ["sourceId", "ra", "dec", "parallax", "parallaxError", "xPc", "yPc", "zPc", "bpRp",
                    "luminosity", "temperature", "luminosityLower", "luminosityUpper", "luminosityIntervalValid",
                    "absoluteG", "apparentG", "galacticLatitude", "correctedExcess", "distanceMedianPc", "distanceLowerPc", "distanceUpperPc"],
        "colorRange": [-1, 6], "magnitudeRange": [-10, 18], "stars": stars,
    }
    serialized = json.dumps(data, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    output = "assets/gaia-bailer-jones-catalog.json"
    (ROOT / output).write_text(serialized+"\n", encoding="utf-8")
    html_path = ROOT / "hrd.html"
    html = html_path.read_text(encoding="utf-8")
    tag = '<script type="application/json" id="wide-catalog">'+serialized+'</script>'
    pattern = r'<script type="application/json" id="wide-catalog">.*?</script>'
    if re.search(pattern, html):
        html = re.sub(pattern, lambda _: tag, html, count=1)
    else:
        marker = '  <script type="application/json" id="hyades-catalog">'
        if marker not in html:
            raise ValueError("Cannot locate catalog insertion point")
        html = html.replace(marker, '  '+tag+'\n'+marker, 1)
    html_path.write_text(html, encoding="utf-8")
    manifest_path = ROOT / "assets/source-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for name, path in paths.items():
        manifest[path.name] = {"endpoint": ENDPOINT, "queryIn": output, "component": name,
                               "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    manifest_path.write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps({"stars": len(stars), "beyond2000pc": sum(row[18]>2000 for row in stars),
                      "regions": counts, "selections": stats}, indent=2))


if __name__ == "__main__":
    main()
