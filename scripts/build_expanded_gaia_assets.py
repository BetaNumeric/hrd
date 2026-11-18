"""Build separate nearby and distant Gaia DR3 teaching samples from real photometry."""
from concurrent.futures import ThreadPoolExecutor
import argparse
import csv
import hashlib
import heapq
import json
import math
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import urlopen

import astropy.units as u
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = "https://gea.esac.esa.int/tap-server/tap/sync"
MAX_ROWS = 500000
SELECTIONS = {
    "local": "g.parallax>10",
    "bright": "g.phot_g_mean_mag<9",  # Nearby reproducibility and distant bright enrichment.
    "distant-field": "g.parallax<=10 AND g.random_index<5000000",
    "distant-white-dwarfs": "g.parallax<=10 AND g.bp_rp<1.5 AND g.phot_g_mean_mag+5*LOG10(g.parallax)-10>9+3*g.bp_rp",
    "distant-luminous": "g.parallax<=10 AND g.phot_g_mean_mag+5*LOG10(g.parallax)-10< -4",
}
REGION_TARGETS = {"white-dwarf-region": 3500, "giant-region": 5000, "luminous-region": 1500}
FIELD_TARGET = 16000
DISTANT_TARGET = 24000
GAIA_SOURCE_COUNT = 1811709771
FIELD_RANDOM_INDEX_LIMIT = 5000000


def corrected_excess(color, excess, magnitude):
    """Riello et al. (2021), Table 2 and Eq. 18; DR3 uses EDR3 photometry."""
    if color < .5:
        expected = 1.154360 + .033772 * color + .032277 * color ** 2
    elif color < 4:
        expected = 1.162004 + .011464 * color + .049255 * color ** 2 - .005879 * color ** 3
    else:
        expected = 1.057572 + .140537 * color
    sigma = .0059898 + 8.817481e-12 * max(magnitude, 4) ** 7.618399
    return excess - expected, sigma


def region(color, magnitude):
    if color < 1.5 and magnitude > 9 + 3 * color:
        return "white-dwarf-region"
    if magnitude < -4:
        return "luminous-region"
    if color > .8 and magnitude < 3:
        return "giant-region"
    if -2 <= magnitude - (3.2 * color + 2) <= 2.5:
        return "main-sequence-region"
    return "other"


def build_rows(chosen):
    """Use the measured sky positions, inverse parallaxes and photometry unchanged."""
    coordinates = SkyCoord(ra=[v["ra"] for _, v in chosen] * u.deg,
                           dec=[v["dec"] for _, v in chosen] * u.deg,
                           distance=[1000/v["parallax"] for _, v in chosen] * u.pc,
                           frame="icrs").galactic
    positions = coordinates.cartesian.xyz.to_value(u.pc).T
    stars, counts = [], {}
    for (source_id, v), xyz in zip(chosen, positions):
        label = region(v["bp_rp"], v["absoluteG"])
        counts[label] = counts.get(label, 0) + 1
        stars.append([source_id, v["ra"], v["dec"], v["parallax"], v["parallax_error"],
                      *[round(float(x), 6) for x in xyz], v["bp_rp"],
                      None, None, None, None, False, v["absoluteG"],
                      v["phot_g_mean_mag"], v["b"], v["cstar"]])
    return stars, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT / ".star-data-cache")
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    template = (ROOT / "scripts/gaia-expanded.adql").read_text()
    queries = {name: template.format(limit=MAX_ROWS, selection=selection).strip() for name, selection in SELECTIONS.items()}

    def obtain(name):
        path = cache / f"gaia-expanded-{name}.csv"
        if not path.exists():
            if not args.download:
                raise ValueError(f"Missing {path}; run with --download")
            print(f"Downloading {name}", flush=True)
            url = ENDPOINT + "?" + urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": MAX_ROWS, "QUERY": queries[name]})
            with urlopen(url, timeout=120) as response:
                raw = response.read()
            if not raw.startswith(b"source_id,"):
                raise ValueError(f"Unexpected archive response for {name}")
            path.write_bytes(raw)
        return path

    with ThreadPoolExecutor(max_workers=2) as pool:
        paths = dict(zip(SELECTIONS, pool.map(obtain, SELECTIONS)))
    pools, stats = {}, {}
    for name, path in paths.items():
        pool = {}
        downloaded, accepted = 0, 0
        for row in csv.DictReader(path.open(encoding="utf-8", newline="")):
            downloaded += 1
            if any(not value.strip() for key, value in row.items() if key != "source_id"):
                continue
            values = {key: float(value) for key, value in row.items() if key != "source_id"}
            if not all(math.isfinite(v) for v in values.values()):
                continue
            p, error = values["parallax"], values["parallax_error"]
            if not (p >= .5 and error > 0 and p / error > 10 and values["ruwe"] < 1.4):
                raise ValueError("Archive row violates astrometric selection")
            if not (p > 10 if name == "local" else values["phot_g_mean_mag"] < 9 if name == "bright" else p <= 10):
                raise ValueError("Cached row does not match its pool selection")
            magnitude = values["phot_g_mean_mag"] + 5 * math.log10(p) - 10
            if not -10 <= magnitude <= 18:
                continue
            cstar, sigma = corrected_excess(values["bp_rp"], values["phot_bp_rp_excess_factor"], values["phot_g_mean_mag"])
            if abs(cstar) > 3 * sigma:
                continue
            values["cstar"] = cstar
            values["absoluteG"] = magnitude
            pool[row["source_id"]] = values
            accepted += 1
        if downloaded >= MAX_ROWS:
            raise ValueError(f"{name} reached the archive row limit; increase MAX_ROWS before sampling")
        stats[name] = {"downloaded": downloaded, "photometryAndRangeAccepted": accepted}
        pools[name] = pool
    # Preserve the existing nearby sample exactly.
    candidates = {**pools["local"], **pools["bright"]}
    # Hashing the exact 64-bit ID makes selection independent of archive row order.
    ordered = sorted(candidates, key=lambda source_id: hashlib.sha256(source_id.encode("ascii")).digest())
    chosen = {source_id: candidates[source_id] for source_id in ordered[:FIELD_TARGET]}
    region_retained = {label: 0 for label in REGION_TARGETS}
    for source_id in ordered:
        values = candidates[source_id]
        label = region(values["bp_rp"], values["absoluteG"])
        if label in REGION_TARGETS and region_retained[label] < REGION_TARGETS[label]:
            chosen[source_id] = values
            region_retained[label] += 1
    chosen = sorted(chosen.items(), key=lambda item: int(item[0]))
    if len(chosen) < 15000:
        raise ValueError("Expanded sample is unexpectedly small")
    stars, counts = build_rows(chosen)
    # Smooth weights avoid artificial horizontal steps from hard region quotas.
    field = pools["distant-field"]
    distant_candidates = {
        **field,
        **{source_id: v for source_id, v in pools["bright"].items() if v["parallax"] <= 10},
        **pools["distant-white-dwarfs"],
        **pools["distant-luminous"],
    }
    field_fraction = FIELD_RANDOM_INDEX_LIMIT / GAIA_SOURCE_COUNT

    def weighted_rank(item):
        source_id, v = item
        color, magnitude = v["bp_rp"], v["absoluteG"]
        full_pool = v["phot_g_mean_mag"] < 9 or region(color, magnitude) == "white-dwarf-region"
        inclusion = 1.0 if full_pool else field_fraction
        wd_boost = 200 / (1 + math.exp((8 + 2.5 * color - magnitude) / .8)) / (1 + math.exp((color - 1.8) / .3))
        importance = 1 + 3 * 10 ** (-.9 * magnitude) + wd_boost
        # A 52-bit midpoint is deterministic and strictly within (0,1).
        mantissa = int.from_bytes(hashlib.sha256(source_id.encode("ascii")).digest()[:8], "big") >> 12
        uniform = (mantissa + .5) / 2 ** 52
        return -math.log(uniform) * inclusion / importance, int(source_id)

    distant_chosen = heapq.nsmallest(DISTANT_TARGET, distant_candidates.items(), key=weighted_rank)
    distant_rows, distant_counts = build_rows(sorted(distant_chosen, key=lambda item: int(item[0])))
    if distant_counts.get("white-dwarf-region", 0) < 1000 or sum(row[14] > 4 for row in distant_rows) < 5000:
        raise ValueError("Deeper sample does not adequately cover the faint diagram")
    asset = {
        "version": 2, "title": "Expanded Gaia DR3 sample", "epoch": "J2016.0",
        "source": "ESA/Gaia/DPAC Gaia DR3, gaia_source",
        "endpoint": ENDPOINT, "queries": queries, "selectionStats": stats, "regionCounts": counts,
        "selection": "Union of local stars (<100 pc) and bright stars (G<9, up to 2 kpc), quality filtered; SHA-256 source-ID sampling of 16000 field stars plus extra white-dwarf, giant and luminous CMD regions; deduplicated; not a census or confirmed type classification",
        "sampling": {"eligibleUnique": len(candidates), "fieldTarget": FIELD_TARGET, "extraRegionTargets": REGION_TARGETS, "extraRegionRetained": region_retained},
        "columns": ["sourceId", "ra", "dec", "parallax", "parallaxError", "xPc", "yPc", "zPc", "bpRp", "luminosity", "temperature", "luminosityLower", "luminosityUpper", "luminosityIntervalValid", "absoluteG", "apparentG", "galacticLatitude", "correctedExcess"],
        "coordinateSystem": "Solar-barycentric Galactic XYZ, J2016.0",
        "distanceMethod": "1000/parallax_mas; positive parallaxes with SNR>10; no zero-point correction; formal errors exclude systematics",
        "extinction": "No extinction correction; color and absolute G are observed values affected by dust",
        "photometryQuality": "abs(corrected BP/RP excess C*) <= 3 sigma, Riello et al. 2021 Table 2 and Eq. 18; sigma evaluated at G>=4",
        "photometryQualityUrl": "https://doi.org/10.1051/0004-6361/202039587",
        "radiusPc": 2000, "colorRange": [-1, 6], "magnitudeRange": [-10, 18],
        "stars": stars,
    }
    html_path = ROOT / "hrd.html"
    html = html_path.read_text(encoding="utf-8")
    asset["parentSelection"] = asset.pop("selection")
    asset["parentSampleCount"] = len(stars)
    partitions = {
        "nearby": {"title": "Nearby Gaia DR3 stars — within 100 pc", "radiusPc": 100,
                   "distanceRangePc": [0, 100], "selection": "parallax > 10 mas (distance < 100 pc); no additional apparent G magnitude cut",
                   "stars": [row for row in stars if row[3] > 10]},
        "distant": {"title": "Distant Gaia DR3 stars — 100–2,000 pc", "radiusPc": 2000,
                    "distanceRangePc": [100, 2000], "selection": "0.5 <= parallax <= 10 mas (100–2,000 pc); no additional apparent-magnitude limit",
                    "stars": distant_rows},
    }
    outputs, report = [], {}
    for name, partition in partitions.items():
        region_counts = {}
        for row in partition["stars"]:
            label = region(row[8], row[14])
            region_counts[label] = region_counts.get(label, 0) + 1
        data = {**asset, **partition, "version": 3, "partition": name, "regionCounts": region_counts}
        query_names = ("local", "bright") if name == "nearby" else ("distant-field", "bright", "distant-white-dwarfs", "distant-luminous")
        data["queries"] = {key: queries[key] for key in query_names}
        data["selectionStats"] = {key: stats[key] for key in query_names}
        if name == "distant":
            data.update({
                "version": 5, "apparentMagnitudeCut": None,
                "parentSelection": "Same astrometric and photometric quality criteria as nearby stars; full-depth all-sky field cohort enriched with bright and white-dwarf query pools over 100–2,000 pc; no additional apparent-magnitude cutoff on the combined sample; smooth importance sampling without replacement",
                "parentSampleCount": len(distant_rows),
                "sampling": {
                    "candidateUnique": len(distant_candidates), "fieldEligible": len(field),
                    "target": DISTANT_TARGET, "fieldRandomIndexLimit": FIELD_RANDOM_INDEX_LIMIT,
                    "gaiaSourceCount": GAIA_SOURCE_COUNT, "fieldInclusionFraction": field_fraction,
                    "method": "Keep the smallest -log(U) * candidateInclusion / importance; U is the 52-bit midpoint of SHA-256(source_id). Candidate inclusion is 1 for G<9 or the full WD-region query, otherwise fieldRandomIndexLimit/gaiaSourceCount.",
                    "importance": "1 + 3*10**(-0.9*M_G) + 200/(1+exp((8+2.5*BP_RP-M_G)/0.8))/(1+exp((BP_RP-1.8)/0.3))",
                    "interpretation": "Smooth weights emphasize luminous and white-dwarf regions without hard display-group quotas; inclusion factors reduce enrichment-pool overrepresentation; this is not a statistically representative density sample",
                },
                "completeness": "A teaching sample with extra stars from sparse diagram regions; Gaia detection and quality limits remain; not a complete census or confirmed type classification",
            })
        serialized = json.dumps(data, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
        filename = f"assets/gaia-{name}-catalog.json"
        (ROOT / filename).write_text(serialized + "\n", encoding="utf-8")
        outputs.append(filename)
        element_id = "nearby-catalog" if name == "nearby" else "distant-catalog"
        tag = f'<script type="application/json" id="{element_id}">' + serialized + '</script>'
        pattern = rf'<script type="application/json" id="{element_id}">.*?</script>'
        if re.search(pattern, html):
            html = re.sub(pattern, lambda _: tag, html, count=1)
        else:
            marker = '  <script type="application/json" id="hyades-catalog">'
            if marker not in html:
                raise ValueError("Cannot locate the catalog insertion point in hrd.html")
            html = html.replace(marker, '  ' + tag + '\n' + marker, 1)
        report[name] = {"stars": len(partition["stars"]), "regions": region_counts}
    html_path.write_text(html, encoding="utf-8")
    manifest_path = ROOT / "assets/source-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for key in list(manifest):
        if key.startswith("gaia-nearby"):
            del manifest[key]
    for name, path in paths.items():
        output = outputs[0] if name in ("local", "bright") else outputs[1]
        manifest[path.name] = {"endpoint": ENDPOINT, "queryIn": output, "outputs": outputs if name == "bright" else [output], "component": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"partitions": report, "selections": stats}, indent=2))


if __name__ == "__main__":
    main()
