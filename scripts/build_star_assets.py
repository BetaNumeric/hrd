"""Build the portable catalog and empirical-PSF assets used by hrd.html.

Requires numpy and astropy. Public source files can be downloaded with --download;
normal rebuilds use the existing cache and make no network requests.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
from pathlib import Path
import re
from urllib.request import urlopen

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
CATALOG_URL = (
    "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
    "-source=J/ApJ/710/1032/pm_stars&-out=F,XREF,YREF,Bmag,Rmag,Q"
    "&-out.max=110000&F=c"
)
PSF_ROOT = (
    "https://www.stsci.edu/files/live/sites/www/files/home/hst/"
    "instrumentation/acs/data-analysis/psf/_documents/"
)
SOURCES = {
    "omega-pm-all.tsv": CATALOG_URL,
    "hst-psf-f435w.fits": PSF_ROOT + "STDPSF_ACSWFC_F435W_SM3.fits",
    "hst-psf-f625w.fits": PSF_ROOT + "STDPSF_ACSWFC_F625W_SM3.fits",
}


def catalog_asset(path: Path) -> dict:
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith("c\t"):
            continue
        fields = line.split("\t")
        x, y, blue, red = map(float, fields[1:5])
        color = blue - red
        # These cuts discard unmeasured-band sentinels and extreme outliers.
        # The source is an observed proper-motion sample, not a complete census.
        if 16 <= blue <= 25.5 and 15.5 <= red <= 24 and -0.5 <= color <= 2.9:
            rows.append([x, y, blue, red])
    if len(rows) < 10000:
        raise ValueError("Catalog response is incomplete or is an archive error page")

    # Preserve the sparse blue/bright sources and select the remaining stars by
    # a stable hash of their actual reference-frame positions. Never move a
    # measurement onto an idealized evolutionary curve.
    rare, ordinary = [], []
    for row in rows:
        (rare if row[3] < 17 or row[2] - row[3] < 0.5 else ordinary).append(row)
    ordinary.sort(key=lambda row: hashlib.sha256(
        f"{row[0]:.2f},{row[1]:.2f}".encode("ascii")).digest())
    chosen = rare + ordinary[:10000 - len(rare)]
    chosen.sort(key=lambda row: (row[0], row[1]))
    return {
        "version": 1,
        "title": "Omega Centauri",
        "source": "Anderson & van der Marel (2010), ApJ 710, 1032",
        "doi": "10.26093/cds/vizier.17101032",
        "table": "J/ApJ/710/1032/pm_stars (central field)",
        "url": CATALOG_URL,
        "columns": ["xRef", "yRef", "f435w", "f625w"],
        "candidateCount": len(rows),
        "retainedRareCount": len(rare),
        "selection": "Retain F625W < 17 or F435W-F625W < 0.5; hash-sample the rest to 10,000 stars.",
        "colorRange": [-0.5, 3],
        "magnitudeRange": [15, 24.5],
        "stars": chosen,
    }


def psf_asset(cache: Path) -> dict:
    planes = [20, 33, 57, 69]
    size = 101
    result = []
    for plane in planes:
        bands = []
        for name in ["hst-psf-f435w.fits", "hst-psf-f625w.fits"]:
            # Retain every native sample of the measured, 4x supersampled ePSF.
            image = fits.getdata(cache / name)[plane].astype(np.float64)
            if image.shape != (size, size):
                raise ValueError(f"Unexpected native PSF dimensions: {image.shape}")
            encoded = np.rint(np.clip((np.log10(np.maximum(image, 1e-8)) + 8) / 8, 0, 1) * 255).astype(np.uint8)
            bands.append(encoded)
        result.append(np.stack(bands, axis=-1))
    packed = np.stack(result).tobytes()
    return {
        "version": 2, "size": size, "count": len(planes),
        "sampling": "Native 101x101 samples of the 4x-supersampled empirical profiles; no spatial downsampling",
        "bands": ["F435W", "F625W"], "planes": planes,
        "source": "STScI ACS/WFC SM3 empirical STDPSF library",
        "url": "https://www.stsci.edu/hst/instrumentation/acs/data-analysis/psf",
        "encoding": "uint8 interleaved blue/red; intensity = 10**(8*byte/255-8); zero byte = zero",
        "data": base64.b64encode(packed).decode("ascii"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT / ".star-data-cache")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--psf-only", action="store_true", help="Rebuild the star profiles without downloading or changing the legacy catalog")
    args = parser.parse_args()
    cache = args.cache_dir.resolve()
    cache.mkdir(parents=True, exist_ok=True)
    for name, url in SOURCES.items():
        if args.psf_only and not name.startswith("hst-psf-"):
            continue
        path = cache / name
        if not path.exists():
            if not args.download:
                raise SystemExit(f"Missing {path}; run with --download to fetch public source data")
            with urlopen(url, timeout=60) as response:
                data = response.read()
            if len(data) < 10000:
                raise ValueError(f"Incomplete source response: {url}")
            path.write_bytes(data)
    assets = {"empirical-psf": ("hubble-acs-psf.json", psf_asset(cache))}
    if not args.psf_only:
        assets["observed-catalog"] = ("omega-centauri-catalog.json", catalog_asset(cache / "omega-pm-all.tsv"))
    html_path = ROOT / "hrd.html"
    html = html_path.read_text(encoding="utf-8")
    for element_id, (name, data) in assets.items():
        serialized = json.dumps(data, separators=(",", ":"), ensure_ascii=True)
        (ROOT / "assets" / name).write_text(serialized + "\n", encoding="utf-8")
        tag = f'<script type="application/json" id="{element_id}">{serialized}</script>'
        pattern = rf'<script type="application/json" id="{element_id}">.*?</script>'
        if re.search(pattern, html):
            html = re.sub(pattern, lambda _: tag, html, count=1)
        else:
            html = html.replace('  <script>\n', f'  {tag}\n  <script>\n', 1)
    html_path.write_text(html, encoding="utf-8")
    manifest_path = ROOT / "assets" / "source-manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.update({name: {"url": url, "sha256": hashlib.sha256((cache / name).read_bytes()).hexdigest()}
                     for name, url in SOURCES.items() if not args.psf_only or name.startswith("hst-psf-")})
    (ROOT / "assets" / "source-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Built four native 101x101 two-band empirical profiles" + ("" if args.psf_only else f" and {len(assets['observed-catalog'][1]['stars']):,} catalog stars"))


if __name__ == "__main__":
    main()
