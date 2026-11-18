"""Download Gaia TAP results only after successful asynchronous completion.

Anonymous Gaia jobs have a short execution limit. Disjoint source-ID sky tiles
keep joins small; a timed-out tile is bisected, never accepted as partial CSV.
"""
from concurrent.futures import ThreadPoolExecutor
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

ENDPOINT = "https://gea.esac.esa.int/tap-server/tap/async"


def completed_csv(query, limit):
    body = urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                      "MAXREC": limit, "QUERY": query, "PHASE": "RUN"}).encode()
    with urlopen(Request(ENDPOINT, data=body), timeout=90) as response:
        job = response.geturl()
    started = time.monotonic()
    while time.monotonic()-started < 240:
        with urlopen(job+"/phase", timeout=90) as response:
            phase = response.read().decode().strip()
        if phase == "COMPLETED":
            with urlopen(job+"/results", timeout=90) as response:
                root = ET.fromstring(response.read())
            result_id = next(iter(root)).attrib["id"]
            # ESA can advertise a duplicated /results/results/ in xlink:href.
            result = job+"/results/"+result_id
            with urlopen(result, timeout=180) as response:
                raw = response.read()
            if not raw.startswith(b"source_id,") or raw.count(b"\n")-1 >= limit:
                raise ValueError("Invalid or possibly row-limited Gaia result")
            return raw
        if phase in ("ERROR", "ABORTED"):
            with urlopen(job+"/error", timeout=90) as response:
                error = response.read().decode()
            raise RuntimeError(error)
        time.sleep(3)
    raise TimeoutError(job)


def tiled_csv(template, selection, extra, limit, cache, name):
    def tile(bounds):
        lo, hi, depth = bounds
        path = cache / f"gaia-bj-{name}-{lo}-{hi}.csv"
        split_path = path.with_suffix(".split")
        if path.exists():
            return path.read_bytes()
        # OFFSET 0 deliberately isolates the Gaia preselection. Repeat the
        # identical ID bounds on the distance table so the planner can avoid a
        # full 1.47-billion-row distance-table scan for each small sky tile.
        query = template.format(selection=f"({selection}) AND source_id >= {lo} AND source_id < {hi}",
                                distance_selection=f"{extra} AND d.source_id >= {lo} AND d.source_id < {hi}", limit=limit)
        def children():
            mid = (lo+hi)//2
            left, right = tile((lo, mid, depth+1)), tile((mid, hi, depth+1))
            return left+right.split(b"\n", 1)[1]
        if split_path.exists():
            raw = children()
        else:
            try:
                raw = completed_csv(query, limit)
            except RuntimeError as error:
                timeout = "Maximum execution time" in str(error) or "canceling statement due to user request" in str(error)
                if not timeout or depth >= 8:
                    raise
                print(f"{name} tile {lo}..{hi} exceeded execution time; splitting", flush=True)
                # Remember costly parent jobs across interrupted builds.
                split_path.write_text("Archive execution limit; query the two disjoint child tiles instead.\n")
                raw = children()
        path.write_bytes(raw)
        row_count = raw.count(b"\n")-1
        print(f"{name} sky tile {lo}..{hi}: {row_count} rows", flush=True)
        return raw

    # Twelve HEALPix base pixels cover the source-ID domain exactly.
    step = 2**59
    with ThreadPoolExecutor(max_workers=2) as pool:
        parts = list(pool.map(tile, [(index*step, (index+1)*step, 0) for index in range(12)]))
    return parts[0]+b"".join(part.split(b"\n", 1)[1] for part in parts[1:])
