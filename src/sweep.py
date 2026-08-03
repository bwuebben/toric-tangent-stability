#!/usr/bin/env python3
"""Sweep tangent-bundle stability over the full classified database of smooth
toric Fano n-folds (data/smooth_toric_fano_<n>d.json) and cross-tabulate the
verdict against the Kahler--Einstein (barycenter = 0) criterion.

Usage:  python3 sweep.py 3              # all 18 smooth toric Fano 3-folds
        python3 sweep.py 4              # all 124 smooth toric Fano 4-folds
        python3 sweep.py 6 <k> <N>      # shard k of N (strided), writes
                                        #   output/sweep_6d_shard<k>.json
        python3 sweep.py 6 merge <N>    # merge N shards -> sweep_6d.json + summary
"""

import json
import os
import sys
from collections import Counter

from toric_stability import analyze, fmt_frac

HERE = os.path.dirname(os.path.abspath(__file__))


def load_db(dim):
    path = os.path.join(HERE, "..", "data", f"smooth_toric_fano_{dim}d.json")
    with open(path) as f:
        return json.load(f)


def summarize(results, out):
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print("=" * 118)
    tab = Counter((r["verdict"], r["kahler_einstein"]) for r in results)
    verd = Counter(r["verdict"] for r in results)
    print(f"\ntotals ({len(results)} varieties):")
    for v in ("stable", "strictly semistable", "unstable"):
        ke = tab.get((v, True), 0)
        nke = tab.get((v, False), 0)
        print(f"  {v:<22} {verd.get(v, 0):>4}   (KE: {ke}, non-KE: {nke})")
    print(f"\nwrote {os.path.normpath(out)}")


def run(dim, entries, out):
    results = []
    for entry in entries:
        rays = [tuple(v) for v in entry["vertices"]]
        res = analyze(rays, name=str(entry.get("id", "?")))
        res["extra"] = entry.get("extra", {})
        results.append(res)
        bary0 = "KE " if res["kahler_einstein"] else "   "
        print(f"#{res['name']:<12} rays {res['num_rays']:>2}  "
              f"(-K)^n = {res['anticanonical_degree']:>4}  "
              f"mu(T) = {fmt_frac(res['mu_TX']):>7}  "
              f"max mu(V) = {fmt_frac(res['max_subsheaf_slope']):>7}  "
              f"{bary0} {res['verdict']}", flush=True)
    summarize(results, out)


def main(dim, argv):
    db = load_db(dim)
    polys = db["polytopes"]
    outdir = os.path.join(HERE, "..", "output")
    if len(argv) == 2 and argv[0] == "merge":
        nshards = int(argv[1])
        results = []
        for k in range(nshards):
            with open(os.path.join(outdir, f"sweep_{dim}d_shard{k}.json")) as f:
                results.extend(json.load(f))
        assert len(results) == len(polys), (len(results), len(polys))
        results.sort(key=lambda r: r["name"])
        summarize(results, os.path.join(outdir, f"sweep_{dim}d.json"))
        return
    if len(argv) == 2:
        k, nshards = int(argv[0]), int(argv[1])
        entries = polys[k::nshards]          # strided: big cases distribute
        print(f"shard {k}/{nshards}: {len(entries)} of {len(polys)} polytopes")
        run(dim, entries, os.path.join(outdir, f"sweep_{dim}d_shard{k}.json"))
        return
    print(f"loaded {len(polys)} smooth Fano polytopes of dim {dim} "
          f"(source: {db.get('source', '?')})")
    print("=" * 118)
    run(dim, polys, os.path.join(outdir, f"sweep_{dim}d.json"))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3, sys.argv[2:])
