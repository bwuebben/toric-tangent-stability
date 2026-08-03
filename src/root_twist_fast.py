#!/usr/bin/env python3
"""Root-twist family cross-check (dims 4-8), fast + logged + resumable.

Same family X(t_1..t_k) as root_twist.py, but:
  * uses toric_fast.analyze_fast (floating-point candidate generation with
    exact evaluation of retained data),
  * logs timestamped per-case and per-level progress to output/root_twist.log,
  * writes each verdict to output/root_twist.json as soon as it is computed, and
    skips cases already present there on restart (resumable).

The exact prism formulas in prism_check.py are the independent certificate
for the reported high-dimensional values.

Usage:  python3 root_twist_fast.py            # run all cases, resume if partial
        python3 root_twist_fast.py --force     # ignore existing json, redo all
"""

import json
import os
import sys
import time

from toric_fast import analyze_fast
from toric_stability import fmt_frac
from root_twist import root_twist_rays

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "output")
LOG = os.path.join(OUT, "root_twist.log")
JSON = os.path.join(OUT, "root_twist.json")

# (label, twists, expected note).  Ordered easy -> hard so early results land fast.
CASES = [
    ("k=2 (a,-a)              dim4 [zero-sum]",    ["a", "-a"]),
    ("k=3 (a,b,-(a+b))        dim5 [= F.5D.0611]", ["a", "b", "-(a+b)"]),
    ("k=3 (0,0,0)             dim5 [product]",     ["0", "0", "0"]),
    ("k=4 (a,-a,b,-b)         dim6 [zero-sum]",    ["a", "-a", "b", "-b"]),
    ("k=4 (a,b,-(a+b),0)      dim6 [one untwisted]", ["a", "b", "-(a+b)", "0"]),
    ("k=4 (a,b,-(a+b),a)      dim6 [sum!=0]",      ["a", "b", "-(a+b)", "a"]),
    ("k=5 (a,b,-(a+b),a,-a)   dim7 [zero-sum]",    ["a", "b", "-(a+b)", "a", "-a"]),
    ("k=6 (a,b,-(a+b))x2      dim8 [zero-sum]",    ["a", "b", "-(a+b)", "a", "b", "-(a+b)"]),
]


def logline(msg):
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{stamp}] {msg}"
    with open(LOG, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def load_done(force):
    if force or not os.path.exists(JSON):
        return {}
    with open(JSON) as f:
        rows = json.load(f)
    for row in rows:
        row.pop("seconds", None)
    return {r["label"]: r for r in rows}


def save(results):
    with open(JSON, "w") as f:
        json.dump(list(results.values()), f, indent=1)


def main(argv):
    force = "--force" in argv
    results = load_done(force)
    logline(f"=== root_twist_fast start (have {len(results)} cached cases) ===")
    for label, twists in CASES:
        if label in results:
            logline(f"skip (cached): {label} -> {results[label]['verdict']}")
            continue
        rays = root_twist_rays(twists)
        t0 = time.time()
        logline(f"BEGIN {label}  ({len(rays)} rays, dim {2 + len(twists)})")
        try:
            res = analyze_fast(rays, label, log=logline)
        except AssertionError as e:
            logline(f"  NOT smooth/reflexive: {e}")
            res = {"verdict": "not-smooth", "error": str(e)}
        dt = time.time() - t0
        res["label"] = label
        results[label] = res
        save(results)                      # persist immediately
        v = res.get("verdict", "?")
        ke = res.get("kahler_einstein")
        deg = res.get("anticanonical_degree", "?")
        mv = res.get("max_subsheaf_slope")
        logline(f"DONE  {label}: {v}  KE={ke}  (-K)^n={deg}  "
                f"max mu(V)={fmt_frac(mv) if mv else '?'}  [{dt:.1f}s]")
    save(results)
    logline("=== root_twist_fast complete ===")
    # human-readable summary
    for label, _ in CASES:
        r = results.get(label, {})
        logline(f"SUMMARY {label}: {r.get('verdict','?')}  "
                f"KE={r.get('kahler_einstein')}")


if __name__ == "__main__":
    main(sys.argv[1:])
