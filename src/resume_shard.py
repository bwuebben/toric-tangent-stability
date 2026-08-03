#!/usr/bin/env python3
"""Crash-recovery helper for the dim-6 (and any dim) sharded sweep.

WHY THIS EXISTS
---------------
`sweep.py 6 <k> <N>` writes its shard JSON (output/sweep_6d_shard<k>.json)
*only once, at the very end*, after all ~1270 strided varieties of the shard
are analyzed.  There is no incremental checkpoint.  So if the machine dies
mid-shard:
  * output/sweep_6d_shard<k>.json does NOT exist yet, and
  * a naive relaunch of `sweep.py 6 <k> <N>` restarts that shard from
    variety 0, re-paying for every heavy variety already computed.

The per-variety verdict LINES do survive, in output/sweep_6d_shard<k>.log
(one `#<name> ... <verdict>` line per finished variety).  This helper reads
those lines to learn which names a shard already finished, and re-runs ONLY
the remaining ones -- writing them to output/sweep_6d_shard<k>.partial.json.

USAGE
-----
  # 1) after a crash, restart each shard's REMAINING work (6 shards, N=6):
  python3 src/resume_shard.py 6 0 6
  python3 src/resume_shard.py 6 1 6
  ...   (through shard 5)  -- run these in parallel, same as the original sweep

  # 2) once every shard is complete (each shard's .log covers all its strided
  #    entries), assemble the final full-detail JSON:
  python3 src/resume_shard.py 6 assemble 6
  #    -> output/sweep_6d.json   (+ prints the verdict tally)

RECOVERY-DETAIL CAVEAT
----------------------
A variety finished BEFORE the crash has full detail (witness subspace,
barycenter, facet degrees) only if it landed in a completed shard JSON or a
.partial.json produced by this helper.  For any variety that exists only as a
pre-crash .log line, `assemble` emits a verdict-only stub
(`"detail": "log-only"`).  That is enough for every COUNT in the paper's
tables (verdict tallies, KE cross-tab), but NOT for the ancillary
witness/degree files.  If you need the definitive full-detail ancillary
files, do one clean uninterrupted `sweep.py` run at the end; the counts from
logs+partials are already trustworthy and can be checked against it.
"""

import json
import os
import re
import sys
from collections import Counter

from toric_stability import analyze  # same dir on sys.path when run from src/
from sweep import load_db, summarize

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "output")

# a finished log line looks like:
#   #F.6D.0803    rays 12  (-K)^n = 25569  ...  KE  stable
#   #F.6D.0000    rays  9  (-K)^n = 247667 ...      unstable
# order matters: match the more-specific suffix first, since
# "unstable".endswith("stable") is True and "strictly semistable" ends "stable"
_VERDICTS = ("strictly semistable", "unstable", "stable")


def done_names_from_log(dim, k):
    """Set of variety names whose analysis completed (one line in the log)."""
    log = os.path.join(OUT, f"sweep_{dim}d_shard{k}.log")
    names = set()
    if not os.path.exists(log):
        return names
    with open(log) as f:
        for line in f:
            if line.startswith("#"):
                names.add(line[1:].split()[0])
    return names


def parse_log_stub(line):
    """Verdict-only record recovered from a single completed log line."""
    name = line[1:].split()[0]
    verdict = next((v for v in _VERDICTS if line.rstrip().endswith(v)), "?")
    ke = " KE " in line
    return {"name": name, "verdict": verdict, "kahler_einstein": ke,
            "detail": "log-only"}


def resume(dim, k, nshards):
    polys = load_db(dim)["polytopes"]
    entries = polys[k::nshards]
    done = done_names_from_log(dim, k)
    remaining = [e for e in entries if str(e.get("id", "?")) not in done]
    print(f"shard {k}/{nshards}: {len(entries)} strided, {len(done)} already "
          f"in log, {len(remaining)} to (re)compute")
    results = []
    for e in remaining:
        rays = [tuple(v) for v in e["vertices"]]
        res = analyze(rays, name=str(e.get("id", "?")))
        res["extra"] = e.get("extra", {})
        results.append(res)
        print(f"#{res['name']:<12} {res['verdict']}", flush=True)
    out = os.path.join(OUT, f"sweep_{dim}d_shard{k}.partial.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {len(results)} newly-computed varieties -> "
          f"{os.path.normpath(out)}")


def assemble(dim, nshards):
    """Combine complete shard JSONs + .partial.json + log-only stubs into the
    final sweep_<dim>d.json, deduplicating by name (full detail wins)."""
    polys = load_db(dim)["polytopes"]
    by_name = {}

    def absorb(path, prefer):
        if not os.path.exists(path):
            return 0
        with open(path) as f:
            recs = json.load(f)
        n = 0
        for r in recs:
            nm = r["name"]
            # full-detail record always wins over a log-only stub
            if nm not in by_name or (prefer == "full"
                                     and by_name[nm].get("detail") == "log-only"):
                by_name[nm] = r
                n += 1
        return n

    for k in range(nshards):
        absorb(os.path.join(OUT, f"sweep_{dim}d_shard{k}.json"), "full")
        absorb(os.path.join(OUT, f"sweep_{dim}d_shard{k}.partial.json"), "full")
    # fill any gaps from log stubs
    stubs = 0
    for k in range(nshards):
        log = os.path.join(OUT, f"sweep_{dim}d_shard{k}.log")
        if not os.path.exists(log):
            continue
        with open(log) as f:
            for line in f:
                if line.startswith("#"):
                    nm = line[1:].split()[0]
                    if nm not in by_name:
                        by_name[nm] = parse_log_stub(line)
                        stubs += 1

    results = sorted(by_name.values(), key=lambda r: r["name"])
    log_only = sum(1 for r in results if r.get("detail") == "log-only")
    print(f"assembled {len(results)} varieties "
          f"({log_only} verdict-only from logs, {stubs} stubs added)")
    missing = len(polys) - len(results)
    if missing:
        print(f"WARNING: {missing} of {len(polys)} varieties still missing "
              f"(some shard neither finished nor logged them)")
    summarize(results, os.path.join(OUT, f"sweep_{dim}d.json"))


def main(argv):
    dim = int(argv[0])
    if argv[1] == "assemble":
        assemble(dim, int(argv[2]))
    else:
        resume(dim, int(argv[1]), int(argv[2]))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
