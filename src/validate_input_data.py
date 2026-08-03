#!/usr/bin/env python3
"""Exact validator for the converted smooth-toric-Fano input files.

For every ray configuration this independently reconstructs the facets of
the fan polytope and checks the defining smooth-Fano conditions: primitive
rays, full dimension, origin strictly inside, every facet a unimodular
simplex at lattice distance one, and every listed ray a genuine vertex.  It
also checks the stored polyDB combinatorial metadata and, where polyDB supplies
``LATTICE_VOLUME``, compares it with the independently stored sweep degree.

Usage:  python3 src/validate_input_data.py          # dimensions 3,4,5,6
        python3 src/validate_input_data.py 3 4      # selected dimensions
"""

import json
import math
import os
import sys

from toric_stability import det, hull_facets, rank


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def validate_dimension(dim):
    data_path = os.path.join(ROOT, "data", f"smooth_toric_fano_{dim}d.json")
    output_path = os.path.join(ROOT, "output", f"sweep_{dim}d.json")
    with open(data_path) as f:
        data = json.load(f)
    with open(output_path) as f:
        degree_by_id = {row["name"]: row["anticanonical_degree"]
                        for row in json.load(f)}

    assert data["dim"] == dim
    assert data["count"] == len(data["polytopes"])
    volume_checks = 0
    for entry in data["polytopes"]:
        name = entry["id"]
        rays = [tuple(v) for v in entry["vertices"]]
        assert rays and all(len(v) == dim for v in rays), name
        assert len(set(rays)) == len(rays), name
        assert rank(rays) == dim, name
        assert all(math.gcd(*[abs(x) for x in ray]) == 1 for ray in rays), name

        facets = hull_facets(rays)
        assert facets, name
        used = set()
        for normal, offset, points in facets:
            # hull_facets orients supporting equations as <normal,x> >= offset.
            assert offset == -1, (name, normal, offset)
            assert len(points) == dim, (name, normal, len(points))
            assert abs(det(points)) == 1, (name, normal, points)
            used.update(points)
        assert used == set(rays), (name, "redundant or missing ray")

        extra = entry["extra"]
        dual = extra["dual_polydb"]
        assert extra["n_vertices"] == len(rays), name
        assert extra["picard_rank"] == len(rays) - dim, name
        assert dual["N_VERTICES"] == len(facets), name
        assert dual["N_FACETS"] == len(rays), name
        if "LATTICE_VOLUME" in dual:
            assert dual["LATTICE_VOLUME"] == degree_by_id[name], name
            volume_checks += 1

    print(f"dim {dim}: {len(data['polytopes'])} smooth reflexive inputs valid; "
          f"independent volume checks {volume_checks}")


def main(argv):
    dims = [int(arg) for arg in argv] if argv else [3, 4, 5, 6]
    assert all(dim in (3, 4, 5, 6) for dim in dims)
    for dim in dims:
        validate_dimension(dim)
    print("INPUT DATA VALIDATION: PASSED")


if __name__ == "__main__":
    main(sys.argv[1:])
