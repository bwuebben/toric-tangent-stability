#!/usr/bin/env python3
"""Verify all 6,281 Picard-number-ten eightfold slope certificates exactly.

The input is the supplied Paffenholz/Obro classification slice. Its upstream
completeness and pairwise nonisomorphism are imported, not regenerated here.
NumPy inverses propose integer matrices only; multiplication over the Python
integers certifies each inverse before use. All supporting inequalities,
localization sums, ranks, coefficients, and product/equality checks are exact.
No other verification or classification module is imported.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import gzip
import hashlib
import json
from math import gcd, prod
import os
from pathlib import Path
import sys
import traceback

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

F = Fraction
PACKAGE = Path(__file__).resolve().parents[1]
EXPECTED_INPUT_SHA256 = "fec018c5d465ebe8ec664b1a698710df7be6cdfef5967c50d1646e57237e58ee"
EXPECTED_COUNTS = {"stable": 64, "strictly semistable": 692, "unstable": 5525}
EXPECTED_COMPONENT_TALLY = {"2": 464, "3": 170, "4": 44, "5": 10, "6": 3, "7": 1}

def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def exact_inverse(matrix):
    """Bareiss elimination on [A|I], followed by exact back substitution."""
    size = len(matrix)
    a = [list(row) + [int(i == j) for j in range(size)]
         for i, row in enumerate(matrix)]
    previous = 1
    for k in range(size - 1):
        pivot_row = next((i for i in range(k, size) if a[i][k]), None)
        require(pivot_row is not None, "singular matrix")
        a[k], a[pivot_row] = a[pivot_row], a[k]
        pivot = a[k][k]
        for i in range(k + 1, size):
            multiplier = a[i][k]
            for j in range(k + 1, 2 * size):
                value = pivot * a[i][j] - multiplier * a[k][j]
                require(value % previous == 0, "Bareiss division is not exact")
                a[i][j] = value // previous
            a[i][k] = 0
        previous = pivot
    require(a[-1][size - 1] != 0, "singular matrix")
    inverse = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for j in range(size):
        for i in range(size - 1, -1, -1):
            inverse[i][j] = (Fraction(a[i][size + j]) - sum(
                a[i][k] * inverse[k][j] for k in range(i + 1, size))) / a[i][i]
    require(all(dot(row, col) == int(i == j)
                for i, row in enumerate(matrix)
                for j, col in enumerate(zip(*inverse))), "incorrect exact inverse")
    return inverse


def rank(rows):
    """Rational forward elimination, without row normalization or RREF."""
    if not rows:
        return 0
    a = [list(map(Fraction, row)) for row in rows]
    height, width, result = len(a), len(a[0]), 0
    for j in range(width):
        pivot = next((i for i in range(result, height) if a[i][j]), None)
        if pivot is None:
            continue
        a[result], a[pivot] = a[pivot], a[result]
        for i in range(result + 1, height):
            if a[i][j]:
                factor = a[i][j] / a[result][j]
                for k in range(j + 1, width):
                    a[i][k] -= factor * a[result][k]
                a[i][j] = 0
        result += 1
        if result == height:
            break
    return result


def geometry(rays, record):
    n, m = len(rays[0]), len(rays)
    facets = record["maximal_cones"]
    require(bool(facets), "empty fan")
    require(all(len(f) == n and len(set(f)) == n
                and all(type(i) is int and 0 <= i < m for i in f)
                for f in facets), "malformed maximal cone")
    keys = [tuple(sorted(f)) for f in facets]
    require(len(set(keys)) == len(keys), "duplicate facet")
    matrices = [[rays[i] for i in f] for f in facets]
    try:
        proposals = np.linalg.inv(np.array(matrices, dtype=float))
    except np.linalg.LinAlgError:
        proposals = [None] * len(facets)
    inverses, ridges = [], defaultdict(list)
    for cone_number, (f, matrix, proposal) in enumerate(zip(facets, matrices, proposals)):
        candidate = None
        if proposal is not None and np.all(np.isfinite(proposal)):
            candidate = [[int(round(float(x))) for x in row] for row in proposal]
            columns = list(zip(*candidate))
            if not all(dot(row, col) == int(i == j)
                       for i, row in enumerate(matrix)
                       for j, col in enumerate(columns)):
                candidate = None
        if candidate is None:
            rational = exact_inverse(matrix)
            require(all(x.denominator == 1 for row in rational for x in row),
                    "facet is not unimodular")
            candidate = [[int(x) for x in row] for row in rational]
        # Integral inverse certified against A proves det(A)=+/-1.
        normal = [sum(row) for row in candidate]
        members = set(f)
        require(all(dot(normal, u) == 1 if i in members else dot(normal, u) < 1
                    for i, u in enumerate(rays)), "facet fails strict exact support")
        for omitted in range(n):
            ridge = tuple(sorted(f[:omitted] + f[omitted + 1:]))
            ridges[ridge].append(cone_number)
        inverses.append(candidate)
    require(all(len(owners) == 2 for owners in ridges.values()), "ridge occurs other than twice")
    adjacency = [set() for _ in facets]
    for left, right in ridges.values():
        adjacency[left].add(right)
        adjacency[right].add(left)
    visited, frontier = {0}, [0]
    while frontier:
        v = frontier.pop()
        for neighbor in adjacency[v] - visited:
            visited.add(neighbor)
            frontier.append(neighbor)
    require(len(visited) == len(facets), "facet adjacency is disconnected")
    require(set().union(*map(set, facets)) == set(range(m)), "ray missing from fan")

    # A vector is accepted only after *all* exact weights are nonzero.  Hashing
    # is a deterministic proposal mechanism, with no random arithmetic claims.
    for attempt in range(10000):
        entropy = hashlib.sha256((record["id"] + ":independent-xi:" + str(attempt)).encode()).digest()
        vector = [int(entropy[j]) % 47 - 23 for j in range(n)]
        if vector == record["localization_vector"]:
            continue
        weights = [[dot(col, vector) for col in zip(*inv)] for inv in inverses]
        if all(all(row) for row in weights):
            break
    else:
        raise ValueError("could not propose a generic independent vector")
    degrees, top = [Fraction(0) for _ in rays], Fraction(0)
    for cone, w in zip(facets, weights):
        first_chern, denominator = sum(w), prod(w)
        common = Fraction(first_chern ** (n - 1), denominator)
        top += first_chern * common
        for ray, weight in zip(cone, w):
            degrees[ray] += weight * common
    require(all(x > 0 and x.denominator == 1 for x in degrees), "nonpositive/nonintegral degree")
    require(top == sum(degrees), "top degree fails sum identity")
    require([int(x) for x in degrees] == record["divisor_degrees"], "independent divisor degrees disagree")
    require(top == record["anticanonical_degree"], "independent anticanonical degree disagrees")
    require(all(sum(d * ray[j] for d, ray in zip(degrees, rays)) == 0
                for j in range(n)), "principal-divisor relation fails")
    return degrees, {"cone_count": len(facets), "ridge_count": len(ridges),
                     "independent_localization_vector": vector,
                     "anticanonical_degree": int(top)}


def witness(rays, delta, proposed, strict):
    members = proposed["members"]
    require(bool(members) and len(members) == len(set(members)) and
            all(type(i) is int and 0 <= i < len(rays) for i in members), "malformed flat witness")
    rows = [rays[i] for i in members]
    dimension = rank(rows)
    require(0 < dimension < len(rays[0]), "witness rank is not proper")
    require(dimension == proposed["rank"], "witness rank disagrees")
    closed = [i for i, ray in enumerate(rays) if rank(rows + [ray]) == dimension]
    require(set(closed) == set(members), "witness is not a closed ray flat")
    excess = sum(delta[i] for i in members) - dimension
    require(excess > 0 if strict else excess == 0, "witness excess has wrong sign")
    require(excess == Fraction(proposed["normalized_excess"]), "witness excess disagrees")


def certificate(rays, degrees, record):
    n, m = len(rays[0]), len(rays)
    delta = [n * d / sum(degrees) for d in degrees]
    require(delta == list(map(Fraction, record["normalized_degrees"])), "normalized degrees disagree")
    status = record["status"]
    if status == "unstable":
        require(record["certificate_type"] == "violating_flat", "wrong unstable certificate type")
        witness(rays, delta, record["witness"], True)
        return {"verified_status": status, "basis_count": 0}
    require(status in ("stable", "strictly semistable"), "unknown classification")
    require(record["certificate_type"] == "positive_rational_basis_decomposition", "wrong inside certificate type")
    bases, coefficients = record["bases"], list(map(Fraction, record["coefficients"]))
    require(len(bases) == len(coefficients) > 0, "basis/coefficient cardinality mismatch")
    require(all(c > 0 for c in coefficients) and sum(coefficients) == 1, "decomposition coefficients fail positivity/sum")
    weighted, arcs = [Fraction(0) for _ in rays], [set() for _ in rays]
    for basis, coefficient in zip(bases, coefficients):
        require(len(basis) == len(set(basis)) == n and
                all(type(i) is int and 0 <= i < m for i in basis), "malformed basis")
        inverse = exact_inverse([rays[i] for i in basis])
        for i in basis:
            weighted[i] += coefficient
        for e in set(range(m)) - set(basis):
            coordinates = [dot(rays[e], col) for col in zip(*inverse)]
            require(all(sum(c * rays[b][j] for b, c in zip(basis, coordinates)) == rays[e][j]
                        for j in range(n)), "basis coordinates do not reconstruct ray")
            arcs[e].update(b for b, c in zip(basis, coordinates) if c)
    require(weighted == delta, "positive rational decomposition does not reconstruct point")
    require(arcs == [set(a) for a in record["exchange_arcs"]], "exchange arcs disagree")
    reach = [[i == j or j in arcs[i] for j in range(m)] for i in range(m)]
    for k in range(m):
        for i in range(m):
            if reach[i][k]:
                for j in range(m):
                    reach[i][j] |= reach[k][j]
    components = {tuple(j for j in range(m) if reach[i][j] and reach[j][i]) for i in range(m)}
    require(components == {tuple(sorted(c)) for c in record["strong_components"]}, "strong components disagree")
    require((len(components) == 1) == (status == "stable"), "classification disagrees with exchange graph")
    if status == "strictly semistable":
        witness(rays, delta, record["witness"], False)
        sinks = [set(c) for c in components if all(arcs[i] <= set(c) for i in c)]
        require(bool(sinks), "exchange graph has no sink component")
        for sink in sinks:
            dimension = rank([rays[i] for i in sink])
            require(sum(delta[i] for i in sink) == dimension, "sink component is not tight")
            require(0 < dimension < n and all(
                (rank([rays[k] for k in sink] + [ray]) == dimension) == (i in sink)
                for i, ray in enumerate(rays)), "sink component is not a proper flat")
        require(set(record["witness"]["members"]) in sinks, "equality witness is not a sink component")
    return {"verified_status": status, "basis_count": len(bases),
            "exchange_arc_count": sum(map(len, arcs)), "strong_component_count": len(components)}



def reduction(matrix):
    """Rational elimination of a complete matrix, returning its pivot columns."""
    a = [list(map(F, row)) for row in matrix]
    if not a:
        return [], []
    cursor, pivots = 0, []
    for column in range(len(a[0])):
        active = next((i for i in range(cursor, len(a)) if a[i][column]), None)
        if active is None:
            continue
        a[cursor], a[active] = a[active], a[cursor]
        divisor = a[cursor][column]
        a[cursor] = [x / divisor for x in a[cursor]]
        for i in range(len(a)):
            if i != cursor and a[i][column]:
                multiple = a[i][column]
                a[i] = [x - multiple * y for x, y in zip(a[i], a[cursor])]
        pivots.append(column)
        cursor += 1
        if cursor == len(a):
            break
    return a[:cursor], pivots


def inverse(matrix):
    n = len(matrix)
    a, pivots = reduction([list(row) + [int(i == j) for j in range(n)]
                           for i, row in enumerate(matrix)])
    require(pivots == list(range(n)), "singular square matrix")
    answer = [row[n:] for row in a]
    require(all(sum(matrix[i][k] * answer[k][j] for k in range(n)) == int(i == j)
               for i in range(n) for j in range(n)), "inverse multiplication failed")
    return answer


def reachability(arcs):
    n = len(arcs)
    rows = [(1 << i) | sum(1 << j for j in arcs[i]) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if rows[i] & (1 << k):
                rows[i] |= rows[k]
    return rows


def strong_components(arcs):
    closure = reachability(arcs)
    answer = {tuple(j for j in range(len(arcs))
                    if closure[i] & (1 << j) and closure[j] & (1 << i))
              for i in range(len(arcs))}
    return sorted(answer)


def nullspace_components(rays):
    """Use fundamental nullspace relations of the full 8-by-18 matrix."""
    n, m = len(rays[0]), len(rays)
    red, pivots = reduction(list(zip(*rays)))
    require(len(pivots) == n, "raw rays do not span")
    arcs = [set() for _ in rays]
    relation_count = 0
    for free in set(range(m)) - set(pivots):
        relation = [F(0) for _ in rays]
        relation[free] = 1
        for row, pivot in zip(red, pivots):
            relation[pivot] = -row[free]
        require(all(sum(relation[i] * rays[i][j] for i in range(m)) == 0
                   for j in range(n)), "nullspace relation failed")
        support = {i for i, x in enumerate(relation) if x}
        for i in support:
            arcs[i] |= support - {i}
        relation_count += 1
    groups = strong_components(arcs)
    dimensions = [rank([rays[i] for i in group]) for group in groups]
    require(relation_count == m - n and sum(dimensions) == n, "components fail direct-sum dimension")
    return groups, dimensions



def structure_certificate(raw, record):
    declared = record["structure"]
    require(declared["id"] == raw["id"], "structural record identity mismatch")
    require(record["status"] == "strictly semistable", "structural certificate is not strictly semistable")
    rays = raw["vertices"]
    n, m = len(rays[0]), len(rays)
    require(record["rays"] == rays and n == 8 and m == 18, "raw rays/dimensions mismatch")
    require(record["input_record_sha256"] == sha(json.dumps(raw, sort_keys=True).encode()),
           "raw record hash mismatch")
    groups, dimensions = nullspace_components(rays)
    actual_groups = [tuple(sorted(g)) for g in declared["toric_product_components"]]
    require(groups == actual_groups, "direct-sum component partition mismatch")
    require(dimensions == declared["component_dimensions"] and len(groups) >= 2, "component dimensions/product claim mismatch")

    cones = {tuple(sorted(cone)) for cone in record["maximal_cones"]}
    require(len(cones) == len(record["maximal_cones"]), "duplicate maximal cones")
    cone_basis = list(min(cones))
    basis_inverse = inverse([rays[i] for i in cone_basis])
    require(all(x.denominator == 1 for row in basis_inverse for x in row), "lattice product basis not unimodular")
    basis_inverse = [[int(x) for x in row] for row in basis_inverse]
    positions = [[j for j, ray in enumerate(cone_basis) if ray in group] for group in groups]
    require([len(p) for p in positions] == dimensions, "lattice basis does not split by dimension")
    for group, allowed in zip(groups, positions):
        for i in group:
            coordinates = [sum(rays[i][k] * basis_inverse[k][j] for k in range(n)) for j in range(n)]
            require(all(coordinates[j] == 0 for j in set(range(n)) - set(allowed)),
                   "ray does not lie in its claimed integral direct summand")
            require(all(sum(coordinates[j] * rays[cone_basis[j]][k] for j in range(n)) == rays[i][k]
                       for k in range(n)), "lattice coordinates do not reconstruct ray")
    projections = []
    for group, dimension in zip(groups, dimensions):
        projection = {tuple(i for i in cone if i in group) for cone in cones}
        require(all(len(cone) == dimension for cone in projection), "projected cone dimension mismatch")
        projections.append(projection)
    factor_counts = list(map(len, projections))
    # The map from a full cone to its component intersections is injective;
    # it already lands in the Cartesian product.  Equal cardinality proves
    # surjectivity without enumerating the Cartesian product.
    require(prod(factor_counts) == len(cones), "fan is not the full Cartesian product")
    require(factor_counts == declared["factor_maximal_cone_counts"] and
           declared["full_cartesian_product_fan_verified"] is True, "factor fan count/flag mismatch")

    delta = [F(n * d, sum(record["divisor_degrees"])) for d in record["divisor_degrees"]]
    require(delta == list(map(F, record["normalized_degrees"])), "normalized degrees mismatch")
    arcs = [set() for _ in rays]
    reconstructed = [F(0) for _ in rays]
    weights = list(map(F, record["coefficients"]))
    require(len(weights) == len(record["bases"]) and all(w > 0 for w in weights) and sum(weights) == 1,
           "positive basis coefficients fail")
    for basis, weight in zip(record["bases"], weights):
        require(len(basis) == len(set(basis)) == n, "basis cardinality mismatch")
        inv = inverse([rays[i] for i in basis])
        for i in basis:
            reconstructed[i] += weight
        for e in set(range(m)) - set(basis):
            coordinates = [sum(rays[e][k] * inv[k][j] for k in range(n)) for j in range(n)]
            arcs[e].update(basis[j] for j, value in enumerate(coordinates) if value)
    require(reconstructed == delta, "positive decomposition identity fails")
    require(arcs == list(map(set, record["exchange_arcs"])), "exchange graph mismatch")
    blocks = strong_components(arcs)
    require({tuple(sorted(g)) for g in record["strong_components"]} == set(blocks), "stored SCC mismatch")
    block_order = [tuple(sorted(g)) for g in declared["equality_blocks"]]
    require(set(blocks) == set(block_order) and len(blocks) == len(block_order), "equality block partition mismatch")
    blocks = block_order
    sizes = [sum(delta[i] for i in block) for block in blocks]
    require(all(x.denominator == 1 and x > 0 for x in sizes), "block slope weights are not positive integers")
    require(list(map(int, sizes)) == declared["block_ranks"], "block quotient dimensions mismatch")

    reduced_classes = []
    free = [i for i in range(m) if i not in cone_basis]
    for block in blocks:
        basis_values = [int(i in block) for i in cone_basis]
        principal_character = [sum(basis_inverse[k][j] * basis_values[j] for j in range(n)) for k in range(n)]
        residual = [int(i in block) - sum(rays[i][k] * principal_character[k] for k in range(n)) for i in range(m)]
        require(all(residual[i] == 0 for i in cone_basis), "Picard quotient elimination failed")
        reduced_classes.append([residual[i] for i in free])
    class_rank = rank(reduced_classes)
    require(class_rank == len(blocks) == declared["block_divisor_rank"], "block divisor classes not independent")

    equality, failures = [], []
    cone_masks = [sum(1 << i for i in cone) for cone in cones]
    for selection in range(1, (1 << len(blocks)) - 1):
        members = {i for j, block in enumerate(blocks) if selection & (1 << j) for i in block}
        if any(arcs[i] - members for i in members):
            continue
        indices = sorted(members)
        dimension = rank([rays[i] for i in indices])
        require(sum(delta[i] for i in indices) == dimension and 0 < dimension < n,
               "outgoing-closed union is not a proper tight set")
        require(all((rank([rays[k] for k in indices] + [ray]) == dimension) == (i in members)
                   for i, ray in enumerate(rays)), "tight union is not a flat")
        mask = sum(1 << i for i in indices)
        if any((mask & cone_mask).bit_count() != dimension for cone_mask in cone_masks):
            failures.append(indices)
        equality.append(indices)
    declared_equalities = [tuple(sorted(g)) for g in declared["proper_equalities"]]
    require(len(declared_equalities) == len(set(declared_equalities)), "duplicate declared equality")
    require({tuple(g) for g in equality} == set(declared_equalities), "complete equality list mismatch")
    require({tuple(g) for g in failures} == {tuple(sorted(g)) for g in declared["condition15_failures"]},
           "condition15 failure list mismatch")
    require(not failures, "condition15 fails")
    return {"id": raw["id"], "passed": True,
            "product_components": list(map(list, groups)), "component_dimensions": dimensions,
            "factor_maximal_cone_counts": factor_counts,
            "lattice_basis_ray_indices": cone_basis, "lattice_basis_inverse": basis_inverse,
            "component_basis_positions": positions, "equality_block_count": len(blocks),
            "proper_equality_count": len(equality), "condition15_failure_count": len(failures),
            "block_divisor_rank": class_rank, "free_divisor_indices": free,
            "block_divisor_classes_in_quotient": reduced_classes,
            "block_quotient_dimensions": list(map(int, sizes)),
            "block_actual_span_dimensions": [rank([rays[i] for i in block]) for block in blocks]}


def verify_one(pair):
    item, record = pair
    name = item["id"]
    result = {"id": name, "certificate_sha256": sha(
        json.dumps(record, sort_keys=True, separators=(",", ":")).encode())}
    try:
        require(record["id"] == name, "record identity mismatch")
        require(record["dimension"] == 8 and record["picard_rank"] == 10,
                "dimension/rank metadata mismatch")
        require(record["input_record_sha256"] == sha(json.dumps(item, sort_keys=True).encode()),
                "raw input-record hash mismatch")
        rays = item["vertices"]
        require(rays == record["rays"], "stored rays differ from raw input")
        require(len(rays) == len(set(map(tuple, rays))) == 18 and
                all(len(u) == 8 and all(type(x) is int for x in u) for u in rays),
                "invalid ray dimensions")
        require(all(gcd(*u) == 1 for u in rays), "nonprimitive ray")
        degrees, geometric_result = geometry(rays, record)
        result.update(geometric_result)
        result.update(certificate(rays, degrees, record))
        result["divisor_degree_count"] = len(degrees)
        if record["status"] == "strictly semistable":
            require("structure" in record, "missing structural certificate")
            result["structure"] = structure_certificate(item, record)
        else:
            require("structure" not in record, "unexpected structural certificate")
        result["passed"] = True
    except Exception:
        result.update(passed=False, traceback=traceback.format_exc())
    return result


def load_inputs(raw_path, bundle_path):
    raw_blob, bundle_blob = raw_path.read_bytes(), bundle_path.read_bytes()
    raw_hash = sha(raw_blob)
    require(raw_hash == EXPECTED_INPUT_SHA256, "raw input file hash differs from the supplied classification slice")
    raw = json.loads(raw_blob)
    bundle = json.loads(gzip.decompress(bundle_blob))
    require(bundle["format"] == "toric-eightfold-slope-certificates-v1", "unknown bundle format")
    require(bundle["input_sha256"] == raw_hash, "bundle/input hash mismatch")
    require((bundle["dimension"], bundle["picard_rank"], bundle["ray_count"], bundle["record_count"])
            == (8, 10, 18, 6281), "bundle population metadata mismatch")
    require(bundle["status_counts"] == EXPECTED_COUNTS, "bundle status-count metadata mismatch")
    all_items = [x for x in raw["polytopes"] if len(x["vertices"]) == 18]
    expected = {x["id"]: x for x in all_items}
    require(len(all_items) == len(expected) == 6281, "raw eighteen-ray population mismatch")
    certificates = bundle["records"]
    records = {x["id"]: x for x in certificates}
    require(len(certificates) == len(records) == 6281 and set(records) == set(expected),
            "certificate population has missing, duplicate, or unexpected identifiers")
    summary = bundle["structural_summary"]
    require(summary == {"record_count": 692, "component_count_tally": EXPECTED_COMPONENT_TALLY,
                        "all_toric_products": True, "condition15_failure_count": 0,
                        "all_block_divisors_independent": True}, "structural summary metadata mismatch")
    return [(item, records[item["id"]]) for item in all_items], raw_hash, sha(bundle_blob)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=PACKAGE / "data/smooth_toric_fano_8d_rho_ge_10.json")
    parser.add_argument("--certificates", type=Path,
                        default=PACKAGE / "data/census8_certificates.json.gz")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--report", type=Path, help="write a deterministic detailed JSON report")
    parser.add_argument("--pilot", action="store_true", help="check twelve deterministic input positions only")
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    pairs, raw_hash, bundle_hash = load_inputs(args.input, args.certificates)
    if args.pilot:
        pairs = [pairs[(len(pairs) - 1) * j // 11] for j in range(12)]
    print(f"Verifying {len(pairs)} certificates with {args.workers} worker(s)", flush=True)
    results = []
    if args.workers == 1:
        iterator = map(verify_one, pairs)
        for number, answer in enumerate(iterator, 1):
            results.append(answer)
            if number <= 5 or number % 250 == 0 or not answer["passed"]:
                print(f"{number}/{len(pairs)} {answer['id']} passed={answer['passed']}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            jobs = [pool.submit(verify_one, pair) for pair in pairs]
            for number, future in enumerate(as_completed(jobs), 1):
                answer = future.result()
                results.append(answer)
                if number <= 5 or number % 250 == 0 or not answer["passed"]:
                    print(f"{number}/{len(pairs)} {answer['id']} passed={answer['passed']}", flush=True)
    results.sort(key=lambda row: row["id"])
    passed = [row for row in results if row["passed"]]
    failures = [row["id"] for row in results if not row["passed"]]
    structures = [row["structure"] for row in passed if "structure" in row]
    counts = dict(Counter(row["verified_status"] for row in passed))
    tally = {str(k): v for k, v in sorted(Counter(len(row["product_components"])
                                               for row in structures).items())}
    report = {
        "scope": "twelve-position pilot" if args.pilot else "all supplied eighteen-ray eightfolds",
        "expected_count": len(pairs), "verified_count": len(passed), "failed_ids": failures,
        "complete_population_verified": len(passed) == 6281 and not args.pilot,
        "status_counts": counts,
        "maximal_cones_checked": sum(row["cone_count"] for row in passed),
        "divisor_degrees_checked": sum(row["divisor_degree_count"] for row in passed),
        "linear_bases_checked": sum(row["basis_count"] for row in passed),
        "strictly_semistable_structures_checked": len(structures),
        "component_count_tally": tally,
        "proper_equalities_checked": sum(row["proper_equality_count"] for row in structures),
        "condition15_failure_count": sum(row["condition15_failure_count"] for row in structures),
        "input_sha256": raw_hash, "certificates_sha256": bundle_hash,
        "verifier_sha256": sha(Path(__file__).read_bytes()),
        "records": results,
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "records"}, sort_keys=True), flush=True)
    if failures:
        for answer in results:
            if not answer["passed"]:
                print(answer["traceback"], file=sys.stderr)
        raise SystemExit(1)
    if not args.pilot:
        require(counts == EXPECTED_COUNTS, "verified status counts disagree")
        require(tally == EXPECTED_COMPONENT_TALLY and len(structures) == 692,
                "verified structural population/counts disagree")
        require(report["maximal_cones_checked"] == 2443089 and
                report["divisor_degrees_checked"] == 113058 and
                report["linear_bases_checked"] == 12329 and
                report["proper_equalities_checked"] == 4208,
                "verified certificate coverage disagrees")


if __name__ == "__main__":
    main()
