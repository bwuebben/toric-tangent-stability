"""Exact finite star comparisons m=1,...,14 from the printed volume formula.

The infinite tail is proved analytically in the manuscript. Standalone,
standard-library-only reproduction of the finite table, not a fan scan.
"""
import argparse
from fractions import Fraction
import json
from math import comb, factorial
from pathlib import Path


def volume_and_facet(m):
    t = m + 1
    volume = Fraction(0)
    facet = Fraction(0)
    for j in range(m + 1):
        for h in range(m - j + 1):
            degree = 2 * m - h
            coefficient = (-1) ** j * comb(m, j) * comb(m - j, h)
            volume += coefficient * Fraction((t - j) ** degree, factorial(degree))
            facet += coefficient * Fraction((t - j) ** (degree - 1), factorial(degree - 1))
    return volume, facet


def records():
    rows = []
    for m in range(1, 15):
        volume, facet = volume_and_facet(m)
        assert volume > 0 and facet > 0
        comparison = 1 - facet / volume
        assert (comparison > 0) == (m <= 5)
        rows.append({"arms": m, "volume": str(volume), "facet_volume": str(facet),
                     "one_minus_central_degree": str(comparison)})
    return rows


def table(rows):
    return "\n".join(
        f"{row['arms']} & ${row['one_minus_central_degree']}$ " + r"\\"
        for row in rows
    ) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-table", type=Path)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--check-json", type=Path)
    args = parser.parse_args()
    rows = records()
    result = {"scope": "Exact finite star comparisons m=1,...,14; analytic tail is in manuscript", "rows": rows}
    if args.check_json:
        assert json.loads(args.check_json.read_text()) == result
        print("PASS: all 14 exact comparisons match the stored record.")
    if args.output_table:
        args.output_table.parent.mkdir(parents=True, exist_ok=True)
        args.output_table.write_text(table(rows))
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(result, indent=2) + "\n")
    if not any((args.output_table, args.output_json, args.check_json)):
        print(table(rows), end="")


if __name__ == "__main__":
    main()
