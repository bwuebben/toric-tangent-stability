"""Reproduce the printed rational data in the uniform cycle proof.

Run from any directory with Python's standard library. The existing exact
research module supplies the independently checked trace and rank gates.
"""
from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "paper" / "anc" / "code"))
import research_graph_targets as cycles


def latex(value):
    x = Fraction(value)
    if x.denominator == 1:
        return str(x.numerator)
    return r"\frac{%d}{%d}" % (x.numerator, x.denominator)


def main():
    rows = cycles.cycle_gates()
    table = []
    for row in rows[:3]:
        values = [str(row["j"])] + ["$" + latex(x) + "$" for x in row["positive_sequences"]]
        table.append(" & ".join(values) + r" \\")
    path = ROOT / "paper" / "tables" / "cycle_seeds.tex"
    path.write_text("\n".join(table) + "\n")
    print("Wrote", path)


if __name__ == "__main__":
    main()
