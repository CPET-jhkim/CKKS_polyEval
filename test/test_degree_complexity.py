import argparse
import json
import sys
from collections import Counter
from pathlib import Path

if __package__ and "." in __package__:
    from ..src import NULL_DECOMP, cal_polyEval, make_all_poly_types
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src import NULL_DECOMP, cal_polyEval, make_all_poly_types


def summarize_degree(degree: int, consider_axi: bool=True, assumption1: bool=True) -> dict:
    poly_types = make_all_poly_types(degree)
    if len(poly_types) != 2 * 3 ** degree:
        raise AssertionError(f"Degree {degree} enumeration is incomplete.")
    initial_cache, intermediate_cache = {}, {}
    distribution, examples = Counter(), {}
    print(f"Degree {degree}: evaluating {len(poly_types)} coefficient types...")
    for coeff_type in poly_types:
        result = cal_polyEval(coeff_type, consider_axi=consider_axi, assumption1=assumption1, initial_cache=initial_cache, intermediate_cache=intermediate_cache)
        name = "".join(coeff_type)
        if result is NULL_DECOMP:
            raise RuntimeError(f"OPD failed for coefficient type {name}.")
        cost = (result.comp.depth, result.comp.cmult, result.comp.pmult, result.comp.add)
        distribution[cost] += 1
        examples.setdefault(cost, name)
    costs = list(distribution)
    metrics = ("depth", "cmult", "pmult", "add")
    return {
        "degree": degree,
        "total": len(poly_types),
        "ranges": {name: {"min": min(cost[i] for cost in costs), "max": max(cost[i] for cost in costs)} for i, name in enumerate(metrics)},
        "distribution": [{"depth": cost[0], "cmult": cost[1], "pmult": cost[2], "add": cost[3], "count": distribution[cost], "example": examples[cost]} for cost in sorted(distribution)],
    }


def print_summary(summary: dict) -> None:
    print(f"\nDegree {summary['degree']}, total {summary['total']}")
    for name, values in summary["ranges"].items():
        print(f"{name}: {values['min']} ~ {values['max']}")
    print("\nDepth | CMult | PMult | Add | Count | Example")
    print("-" * 55)
    for row in summary["distribution"]:
        print(f"{row['depth']:>5} | {row['cmult']:>5} | {row['pmult']:>5} | {row['add']:>3} | {row['count']:>5} | {row['example']}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize OPD complexity for every coefficient type of one degree.")
    parser.add_argument("degree", type=int)
    parser.add_argument("--consider-axi", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--assumption1", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--output", type=Path, help="Optional JSON summary path.")
    args = parser.parse_args()
    summary = summarize_degree(args.degree, args.consider_axi, args.assumption1)
    print_summary(summary)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8") as file:
            json.dump(summary, file, ensure_ascii=False, indent=2)
        print(f"Saved summary to {args.output}")


if __name__ == "__main__":
    main()
