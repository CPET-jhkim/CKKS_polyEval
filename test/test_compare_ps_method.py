import argparse
import sys
from pathlib import Path

if __package__ and "." in __package__:
    from ..src import NULL_DECOMP, Poly, cal_polyEval, cal_PSMethod, make_all_poly_types
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src import NULL_DECOMP, Poly, cal_polyEval, cal_PSMethod, make_all_poly_types


def compare_ps_method(min_degree: int=2, max_degree: int=8, consider_axi: bool=True, assumption1: bool=True) -> list[dict]:
    if min_degree < 1 or max_degree < min_degree:
        raise ValueError("Use 1 <= min_degree <= max_degree.")
    initial_cache, intermediate_cache, rows = {}, {}, []
    for degree in range(min_degree, max_degree + 1):
        poly_types = make_all_poly_types(degree)
        if len(poly_types) != 2 * 3 ** degree:
            raise AssertionError(f"Degree {degree} enumeration is incomplete.")
        counts = {"total": len(poly_types), "improved": 0, "depth": 0, "cmult": 0, "pmult": 0, "multiple": 0}
        print(f"Degree {degree}: comparing {len(poly_types)} coefficient types...")
        for coeff_type in poly_types:
            poly = Poly(coeff_type)
            baseline = cal_PSMethod(poly)
            result = cal_polyEval(poly, consider_axi=consider_axi, assumption1=assumption1, initial_cache=initial_cache, intermediate_cache=intermediate_cache)
            name = "".join(coeff_type)
            if baseline is NULL_DECOMP or result is NULL_DECOMP:
                raise RuntimeError(f"Evaluation failed for degree {degree}, type {name}.")
            if result.comp > baseline.comp:
                raise AssertionError(f"OPD is worse than PS for {name}: {result.comp.return_params()} > {baseline.comp.return_params()}")
            reductions = (baseline.comp.depth > result.comp.depth, baseline.comp.cmult > result.comp.cmult, baseline.comp.pmult > result.comp.pmult)
            counts["improved"] += any(reductions)
            counts["depth"] += reductions[0]
            counts["cmult"] += reductions[1]
            counts["pmult"] += reductions[2]
            counts["multiple"] += sum(reductions) >= 2
        rows.append({"degree": degree, **counts})
    return rows


def print_summary(rows: list[dict]) -> None:
    headers = ("degree", "total", "improved", "depth", "cmult", "pmult", "multiple")
    print("\n" + " | ".join(headers))
    print("-" * 70)
    for row in rows:
        print(" | ".join(str(row[name]) for name in headers))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare OPD with the PS-Method for every coefficient type.")
    parser.add_argument("--min-degree", type=int, default=2)
    parser.add_argument("--max-degree", type=int, default=8)
    parser.add_argument("--consider-axi", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--assumption1", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    print_summary(compare_ps_method(args.min_degree, args.max_degree, args.consider_axi, args.assumption1))


if __name__ == "__main__":
    main()
