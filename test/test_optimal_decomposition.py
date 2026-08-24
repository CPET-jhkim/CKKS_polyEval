import argparse
import sys
from pathlib import Path

if __package__ and "." in __package__:
    from ..src import NULL_DECOMP, Poly, cal_polyEval
    from ..src.poly_util import print_dcmp_detail
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src import NULL_DECOMP, Poly, cal_polyEval
    from src.poly_util import print_dcmp_detail


def parse_coefficients(raw: str) -> list[float]:
    values = [float(value.strip()) for value in raw.split(",") if value.strip()]
    if not values:
        raise ValueError("At least one coefficient is required.")
    return values


def parse_coefficient_types(raw: str) -> list[str]:
    values = [value.strip().upper() for value in raw.split(",")] if "," in raw else list(raw.strip().upper())
    if not values or any(value not in {"0", "I", "F"} for value in values):
        raise ValueError("Coefficient types must contain only 0, I, and F.")
    return values


def calculate_optimal(poly_input: list[float] | list[str], consider_axi: bool=True, assumption1: bool=True):
    result = cal_polyEval(Poly(poly_input), consider_axi=consider_axi, assumption1=assumption1)
    if result is NULL_DECOMP:
        raise RuntimeError("No decomposition was found.")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate the optimal decomposition for coefficients or coefficient types.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--coeff", help="Comma-separated coefficients ordered from the constant term, e.g. 1,2.5,0,4")
    group.add_argument("--type", help="Coefficient types, e.g. I,F,0,I or IF0I")
    parser.add_argument("--consider-axi", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--assumption1", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    poly_input = parse_coefficients(args.coeff) if args.coeff is not None else parse_coefficient_types(args.type)
    result = calculate_optimal(poly_input, args.consider_axi, args.assumption1)
    print(f"Input: {poly_input}")
    print(f"Optimal decomposition: {result.restore_dcmp()}")
    print(f"Complexity (Depth|CMult|PMult|Add): {result.comp.return_params()}")
    print_dcmp_detail(result, 0)


if __name__ == "__main__":
    main()
