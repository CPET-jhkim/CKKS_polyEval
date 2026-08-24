import argparse
import sys
from pathlib import Path

if __package__ and "." in __package__:
    from ..src import NULL_DECOMP, cal_polyEval, make_all_poly_types
    from ..src.save import get_poly_type_key, load_cache, save_cache, serialize_decomp
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src import NULL_DECOMP, cal_polyEval, make_all_poly_types
    from src.save import get_poly_type_key, load_cache, save_cache, serialize_decomp


PROJECT_DIR = Path(__file__).resolve().parents[1]


def save_all_decompositions(output: Path, min_degree: int=2, max_degree: int=8, consider_axi: bool=True, assumption1: bool=True, resume: bool=True) -> dict:
    if min_degree < 1 or max_degree < min_degree:
        raise ValueError("Use 1 <= min_degree <= max_degree.")
    cache = load_cache(str(output)) if resume else {}
    initial_cache, intermediate_cache = {}, {}
    for degree in range(min_degree, max_degree + 1):
        poly_types = make_all_poly_types(degree)
        if len(poly_types) != 2 * 3 ** degree:
            raise AssertionError(f"Degree {degree} enumeration is incomplete.")
        completed = 0
        print(f"Degree {degree}: evaluating {len(poly_types)} coefficient types...")
        for coeff_type in poly_types:
            key = get_poly_type_key(coeff_type, assumption1)
            if key not in cache:
                result = cal_polyEval(coeff_type, consider_axi=consider_axi, assumption1=assumption1, initial_cache=initial_cache, intermediate_cache=intermediate_cache)
                if result is NULL_DECOMP:
                    raise RuntimeError(f"OPD failed for coefficient type {key}.")
                cache[key] = serialize_decomp(result)
            completed += 1
        save_cache(cache, str(output))
        print(f"Degree {degree}: saved {completed}/{len(poly_types)} cases to {output}")
    return cache


def main() -> None:
    parser = argparse.ArgumentParser(description="Save optimal decompositions for every degree 2-8 coefficient type as JSON.")
    parser.add_argument("--min-degree", type=int, default=2)
    parser.add_argument("--max-degree", type=int, default=8)
    parser.add_argument("--consider-axi", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--assumption1", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or PROJECT_DIR / "data" / f"decomp_cache_v3_axi{int(args.consider_axi)}_A{int(args.assumption1)}.json"
    cache = save_all_decompositions(output, args.min_degree, args.max_degree, args.consider_axi, args.assumption1, args.resume)
    print(f"Saved {len(cache)} total coefficient types. The JSON root contains coefficient-type keys directly.")


if __name__ == "__main__":
    main()
