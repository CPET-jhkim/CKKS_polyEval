from poly_util import *
import os
import sys
from math import sqrt, ceil
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from algorithm import cal_polyEval
    from basic_class import Decomp, Poly, NULL_DECOMP
    from print import print_poly
except ImportError:
    from .algorithm import cal_polyEval
    from .basic_class import Decomp, Poly, NULL_DECOMP
    from .print import print_poly


MULTI_THREAD = True
MAX_WORKERS = os.cpu_count() or 1
FUTURES_CHUNK_SIZE = 5000


def _eval_coeff_worker(coeff) -> Decomp:
    """
    개별 coeff에 대한 Decomp 계산.
    반환값은 Decomp 객체.
    """
    _ = Poly(coeff)
    dcmp: Decomp = cal_polyEval(coeff)
    if dcmp is NULL_DECOMP:
        print(f"DCMP ERROR! coeff: {coeff}")
        sys.exit()
    return dcmp


def _iter_chunks(seq, chunk_size: int):
    for i in range(0, len(seq), chunk_size):
        yield seq[i:i + chunk_size]


def eval_degree(
    deg: int,
    multiThread: bool = True,
    max_workers: int | None = None,
    futures_chunk_size: int = FUTURES_CHUNK_SIZE
) -> list[Decomp]:
    """
    deg에 해당하는 모든 다항식에 대해 Decomp 계산.
    반환값은 모든 Decomp 객체 리스트.
    """
    coeffs: list = list(make_all_polys(deg))
    dcmps: list[Decomp] = []

    if len(coeffs) == 0:
        return dcmps

    if not multiThread:
        for coeff in tqdm(coeffs, total=len(coeffs), desc=f"Degree {deg}", ascii=True):
            try:
                dcmp = _eval_coeff_worker(coeff)
                dcmps.append(dcmp)
            except Exception as e:
                print(f"ERROR!, {coeff}")
                print(repr(e))
                sys.exit(1)

        return dcmps

    worker_count = max_workers or (os.cpu_count() or 1)
    worker_count = min(worker_count, len(coeffs))

    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        with tqdm(total=len(coeffs), desc=f"Degree {deg}", ascii=True) as pbar:
            for sub_coeffs in _iter_chunks(coeffs, futures_chunk_size):
                futures = {
                    executor.submit(_eval_coeff_worker, coeff): coeff
                    for coeff in sub_coeffs
                }

                for future in as_completed(futures):
                    coeff = futures[future]

                    try:
                        dcmp: Decomp = future.result()
                        dcmps.append(dcmp)
                    except Exception as e:
                        print(f"ERROR!, {coeff}")
                        print(repr(e))
                        sys.exit(1)

                    pbar.update(1)
        
    return dcmps


def print_decomp(dcmp: Decomp):
    poly = Poly(dcmp.coeff)
    poly.print("type")
    poly.print()
    print(f"DCMP:\t{dcmp.restore_dcmp()}")
    dcmp.comp.print_params()
    print()


if __name__ == '__main__':

    for deg in range(2, 9):
        dcmps: list[Decomp] = eval_degree(
            deg,
            multiThread=MULTI_THREAD,
            max_workers=MAX_WORKERS,
            futures_chunk_size=FUTURES_CHUNK_SIZE
        )

        if len(dcmps) == 0:
            print(f"DEGREE {deg}")
            print("No polynomial generated.")
            print("#" * 20, end="\n\n")
            continue

        depths, cmults, pmults, adds = zip(
            *((d.comp.depth, d.comp.cmult, d.comp.pmult, d.comp.add) for d in dcmps)
        )

        min_depth, max_depth = min(depths), max(depths)
        min_cmult, max_cmult = min(cmults), max(cmults)
        min_pmult, max_pmult = min(pmults), max(pmults)
        min_add, max_add = min(adds), max(adds)
    
        # 출력
        print(f"DEGREE {deg}")
        # est_max_cmult = 2 * ceil(sqrt(deg))
        print(f"Depth: {min_depth} ~ {max_depth}")
        print(f"CMult: {min_cmult} ~ {max_cmult}")
        print(f"PMult: {min_pmult} ~ {max_pmult}")
        print(f"Add: {min_add} ~ {max_add}")

        # dcmp.coeff 내부에서 is_integer()가 False인 계수 개수
        # 이 값이 클수록 정렬상 뒤로 가도록 마지막 기준에 둔다.
        def non_integer_count(dcmp):
            return sum(not c.is_integer() for c in dcmp.coeff)

        # DEPTH 기준: depth -> cmult -> pmult -> add -> non_integer_count
        # CMULT 기준: cmult -> depth -> pmult -> add -> non_integer_count
        # PMULT 기준: pmult -> depth -> cmult -> add -> non_integer_count
        # ADD 기준: add -> depth -> cmult -> pmult -> non_integer_count

        key_depth = lambda d: (
            d.comp.depth,
            d.comp.cmult,
            d.comp.pmult,
            d.comp.add,
            non_integer_count(d),
        )

        key_cmult = lambda d: (
            d.comp.cmult,
            d.comp.depth,
            d.comp.pmult,
            d.comp.add,
            non_integer_count(d),
        )

        key_pmult = lambda d: (
            d.comp.pmult,
            d.comp.depth,
            d.comp.cmult,
            d.comp.add,
            non_integer_count(d),
        )

        key_add = lambda d: (
            d.comp.add,
            d.comp.depth,
            d.comp.cmult,
            d.comp.pmult,
            non_integer_count(d),
        )

        target_dcmps = [
            ("MIN_DEPTH", min(dcmps, key=key_depth)),
            ("MAX_DEPTH", max(reversed(dcmps), key=key_depth)),

            ("MIN_CMULT", min(dcmps, key=key_cmult)),
            ("MAX_CMULT", max(reversed(dcmps), key=key_cmult)),

            ("MIN_PMULT", min(dcmps, key=key_pmult)),
            ("MAX_PMULT", max(reversed(dcmps), key=key_pmult)),

            ("MIN_ADD", min(dcmps, key=key_add)),
            ("MAX_ADD", max(reversed(dcmps), key=key_add)),
        ]

        for label, dcmp in target_dcmps:
            print(f"\n[{label}]")

            poly = Poly(dcmp.coeff)
            poly.print("type")
            poly.print("poly")

            print(f"DCMP:\t{dcmp.restore_dcmp()}")
            dcmp.comp.print_params()