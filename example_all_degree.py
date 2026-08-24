from poly_util import *
import sys
from math import sqrt, ceil
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from algorithm import cal_polyEval
    from basic_class import Decomp, Poly
    from print import print_poly
except:
    from .algorithm import cal_polyEval
    from .basic_class import Decomp, Poly
    from .print import print_poly

if __name__ == '__main__': 

    for deg in range(2, 9):
        coeffs: list = make_all_polys(deg)
        dcmps: list[Decomp] = []
        for coeff in coeffs:
            try:
                poly = Poly(coeff)
                dcmp = cal_polyEval(coeff)
                dcmps.append(dcmp)
            except:
                print(f"ERROR!, {coeff}")
                sys.exit()
            # poly.print("type")
            # poly.print("poly")
            # print(f"DCMP:\t{dcmp.restore_dcmp()}")
            # dcmp.comp.print_params()
            # print()
                
        # max depth/cmult계산
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
        # print(f"Estimated Depth: {ceil(log2(deg))} ~ {ceil(log2(deg+1))}")
        # print(f"Estimated CMult: {ceil(log2(deg))} ~ {ceil(log2(deg+1))}")
        print("#"*20, end="\n\n")
