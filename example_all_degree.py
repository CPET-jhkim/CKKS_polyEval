from poly_util import *
import sys
from math import sqrt, ceil
try:
    from algorithm import cal_polyEval
    from basic_class import Decomp, Poly
    from print import print_poly
except:
    from .algorithm import cal_polyEval
    from .basic_class import Decomp, Poly
    from .print import print_poly

if __name__ == '__main__': 

    deg = 7
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
        poly.print("type")
        poly.print("poly")
        print(f"DCMP:\t{dcmp.restore_dcmp()}\n")
            
    # max depth/cmult계산
    max_depth = max([d.comp.depth for d in dcmps])
    max_cmult = max([d.comp.cmult for d in dcmps])
    
    # 출력
    print(f"DEGREE {deg}")
    est_max_cmult = 2 * ceil(sqrt(deg))
    print(f"Depth: Max {max_depth}\tCMult: Max {max_cmult} <= {est_max_cmult}")
    print("#"*20, end="\n\n")
