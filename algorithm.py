from complexity import Complexity, attach
from polynomial import Poly, Decomp
from math import log2, ceil
from util import *

def calculate(poly: Poly, made_powers: set[int], is_root: bool = True):
    max_deg = poly.deg
    ct = poly.coeff_type
    
    # 0. 분해가 필요없는 기초다항식 처리
    if max_deg <= 0:
        return 0, Decomp(poly.coeff, Complexity(), made_powers)
    
    elif max_deg == 1:
        comp_res = Complexity()
        comp_res.depth = 0 if ct[-1] == "I" else 1
        comp_res.pmult = 0 if ct[-1] == "I" else 1
        comp_res.add = 1 if ct[0] != "0" else 0
        res = Decomp(poly.coeff, comp_res, made_powers)
        return 1, res

    # 다항식 분해 후 복잡도 연산.
    results = []
    
    # Case 1: multA=True
    if ct[max_deg] == "F" and is_multA_required(max_deg):
        results.extend(process_decomposition(poly, 1, True, made_powers, (1, 0)))

    # Case 2: multA=False
    routes = solve_xn_operation(False, max_deg)
    for mp, ops_list in routes:
        for xi in mp:
            if xi <= max_deg:
                cost = ops_list[xi][1] 
                results.extend(process_decomposition(poly, xi, False, mp, cost))

    # 정렬 후 반환
    results.sort(key=lambda x: x[1])
    if is_root: # debug
        for i, dcmp in results:
            print(f"{'dcmp:':<8}{dcmp.restore_dcmp()}")
            dcmp.comp.print_params()
            print(f"{'multA:':<8}{dcmp.multA}")
            print(f"{'depth:':<8}{dcmp.check_depth()}")
            
    return results[0]

def process_decomposition(poly: Poly, xi: int, multA: bool, made_powers: set[int], comp: tuple[int, int]):
    max_deg = poly.deg
    ct = poly.coeff_type
    results = []
    val = 1 if ct[max_deg] == "F" else 0
    target_depth = ceil(log2(max_deg + val))
    
    poly_p, poly_q = poly.seperate(xi, multA)

    next_powers = made_powers if multA else made_powers
    
    _, decomp_p = calculate(poly_p, next_powers, is_root=False)
    _, decomp_q = calculate(poly_q, next_powers, is_root=False)

    comp_i = Complexity()
    pmult = 1 if val and multA else 0
    comp_i.insert_value(comp[0], comp[1], pmult, 0)

    comp_pi = attach(comp_i, decomp_p.comp, 'x')

    if decomp_q.coeff != []:
        comp_piq = attach(comp_pi, decomp_q.comp, '+')
    else:
        comp_piq = comp_pi
    
    if comp_piq.depth == target_depth:
        dcmp = Decomp(poly.coeff, comp_piq, decomp_p.mp | decomp_q.mp)
        dcmp.update(multA, xi, decomp_p, decomp_q)
        results.append((xi, dcmp))

    return results