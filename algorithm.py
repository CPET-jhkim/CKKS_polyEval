try:
    from .basic_class import Complexity, Poly, Decomp, XI, NULL_DECOMP
    from .poly_util import check_without_dcmp, solve_xn_routes, attach
except ImportError as e:
    from basic_class import Complexity, Poly, Decomp, XI, NULL_DECOMP
    from poly_util import check_without_dcmp, solve_xn_routes, attach
    
from math import log2, ceil
import sys

def cal_polyEval(poly: Poly | list[float], made_powers: set[int] = {0, 1}, is_root: bool = True) -> Decomp:
    if type(poly) != Poly:
        poly = Poly(poly)
    max_deg = poly.deg
    ct = poly.coeff_type
    
    # # 0. Base Case(Without Decomposition)
    # if max_deg <= 0:
    #     res = Decomp(poly.coeff, Complexity())
    #     res.made_powers = {0}
    #     return res
    
    # elif max_deg == 1:
    #     comp_res = Complexity()
    #     comp_res.depth = 0 if ct[-1] == "I" else 1
    #     comp_res.pmult = 0 if ct[-1] == "I" else 1
    #     comp_res.add = 1 if ct[0] != "0" else 0
    #     res = Decomp(poly.coeff, comp_res)
    #     res.made_powers = {0, 1}
    #     return res
    
    results: list[Decomp] = []
    
    # 1. Check evaluation availability without decomposition. = All x^i is memoized.
    if check_without_dcmp(made_powers, poly.coeff):
        comp_res = Complexity()
        val = 1 if poly.coeff_type[-1] == "F" else 0
        comp_res.depth = ceil(log2(max(made_powers) + val))
        comp_res.cmult = 0
        comp_res.pmult = poly.coeff_type[1:].count("F")
        comp_res.add = len(poly.coeff_type) - poly.coeff_type.count("0") - 1
        res = Decomp(poly.coeff, comp_res)
        res.made_powers = made_powers
        return res
        
    # 2. Calculate complexity with decomposition.
    # x^i will be selected only with non-zero coefficients'.
    # def find_pi(coeff_type: list[str]):
    #     return [i for i, ct in enumerate(poly.coeff_type) if ct != "0" and i != 0]
    # possible_i = find_pi(poly.coeff_type)
    possible_i = [i for i in range(1, poly.deg+1)]
    # possible_i = [2]
    for xi in possible_i:
        # Case 1: multA=True
        if ct[max_deg] == "F": #and cal_multA(max_deg, xi):
            XIs = solve_xn_routes(True, xi)
            for class_xi in XIs:
                if class_xi.n in made_powers:
                    class_xi.add_count = 0
                res = process_decomposition(poly, class_xi)
                if res is not None:
                    results.append(res)

        # Case 2: multA=False
        XIs = solve_xn_routes(False, xi)
        for class_xi in XIs:
            if class_xi.n in made_powers:
                class_xi.add_count = 0
            res = process_decomposition(poly, class_xi)
            if res is not None:
                results.append(res)

    # Sort and return
    if len(results) == 0:
        return NULL_DECOMP()
    
    results.sort()
    best = results[0]
    best.made_powers = best.merge_mp()

    return best

def process_decomposition(poly: Poly, xi: XI) -> Decomp | None:
    max_deg = poly.deg
    ct = poly.coeff_type
    val = 1 if ct[max_deg] == "F" else 0
    
    # Depth, CMult values are fixed by the degree of polynomial.
    # If the decomposition's complexity exceeds the target, we skip it.
    target_depth = ceil(log2(max_deg + val))
    target_cmult = max_deg - 1
    
    poly_p, poly_q = poly.seperate(xi.n, xi.multA)
    comp_i = Complexity(xi)
    
    # Condition 1: poly_p, poly_q are all empty.
    if poly_p.is_empty() and poly_q.is_empty():
        if comp_i.depth <= target_depth and comp_i.cmult <= target_cmult:
            result = Decomp(poly.coeff, comp_i)
            result.update(xi, None, None)
            return result        
        return None
    
    # Condition 2: poly_q is empty
    if not poly_p.is_empty() and poly_q.is_empty():
        decomp_p: Decomp = cal_polyEval(poly_p, xi.made_powers, is_root=False)
        comp_pi = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
        if comp_pi.depth <= target_depth and comp_pi.cmult <= target_cmult:
            result = Decomp(poly.coeff, comp_pi)
            result.update(xi, decomp_p, None)
            return result
        return None
    
    # Condition 3: poly_p is empty
    if poly_p.is_empty() and not poly_q.is_empty():
        decomp_q: Decomp = cal_polyEval(poly_q, xi.made_powers, is_root=False)
        comp_iq = attach(xi, comp_i, poly_q, decomp_q.comp, '+')
        if comp_iq.depth <= target_depth and comp_iq.cmult <= target_cmult:
            result = Decomp(poly.coeff, comp_iq)
            result.update(xi, None, decomp_q)
            return result
        return None

    # Condition 4: all poly is not empty
    if not poly_p.is_empty() and not poly_q.is_empty():
        if len(poly_p.coeff) <= len(poly_q.coeff):
            decomp_p = cal_polyEval(poly_p, xi.made_powers, is_root=False)
            decomp_q = cal_polyEval(poly_q, decomp_p.made_powers, is_root=False)
        else:
            decomp_q = cal_polyEval(poly_q, xi.made_powers, is_root=False)
            decomp_p = cal_polyEval(poly_p, decomp_q.made_powers, is_root=False)
        
        if decomp_p.is_empty() or decomp_q.is_empty():
            print(f"Failed coeff: {poly.coeff}")
            sys.exit()
        
        comp_pi = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
        comp_piq = attach(None, comp_pi, poly_q, decomp_q.comp, '+')

        if comp_piq.depth == target_depth and comp_piq.cmult <= target_cmult:
            result = Decomp(poly.coeff, comp_piq)
            result.update(xi, decomp_p, decomp_q)
            return result
        return None

    
def process_decomposition_horner(poly: list | Poly, is_root=False) -> Decomp:
    if not isinstance(poly, Poly):
        poly = Poly(poly)
        
    if poly.deg == 0:
        return Decomp(poly.coeff, Complexity(), XI())
    
    if poly.deg % 2 == 1:
        xi = XI(False, 1)
    else:
        xi = XI(False, 2)
        
    xi.add_count = 1 if is_root else 0
    poly_p, poly_q = poly.seperate(xi.n)
    xi.made_powers = {0, 1, 2}
    comp_i = Complexity()
    comp_i.insert_value(xi.depth, xi.add_count, xi.pmult, 0)
    
    dcmp_p = process_decomposition_horner(poly_p)
    dcmp_q = process_decomposition_horner(poly_q)
    
    comp_pi = attach(xi, comp_i, poly_p, dcmp_p.comp, 'x')
    if poly_q.coeff != []:
        comp_piq = attach(None, comp_pi, poly_q, dcmp_q.comp, '+')
    
    res = Decomp(poly.coeff, comp_piq, xi)
    return res