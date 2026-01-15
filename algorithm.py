from basic_class import Poly, Decomp
from math import log2, ceil
from util import *

def calculate(poly: Poly, made_powers: set[int], is_root: bool = True) -> Decomp | bool:
    max_deg = poly.deg
    ct = poly.coeff_type
    
    # 0. 분해가 필요없는 기초다항식 처리
    if max_deg <= 0:
        res = Decomp(poly.coeff, Complexity())
        res.made_powers = {0}
        return res
    
    elif max_deg == 1:
        comp_res = Complexity()
        comp_res.depth = 0 if ct[-1] == "I" else 1
        comp_res.pmult = 0 if ct[-1] == "I" else 1
        comp_res.add = 1 if ct[0] != "0" else 0
        res = Decomp(poly.coeff, comp_res)
        res.made_powers = {0, 1}
        return res
    
    results: list[Decomp] = []
    
    # 1. 분해 없이 연산 가능한지 확인
    if check_without_dcmp(made_powers, poly.coeff):
        comp_res = Complexity()
        val = 1 if poly.coeff_type[-1] == "F" else 0
        comp_res.depth = ceil(log2(max(made_powers) + val))
        comp_res.cmult = 0
        comp_res.pmult = poly.coeff_type[1:].count("F")
        comp_res.add = len(poly.coeff_type) - poly.coeff_type.count("0") - 1
        res = Decomp(poly.coeff, comp_res)
        res.made_powers = made_powers
        results.append(res)
        
    # 2. 다항식 분해 후 복잡도 연산.
    # i는 계수 중 0이 아닌 것들로만 구성.
    def find_pi(coeff_type: list[str]):
        return [i for i, ct in enumerate(poly.coeff_type) if ct != "0" and i != 0]
    possible_i = find_pi(poly.coeff_type)
    for xi in possible_i:
        # Case 1: multA=True
        if ct[max_deg] == "F": #and cal_multA(max_deg, xi):
            XIs = solve_xn_routes(True, xi)
            for class_xi in XIs:
                if class_xi.n in made_powers:
                    class_xi.add_count = 0
                res = process_decomposition(poly, class_xi)
                if type(res) == list:
                    results.extend(res)

        # Case 2: multA=False
        XIs = solve_xn_routes(False, xi)
        for class_xi in XIs:
            if class_xi.n in made_powers:
                class_xi.add_count = 0
            res = process_decomposition(poly, class_xi)
            if type(res) == list:
                results.extend(res)

    # 정렬 후 반환
    if len(results) == 0:
        return False
    
    results.sort()
    if is_root: # debug
        for dcmp in results[:1]:
            print(f"{'dcmp:':<8}{dcmp.restore_dcmp()}")
            dcmp.comp.print_params()
            # print(f"{'multA:':<8}{dcmp.xi.multA}")
            # print(f"{'depth:':<8}{dcmp.check_depth()}")
            print("-"*10)
    
    results[0].made_powers = results[0].merge_mp()
    return results[0]

def process_decomposition(poly: Poly, xi: XI) -> list[Decomp] | bool:
    max_deg = poly.deg
    ct = poly.coeff_type
    results = []
    val = 1 if ct[max_deg] == "F" else 0
    target_depth = ceil(log2(max_deg + val))
    target_cmult = max_deg - 1
    
    poly_p, poly_q = poly.seperate(xi.n, xi.multA)
    
    decomp_p = calculate(poly_p, xi.made_powers, is_root=False)
    if type(decomp_p) == Decomp: 
        decomp_q = calculate(poly_q, decomp_p.made_powers, is_root=False)

        if decomp_q is False:
            return False
        
        comp_i = Complexity()
        comp_i.insert_value(xi.depth, xi.add_count, xi.pmult, 0)

        comp_pi = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
        comp_piq = attach(None, comp_pi, poly_q, decomp_q.comp, '+')

        if comp_piq.depth == target_depth and comp_piq.cmult <= target_cmult:
            dcmp = Decomp(poly.coeff, comp_piq)
            dcmp.update(xi, decomp_p, decomp_q)
            results.append(dcmp)

        return results