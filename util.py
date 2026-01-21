# util.py

import random
from itertools import product
from math import log2, ceil
from basic_class import Complexity, Poly, XI

def make_type_poly(poly_type: list[str]) -> list[float]:
    '''
    Create random polynomial by input coefficient type
    '''
    assert all(x in {"F", "I", "0"} for x in poly_type)
    
    poly = []
    for ctype in poly_type:
        if ctype == "I":
            poly.append(random.randint(1, 9))
        elif ctype == "F":
            poly.append(float(random.randint(1, 9)) + random.randint(1, 9) / 10.0)
        elif ctype == "0":
            poly.append(0)
    return poly

def make_all_polys(max_deg: int) -> list[list[float]]:
    '''
    Create every possible polynomial of degree n.
    Coefficient type can be 0/I/F.(Max degree can't be 0.)
    '''
    
    gen_zero = lambda: 0.0
    gen_int  = lambda: float(random.randint(1, 10))
    
    gen_float = lambda: float(random.randint(1, 9)) + random.randint(1, 9) / 10.0

    generators_middle = [gen_zero, gen_int, gen_float]
    generators_highest = [gen_int, gen_float]

    iterables = [generators_middle] * max_deg + [generators_highest]
    
    generator_combinations = product(*iterables)

    all_polys = []
    for combo in generator_combinations:
        poly_values = [func() for func in combo]
        all_polys.append(poly_values)
        
    return all_polys
        
def is_multA_required(n: int) -> bool:
    """
    Check whether coefficient can be attached in x^n.
    Depth shouldn't be increased.
    """
    if n == 0:
        return False
        
    return (n & (n - 1)) != 0

def cal_multA(max_deg: int, n: int) -> bool:
    if max_deg == n:
        return False
    
    xi_diff = ceil(log2(n)) < ceil(log2(n + 1))
    p_diff = ceil(log2(max_deg - n)) < ceil(log2(max_deg - n + 1))
    if xi_diff is False and p_diff is True:
        return True

    return False

def solve_xn_routes(multA, n) -> list[XI]:
    """
    Find optimal route to make x^n.
    """

    # ---------------------------------------------------------
    # 1. DP: Compute optimal Depth (parallel time) for all n
    # ---------------------------------------------------------
    # min_depths[i] = (pure_depth, coeff_depth)
    min_depths = {1: (0, 0)} 
    
    for i in range(2, n + 1):
        d_pure, d_coeff = float('inf'), float('inf')
        for j in range(1, i // 2 + 1):
            k = i - j
            
            # Pure (x^j * x^k = x^i)
            # max time when both operands are ready + 1
            p_val = max(min_depths[j][0], min_depths[k][0]) + 1
            if p_val < d_pure: 
                d_pure = p_val
            
            # Coeff (when multA=True, generate ax^i)
            if multA:
                # Case 1: ax^j * x^k -> ax^i
                c_val1 = max(min_depths[j][1], min_depths[k][0]) + 1
                if c_val1 < d_coeff: 
                    d_coeff = c_val1
                
                # Case 2: x^j * ax^k -> ax^i
                c_val2 = max(min_depths[k][1], min_depths[j][0]) + 1
                if c_val2 < d_coeff: 
                    d_coeff = c_val2
                
        min_depths[i] = (d_pure, d_coeff)

    # Set target optimal Depth
    target_opt_depth = min_depths[n][1] if multA else min_depths[n][0]

    # ---------------------------------------------------------
    # 2. Recursive Search: Backtracking and route construction
    # ---------------------------------------------------------
    # Return value: list of dict -> [{'ops': set((a,b)...), 'depth': int}, ...]
    # ops is an unordered set of operations; later sorted to build the route
    
    memo = {}

    def find_paths(target, has_coeff):
        state_key = (target, has_coeff)
        if state_key in memo:
            return memo[state_key]

        # Base Case
        if target == 1:
            return [{'ops': set(), 'depth': 0}]
            
        res = []
        for j in range(1, target // 2 + 1):
            k = target - j
            
            # Case 1: Pure generation (x^target)
            # Only performed when has_coeff is False
            if not has_coeff:
                # Check critical path condition (pruning)
                if max(min_depths[j][0], min_depths[k][0]) + 1 <= min_depths[target][0]:
                    l_list = find_paths(j, False)
                    r_list = find_paths(k, False)
                    
                    for l in l_list:
                        for r in r_list:
                            nd = max(l['depth'], r['depth']) + 1
                            # Collect only paths that exactly match optimal Depth
                            if nd == min_depths[target][0]:
                                new_ops = l['ops'] | r['ops']
                                # Add current operation (tuple ordered as (smaller, larger))
                                new_ops.add(tuple(sorted((j, k))))
                                res.append({
                                    'ops': new_ops,
                                    'depth': nd
                                })
            
            # Case 2: Coeff generation (ax^target)
            # Performed when has_coeff is True
            else:
                # 2-1: ax^j * x^k
                if max(min_depths[j][1], min_depths[k][0]) + 1 <= target_opt_depth:
                    c_list = find_paths(j, True)   # j side has coeff
                    p_list = find_paths(k, False)  # k side is pure
                    
                    for c in c_list:
                        for p in p_list:
                            nd = max(c['depth'], p['depth']) + 1
                            if nd == target_opt_depth:
                                new_ops = c['ops'] | p['ops']
                                new_ops.add(tuple(sorted((j, k))))
                                res.append({'ops': new_ops, 'depth': nd})
                
                # 2-2: ax^k * x^j (only when j != k, to avoid duplicates)
                if j != k:
                    if max(min_depths[k][1], min_depths[j][0]) + 1 <= target_opt_depth:
                        c_list = find_paths(k, True)   # k side has coeff
                        p_list = find_paths(j, False)  # j side is pure
                        
                        for c in c_list:
                            for p in p_list:
                                nd = max(c['depth'], p['depth']) + 1
                                if nd == target_opt_depth:
                                    new_ops = c['ops'] | p['ops']
                                    new_ops.add(tuple(sorted((j, k))))
                                    res.append({'ops': new_ops, 'depth': nd})
        
        memo[state_key] = res
        return res

    candidates = find_paths(n, multA)

    # ---------------------------------------------------------
    # 3. Result Processing -> Convert to XI objects
    # ---------------------------------------------------------
    if not candidates:
        return []

    # Filter only candidates with the minimum number of operations (add_count)
    min_ops_count = min(len(c['ops']) for c in candidates)
    
    final_xi_list = []
    seen_sets = set()  # Set to remove duplicate configurations (based on made_powers)

    for cand in candidates:
        ops_set = cand['ops']
        if len(ops_set) > min_ops_count:
            continue
        
        # Sort ops to create the route list
        # Sorting rule: ascending by result value (sum)
        # Since x^a * x^b = x^(a+b), sorting by a+b gives a topological order
        sorted_route = sorted(list(ops_set), key=lambda x: x[0] + x[1])
        
        # Build made_powers: set of 1 and all result powers from the route
        made_powers = {0, 1}
        for op in sorted_route:
            made_powers.add(op[0] + op[1])
            
        # Duplicate check (skip if configuration is identical)
        power_sig = tuple(sorted(list(made_powers)))
        if power_sig in seen_sets:
            continue
        seen_sets.add(power_sig)
        
        # Create XI object and append
        xi_obj = XI(multA, n)
        xi_obj.add_routes(sorted_route, made_powers)
        final_xi_list.append(xi_obj)
        
    return final_xi_list

def check_without_dcmp(made_powers: set[int], coeff: list[float]) -> bool:
    '''
    Check whether polynomial can be evaluated without decompositon.
    '''
    flag = True
    for i, coef in enumerate(coeff):
        if coef != 0:
            if i not in made_powers:
                flag = False
    
    return flag

def attach(d1: Poly | XI | None, c1: Complexity, d2: Poly, c2: Complexity, attach_type: str) -> Complexity:
    '''
    Attach calculation complexity
    '''
    res = Complexity()
    if attach_type == 'x': # Mult
        # d2 = constant
        if d2.deg == 0 and d2.coeff_type[0] == "F":
            res = c1
            res.depth += 1
            res.pmult += 1

        # normal poly*poly
        elif type(d1) == Poly:
            res.depth = max(c1.depth, c2.depth) + 1
            res.cmult = c1.cmult + c2.cmult + 1
            res.pmult = c1.pmult + c2.pmult
            res.add = c1.add + c2.add
        
        # (x^i) * poly
        elif type(d1) == XI:
            assert d1.n > 0
            res.depth = max(c1.depth, c2.depth) + 1
            res.cmult = c1.cmult + c2.cmult + 1
            res.pmult = c1.pmult + c2.pmult
            res.add = c1.add + c2.add
            
        
    elif attach_type == '+': # Add
        if d2.deg == 0:
            return c1
        res.depth = max(c1.depth, c2.depth)
        res.cmult = c1.cmult + c2.cmult
        res.pmult = c1.pmult + c2.pmult
        res.add = c1.add + c2.add + 1
        
    return res