# def is_multA_required(n: int) -> bool:
#     """
#     Check whether coefficient can be attached in x^n.
#     Depth shouldn't be increased.
#     """
#     if n == 0:
#         return False
#     return (n & (n - 1)) != 0

# def cal_multA(max_deg: int, n: int) -> bool:
#     if max_deg == n:
#         return False
#     xi_diff = ceil(log2(n)) < ceil(log2(n + 1))
#     p_diff = ceil(log2(max_deg - n)) < ceil(log2(max_deg - n + 1))
#     if xi_diff is False and p_diff is True:
#         return True
#     return False

# def solve_xn_routes(multA, n) -> list[XI]:
#     """
#     Find optimal route to make x^n.
#     """

#     # ---------------------------------------------------------
#     # 1. DP: Compute optimal Depth (parallel time) for all n
#     # ---------------------------------------------------------
#     # min_depths[i] = (pure_depth, coeff_depth)
#     min_depths = {1: (0, 0)} 
    
#     for i in range(2, n + 1):
#         d_pure, d_coeff = float('inf'), float('inf')
#         for j in range(1, i // 2 + 1):
#             k = i - j
            
#             # Pure (x^j * x^k = x^i)
#             # max time when both operands are ready + 1
#             p_val = max(min_depths[j][0], min_depths[k][0]) + 1
#             if p_val < d_pure: 
#                 d_pure = p_val
            
#             # Coeff (when multA=True, generate ax^i)
#             if multA:
#                 # Case 1: ax^j * x^k -> ax^i
#                 c_val1 = max(min_depths[j][1], min_depths[k][0]) + 1
#                 if c_val1 < d_coeff: 
#                     d_coeff = c_val1
                
#                 # Case 2: x^j * ax^k -> ax^i
#                 c_val2 = max(min_depths[k][1], min_depths[j][0]) + 1
#                 if c_val2 < d_coeff: 
#                     d_coeff = c_val2
                
#         min_depths[i] = (d_pure, d_coeff)

#     # Set target optimal Depth
#     target_opt_depth = min_depths[n][1] if multA else min_depths[n][0]

#     # ---------------------------------------------------------
#     # 2. Recursive Search: Backtracking and route construction
#     # ---------------------------------------------------------
#     # Return value: list of dict -> [{'ops': set((a,b)...), 'depth': int}, ...]
#     # ops is an unordered set of operations; later sorted to build the route
    
#     memo = {}

#     def find_paths(target, has_coeff):
#         state_key = (target, has_coeff)
#         if state_key in memo:
#             return memo[state_key]

#         # Base Case
#         if target == 1:
#             return [{'ops': set(), 'depth': 0}]
            
#         res = []
#         for j in range(1, target // 2 + 1):
#             k = target - j
            
#             # Case 1: Pure generation (x^target)
#             # Only performed when has_coeff is False
#             if not has_coeff:
#                 # Check critical path condition (pruning)
#                 if max(min_depths[j][0], min_depths[k][0]) + 1 <= min_depths[target][0]:
#                     l_list = find_paths(j, False)
#                     r_list = find_paths(k, False)
                    
#                     for l in l_list:
#                         for r in r_list:
#                             nd = max(l['depth'], r['depth']) + 1
#                             # Collect only paths that exactly match optimal Depth
#                             if nd == min_depths[target][0]:
#                                 new_ops = l['ops'] | r['ops']
#                                 # Add current operation (tuple ordered as (smaller, larger))
#                                 new_ops.add(tuple(sorted((j, k))))
#                                 res.append({
#                                     'ops': new_ops,
#                                     'depth': nd
#                                 })
            
#             # Case 2: Coeff generation (ax^target)
#             # Performed when has_coeff is True
#             else:
#                 # 2-1: ax^j * x^k
#                 if max(min_depths[j][1], min_depths[k][0]) + 1 <= target_opt_depth:
#                     c_list = find_paths(j, True)   # j side has coeff
#                     p_list = find_paths(k, False)  # k side is pure
                    
#                     for c in c_list:
#                         for p in p_list:
#                             nd = max(c['depth'], p['depth']) + 1
#                             if nd == target_opt_depth:
#                                 new_ops = c['ops'] | p['ops']
#                                 new_ops.add(tuple(sorted((j, k))))
#                                 res.append({'ops': new_ops, 'depth': nd})
                
#                 # 2-2: ax^k * x^j (only when j != k, to avoid duplicates)
#                 if j != k:
#                     if max(min_depths[k][1], min_depths[j][0]) + 1 <= target_opt_depth:
#                         c_list = find_paths(k, True)   # k side has coeff
#                         p_list = find_paths(j, False)  # j side is pure
                        
#                         for c in c_list:
#                             for p in p_list:
#                                 nd = max(c['depth'], p['depth']) + 1
#                                 if nd == target_opt_depth:
#                                     new_ops = c['ops'] | p['ops']
#                                     new_ops.add(tuple(sorted((j, k))))
#                                     res.append({'ops': new_ops, 'depth': nd})
        
#         memo[state_key] = res
#         return res

#     candidates = find_paths(n, multA)

#     # ---------------------------------------------------------
#     # 3. Result Processing -> Convert to XI objects
#     # ---------------------------------------------------------
#     if not candidates:
#         return []

#     # Filter only candidates with the minimum number of operations (add_count)
#     min_ops_count = min(len(c['ops']) for c in candidates)
    
#     final_xi_list = []
#     seen_sets = set()  # Set to remove duplicate configurations (based on made_powers)

#     for cand in candidates:
#         ops_set = cand['ops']
#         if len(ops_set) > min_ops_count:
#             continue
        
#         # Sort ops to create the route list
#         # Sorting rule: ascending by result value (sum)
#         # Since x^a * x^b = x^(a+b), sorting by a+b gives a topological order
#         sorted_route = sorted(list(ops_set), key=lambda x: x[0] + x[1])
        
#         # Build made_powers: set of 1 and all result powers from the route
#         made_powers = {0, 1}
#         for op in sorted_route:
#             made_powers.add(op[0] + op[1])
            
#         # Duplicate check (skip if configuration is identical)
#         power_sig = tuple(sorted(list(made_powers)))
#         if power_sig in seen_sets:
#             continue
#         seen_sets.add(power_sig)
        
#         # Create XI object and append
#         xi_obj = XI(multA, n)
#         xi_obj.add_routes(sorted_route, made_powers)
#         final_xi_list.append(xi_obj)
        
#     return final_xi_list

# def solve_xn_routes(multA, n, made_powers: set[int]={0, 1}) -> list[XI]:
#     """
#     Find optimal route to make x^n or ax^n.

#     made_powers:
#         Existing powers that have already been generated.

#     Design rule:
#         - returned route contains only newly required operations
#         - returned made_powers contains all input made_powers
#         - if n already exists in made_powers, route is empty
#     """

#     # ---------------------------------------------------------
#     # 0. Existing made powers
#     # ---------------------------------------------------------
#     input_made_powers = set(made_powers or set())

#     # Internally, 0 and 1 are always considered available.
#     # Output made_powers must still contain every input value.
#     base_made_powers = {0, 1} | input_made_powers

#     # Only powers in [1, n] are useful for constructing x^n.
#     available_for_dp = {
#         p for p in base_made_powers
#         if 1 <= p <= n
#     }

#     # If target already exists, no additional route is required.
#     if n in base_made_powers:
#         xi_obj = XI(multA, n)
#         xi_obj.add_routes([], set(base_made_powers))
#         return [xi_obj]

#     # ---------------------------------------------------------
#     # 1. DP: Compute optimal Depth using existing made_powers
#     # ---------------------------------------------------------
#     # min_depths[i] = (pure_depth, coeff_depth)
#     min_depths = {}

#     for i in range(1, n + 1):
#         if i in available_for_dp:
#             min_depths[i] = (0, 0)
#             continue

#         d_pure = float("inf")
#         d_coeff = float("inf")

#         for j in range(1, i // 2 + 1):
#             k = i - j

#             # Pure: x^j * x^k -> x^i
#             p_val = max(min_depths[j][0], min_depths[k][0]) + 1
#             if p_val < d_pure:
#                 d_pure = p_val

#             # Coeff: ax^j * x^k -> ax^i
#             # or     x^j * ax^k -> ax^i
#             if multA:
#                 c_val1 = max(min_depths[j][1], min_depths[k][0]) + 1
#                 if c_val1 < d_coeff:
#                     d_coeff = c_val1

#                 c_val2 = max(min_depths[k][1], min_depths[j][0]) + 1
#                 if c_val2 < d_coeff:
#                     d_coeff = c_val2

#         min_depths[i] = (d_pure, d_coeff)

#     # ---------------------------------------------------------
#     # 2. Recursive Search: route contains only new operations
#     # ---------------------------------------------------------
#     memo = {}

#     def find_paths(target, has_coeff):
#         state_key = (target, has_coeff)

#         if state_key in memo:
#             return memo[state_key]

#         # Existing power is reused directly.
#         # No route is added for it.
#         if target in available_for_dp:
#             return [{"ops": set(), "depth": 0}]

#         res = []

#         for j in range(1, target // 2 + 1):
#             k = target - j

#             if not has_coeff:
#                 target_depth = min_depths[target][0]

#                 if max(min_depths[j][0], min_depths[k][0]) + 1 <= target_depth:
#                     l_list = find_paths(j, False)
#                     r_list = find_paths(k, False)

#                     for l in l_list:
#                         for r in r_list:
#                             nd = max(l["depth"], r["depth"]) + 1

#                             if nd == target_depth:
#                                 new_ops = l["ops"] | r["ops"]

#                                 # Add only the operation that creates this target.
#                                 # Since target is not in available_for_dp here,
#                                 # this route does not overlap with existing powers.
#                                 new_ops.add(tuple(sorted((j, k))))

#                                 res.append({
#                                     "ops": new_ops,
#                                     "depth": nd,
#                                 })

#             else:
#                 target_depth = min_depths[target][1]

#                 # ax^j * x^k -> ax^target
#                 if max(min_depths[j][1], min_depths[k][0]) + 1 <= target_depth:
#                     c_list = find_paths(j, True)
#                     p_list = find_paths(k, False)

#                     for c in c_list:
#                         for p in p_list:
#                             nd = max(c["depth"], p["depth"]) + 1

#                             if nd == target_depth:
#                                 new_ops = c["ops"] | p["ops"]
#                                 new_ops.add(tuple(sorted((j, k))))

#                                 res.append({
#                                     "ops": new_ops,
#                                     "depth": nd,
#                                 })

#                 # ax^k * x^j -> ax^target
#                 if j != k:
#                     if max(min_depths[k][1], min_depths[j][0]) + 1 <= target_depth:
#                         c_list = find_paths(k, True)
#                         p_list = find_paths(j, False)

#                         for c in c_list:
#                             for p in p_list:
#                                 nd = max(c["depth"], p["depth"]) + 1

#                                 if nd == target_depth:
#                                     new_ops = c["ops"] | p["ops"]
#                                     new_ops.add(tuple(sorted((j, k))))

#                                     res.append({
#                                         "ops": new_ops,
#                                         "depth": nd,
#                                     })

#         memo[state_key] = res
#         return res

#     candidates = find_paths(n, multA)

#     # ---------------------------------------------------------
#     # 3. Result Processing -> Convert to XI objects
#     # ---------------------------------------------------------
#     if not candidates:
#         return []

#     # Keep only candidates with the minimum number of new operations.
#     min_ops_count = min(len(c["ops"]) for c in candidates)

#     final_xi_list = []
#     seen_sets = set()

#     for cand in candidates:
#         ops_set = cand["ops"]

#         if len(ops_set) > min_ops_count:
#             continue

#         # Route contains only new operations.
#         sorted_route = sorted(
#             list(ops_set),
#             key=lambda op: op[0] + op[1],
#         )

#         # Output made_powers includes all input powers and all newly generated powers.
#         updated_made_powers = set(base_made_powers)

#         for op in sorted_route:
#             updated_made_powers.add(op[0] + op[1])

#         power_sig = tuple(sorted(updated_made_powers))

#         if power_sig in seen_sets:
#             continue

#         seen_sets.add(power_sig)

#         xi_obj = XI(multA, n)
#         xi_obj.add_routes(sorted_route, updated_made_powers)
#         final_xi_list.append(xi_obj)

#     print(f"target: x^{n}, multA: {multA}, mp: {made_powers} -> {updated_made_powers}")
#     print(f"Found routes:")
#     for i, xi in enumerate(final_xi_list):
#         route = xi.route
#         print(f"{i}\t{' -> '.join([f"{a[0]}+{a[1]}" for a in route])}")
#     print()
#     return final_xi_list


def solve_xn_routes(multA, n, made_powers=None) -> list[XI]:
    """
    Find optimal route to make x^n or ax^n.

    made_powers:
        Existing pure powers x^k that have already been generated.

    Design rule:
        - made_powers tracks only pure powers x^k
        - ax^2 is treated as the initial coefficient-bearing term
        - coefficient operations do not add their target exponent to made_powers
        - returned route contains only newly required operations
        - returned made_powers contains all input made_powers
    """

    BASE_COEFF_POWER = 2

    # ---------------------------------------------------------
    # 0. Existing made powers
    # ---------------------------------------------------------
    input_made_powers = set(made_powers or set())

    # Internally, 0 and 1 are always considered available.
    base_made_powers = {0, 1} | input_made_powers

    # Pure powers available for DP.
    available_pure = {
        p for p in base_made_powers
        if 1 <= p <= n
    }

    # Coefficient-bearing powers available for DP.
    # Do NOT treat every made power as ax^p.
    available_coeff = set()
    if multA and BASE_COEFF_POWER in base_made_powers and BASE_COEFF_POWER <= n:
        available_coeff.add(BASE_COEFF_POWER)

    # Pure target already exists.
    if not multA and n in base_made_powers:
        xi_obj = XI(multA, n)
        xi_obj.add_routes([], set(base_made_powers))
        return [xi_obj]

    # Coefficient target already exists.
    if multA and n in available_coeff:
        xi_obj = XI(multA, n)
        xi_obj.add_routes([], set(base_made_powers))
        return [xi_obj]

    # ---------------------------------------------------------
    # 1. DP: Compute optimal depth
    # ---------------------------------------------------------
    # min_depths[i] = (pure_depth, coeff_depth)
    min_depths = {}

    for i in range(1, n + 1):
        d_pure = 0 if i in available_pure else float("inf")
        d_coeff = 0 if i in available_coeff else float("inf")

        # Pure: x^j * x^k -> x^i
        if d_pure != 0:
            for j in range(1, i // 2 + 1):
                k = i - j

                p_val = max(min_depths[j][0], min_depths[k][0]) + 1
                if p_val < d_pure:
                    d_pure = p_val

        # Coeff:
        # ax^j * x^k -> ax^i
        # x^j * ax^k -> ax^i
        if multA and d_coeff != 0:
            for j in range(1, i // 2 + 1):
                k = i - j

                c_val1 = max(min_depths[j][1], min_depths[k][0]) + 1
                if c_val1 < d_coeff:
                    d_coeff = c_val1

                if j != k:
                    c_val2 = max(min_depths[k][1], min_depths[j][0]) + 1
                    if c_val2 < d_coeff:
                        d_coeff = c_val2

        min_depths[i] = (d_pure, d_coeff)

    # ---------------------------------------------------------
    # 2. Recursive Search
    # ---------------------------------------------------------
    memo = {}

    def find_paths(target, has_coeff):
        state_key = (target, has_coeff)

        if state_key in memo:
            return memo[state_key]

        # Existing pure power is reused directly.
        if not has_coeff and target in available_pure:
            return [{"ops": set(), "depth": 0}]

        # Existing coefficient-bearing power is reused directly.
        if has_coeff and target in available_coeff:
            return [{"ops": set(), "depth": 0}]

        target_depth = min_depths[target][1 if has_coeff else 0]

        if target_depth == float("inf"):
            memo[state_key] = []
            return []

        res = []

        for j in range(1, target // 2 + 1):
            k = target - j

            if not has_coeff:
                # x^j * x^k -> x^target
                if max(min_depths[j][0], min_depths[k][0]) + 1 <= target_depth:
                    l_list = find_paths(j, False)
                    r_list = find_paths(k, False)

                    for l in l_list:
                        for r in r_list:
                            nd = max(l["depth"], r["depth"]) + 1

                            if nd == target_depth:
                                new_ops = l["ops"] | r["ops"]

                                # Internal route keeps operation type.
                                new_ops.add(("pure", *tuple(sorted((j, k)))))

                                res.append({
                                    "ops": new_ops,
                                    "depth": nd,
                                })

            else:
                # ax^j * x^k -> ax^target
                if max(min_depths[j][1], min_depths[k][0]) + 1 <= target_depth:
                    c_list = find_paths(j, True)
                    p_list = find_paths(k, False)

                    for c in c_list:
                        for p in p_list:
                            nd = max(c["depth"], p["depth"]) + 1

                            if nd == target_depth:
                                new_ops = c["ops"] | p["ops"]
                                new_ops.add(("coeff", *tuple(sorted((j, k)))))

                                res.append({
                                    "ops": new_ops,
                                    "depth": nd,
                                })

                # x^j * ax^k -> ax^target
                if j != k:
                    if max(min_depths[k][1], min_depths[j][0]) + 1 <= target_depth:
                        c_list = find_paths(k, True)
                        p_list = find_paths(j, False)

                        for c in c_list:
                            for p in p_list:
                                nd = max(c["depth"], p["depth"]) + 1

                                if nd == target_depth:
                                    new_ops = c["ops"] | p["ops"]
                                    new_ops.add(("coeff", *tuple(sorted((j, k)))))

                                    res.append({
                                        "ops": new_ops,
                                        "depth": nd,
                                    })

        memo[state_key] = res
        return res

    candidates = find_paths(n, multA)

    # ---------------------------------------------------------
    # 3. Result Processing -> Convert to XI objects
    # ---------------------------------------------------------
    if not candidates:
        return []

    # Keep only candidates with the minimum number of new operations.
    min_ops_count = min(len(c["ops"]) for c in candidates)

    final_xi_list = []
    seen = set()

    for cand in candidates:
        ops_set = cand["ops"]

        if len(ops_set) > min_ops_count:
            continue

        sorted_internal_route = sorted(
            list(ops_set),
            key=lambda op: (op[1] + op[2], op[0], op[1], op[2]),
        )

        # Output made_powers includes input powers and newly generated pure powers only.
        updated_made_powers = set(base_made_powers)

        for kind, a, b in sorted_internal_route:
            if kind == "pure":
                updated_made_powers.add(a + b)

        # XI.route keeps the original external format: [(a, b), ...]
        sorted_route = [
            (a, b)
            for _, a, b in sorted_internal_route
        ]

        route_sig = tuple(sorted_internal_route)
        power_sig = tuple(sorted(updated_made_powers))
        sig = (route_sig, power_sig)

        if sig in seen:
            continue

        seen.add(sig)

        xi_obj = XI(multA, n)
        xi_obj.add_routes(sorted_route, updated_made_powers)
        final_xi_list.append(xi_obj)
        
    print(f"target: x^{n}, multA: {multA}, mp: {made_powers} -> {updated_made_powers}")
    print(f"Found routes:")
    for i, xi in enumerate(final_xi_list):
        route = xi.route
        print(f"{i}\t{' -> '.join([f"{a[0]}+{a[1]}" for a in route])}")
    print()
    
    return final_xi_list