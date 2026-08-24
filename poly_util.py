# util.py
import random
from itertools import product
from math import log2, ceil
try:
    from basic_class import Complexity, Poly, XI, Decomp
except:
    from .basic_class import Complexity, Poly, XI, Decomp

def _make_int_coef(i: int) -> float:
    return float(i + 1)


def _make_float_coef(i: int) -> float:
    return float(i + 1) + random.randint(1, 9) / 10.0


def make_type_poly(poly_type: list[str]) -> list[float]:
    '''
    Create random polynomial by input coefficient type
    '''
    assert all(x in {"F", "I", "0"} for x in poly_type)

    poly = []
    for i, ctype in enumerate(poly_type):
        if ctype == "I":
            poly.append(_make_int_coef(i))
        elif ctype == "F":
            poly.append(_make_float_coef(i))
        elif ctype == "0":
            poly.append(0.0)

    return poly

def make_all_polys(max_deg: int) -> list[list[float]]:
    '''
    Create every possible polynomial of degree n.
    Coefficient type can be 0/I/F.(Max degree can't be 0.)
    '''
    assert max_deg >= 1

    types_middle = ["0", "I", "F"]
    types_highest = ["I", "F"]

    iterables = [types_middle] * max_deg + [types_highest]

    type_combinations = product(*iterables)

    all_polys = []
    for combo in type_combinations:
        poly_values = []

        for i, ctype in enumerate(combo):
            if ctype == "0":
                poly_values.append(0.0)
            elif ctype == "I":
                poly_values.append(_make_int_coef(i))
            elif ctype == "F":
                poly_values.append(_make_float_coef(i))

        all_polys.append(poly_values)

    return all_polys

def make_deg_poly(deg: int, type: str) -> list[float]:
    """
    deg차 다항식 생성
    계수 순서: [상수항, x항, x^2항, ..., x^deg항]
    """

    assert deg >= 1
    assert type in {"all", "even", "odd", "random"}

    poly = []

    if type == "all":
        poly = [_make_float_coef(i) for i in range(deg + 1)]

    elif type == "even":
        if deg % 2 != 0:
            raise ValueError("짝수다항식의 최고차항 deg는 짝수여야 합니다.")

        for i in range(deg + 1):
            if i % 2 == 0:
                poly.append(_make_float_coef(i))
            else:
                poly.append(0.0)

    elif type == "odd":
        if deg % 2 == 0:
            raise ValueError("홀수다항식의 최고차항 deg는 홀수여야 합니다.")

        for i in range(deg + 1):
            if i % 2 == 1:
                poly.append(_make_float_coef(i))
            else:
                poly.append(0.0)

    elif type == "random":
        for i in range(deg):
            ctype = random.choice(["0", "I", "F"])

            if ctype == "0":
                poly.append(0.0)
            elif ctype == "I":
                poly.append(_make_int_coef(i))
            elif ctype == "F":
                poly.append(_make_float_coef(i))

        highest_type = random.choice(["I", "F"])

        if highest_type == "I":
            poly.append(_make_int_coef(deg))
        elif highest_type == "F":
            poly.append(_make_float_coef(deg))

    return poly

def print_dcmp_detail(dcmp: Decomp, step: int):
    print(f"Step {step}")
    print(f"\t"*step, end='')
    print(f"dcmp: {dcmp.restore_dcmp()}")
    print(f"\t"*step, end='')
    print(f"mp: {dcmp.made_powers}")
    print(f"\t"*step, end='')
    dcmp.comp.print_params()
    if dcmp.xi.n != 0:
        print(f"\t"*step, end='')
        dcmp.xi.print_params()
    
    if dcmp.dcmp_p is not None:
        print_dcmp_detail(dcmp.dcmp_p, step+1)
    if dcmp.dcmp_q is not None:
        print_dcmp_detail(dcmp.dcmp_q, step+1)


def solve_xn_routes(multA, n, made_powers=None) -> list[XI]:
    """
    기존에 생성한 x^i들을 활용하여 (a)x^n을 구성하는 최소 경로들을 탐색한다.

    Rules:
        - made_powers는 순수 power x^k만 추적한다.
        - coefficient route는 순수 x^(a+b)를 만들지 않으므로 made_powers에 추가하지 않는다.
        - multA=True이고 x^n이 이미 존재하면 a * x^n도 고려한다.
        - multA=True이고 x^n이 아직 없더라도, x^n을 만든 뒤 a * x^n 하는 fallback도 고려한다.
        - 최종적으로 같은 route/depth/pmult/cmult 후보가 있으면 made_powers가 더 큰 후보만 남긴다.
    """
    from math import ceil, log2

    BASE_COEFF_POWER = 2

    def ceil_log2(v: int) -> int:
        if v <= 1:
            return 0
        return ceil(log2(v))

    def pure_depth(k: int) -> int:
        return ceil_log2(k)

    def direct_coeff_depth(k: int) -> int:
        # a * x^k
        return ceil_log2(k) + 1

    def base_coeff_depth(k: int) -> int:
        # ax^k as coefficient-bearing base depth
        return ceil_log2(k + 1)

    # ---------------------------------------------------------
    # 0. Existing made powers
    # ---------------------------------------------------------
    input_made_powers = set(made_powers or set())
    base_made_powers = {0, 1} | input_made_powers

    # 재사용 가능한 순수 x^i 집합.
    available_pure = {
        p for p in base_made_powers
        if 1 <= p <= n
    }

    # 사용 가능한 coefficient-bearing base 집합.
    available_coeff_base = set()
    if multA and BASE_COEFF_POWER in base_made_powers and BASE_COEFF_POWER <= n:
        available_coeff_base.add(BASE_COEFF_POWER)

    # 이미 x^n이 있고, 순수 x^n만 필요하면 바로 반환.
    if not multA and n in base_made_powers:
        xi_obj = XI(multA, n)
        xi_obj.add_routes(
            [],
            pure_depth(n),
            0,
            0,
            set(base_made_powers),
            route_ops=[],
        )
        return [xi_obj]

    # ---------------------------------------------------------
    # 1. DP: 최소 depth 계산
    # ---------------------------------------------------------
    # min_depths[i] = (pure_depth: x^i 비용, coeff_depth: ax^i 비용)
    min_depths = {}

    for i in range(1, n + 1):
        # Pure x^i
        d_pure = pure_depth(i) if i in available_pure else float("inf")

        for j in range(1, i // 2 + 1):
            k = i - j

            if j in min_depths and k in min_depths:
                p_left = min_depths[j][0]
                p_right = min_depths[k][0]

                if p_left != float("inf") and p_right != float("inf"):
                    p_val = max(p_left, p_right) + 1
                    if p_val < d_pure:
                        d_pure = p_val

        # Coefficient ax^i
        d_coeff = float("inf")

        if multA:
            # ax = a * x
            if i == 1:
                d_coeff = min(d_coeff, 1)

            # 이미 순수 x^2가 만들어져 있으면 ax^2를 base로 사용할 수 있다.
            if i in available_coeff_base:
                d_coeff = min(d_coeff, base_coeff_depth(i))

            # 목표 x^n이 이미 있으면 a * x^n을 고려한다.
            if i == n and i in available_pure:
                d_coeff = min(d_coeff, direct_coeff_depth(i))

            # 목표 x^n이 아직 없어도, 순수 x^n을 만든 뒤 a * x^n을 수행할 수 있다.
            if i == n and i not in available_pure and d_pure != float("inf"):
                d_coeff = min(d_coeff, d_pure + 1)

            for j in range(1, i // 2 + 1):
                k = i - j

                if j in min_depths and k in min_depths:
                    # ax^j * x^k -> ax^i
                    c_left = min_depths[j][1]
                    p_right = min_depths[k][0]

                    if c_left != float("inf") and p_right != float("inf"):
                        c_val = max(c_left, p_right) + 1
                        if c_val < d_coeff:
                            d_coeff = c_val

                    # ax^k * x^j -> ax^i
                    if j != k:
                        c_right = min_depths[k][1]
                        p_left = min_depths[j][0]

                        if c_right != float("inf") and p_left != float("inf"):
                            c_val = max(c_right, p_left) + 1
                            if c_val < d_coeff:
                                d_coeff = c_val

        min_depths[i] = (d_pure, d_coeff)

    # ---------------------------------------------------------
    # 2. Recursive Search: 최소 depth를 만족하는 route 복원
    # ---------------------------------------------------------
    memo = {}

    def find_paths(target, has_coeff):
        state_key = (target, has_coeff)

        if state_key in memo:
            return memo[state_key]

        target_depth = min_depths[target][1 if has_coeff else 0]

        if target_depth == float("inf"):
            memo[state_key] = []
            return []

        res = []

        if not has_coeff:
            # 이미 존재하는 순수 power는 재사용한다.
            if target in available_pure:
                d = pure_depth(target)

                if d == target_depth:
                    res.append({
                        "ops": set(),
                        "depth": d,
                    })

                memo[state_key] = res
                return res

            # x^j * x^k -> x^target
            for j in range(1, target // 2 + 1):
                k = target - j

                if max(min_depths[j][0], min_depths[k][0]) + 1 > target_depth:
                    continue

                l_list = find_paths(j, False)
                r_list = find_paths(k, False)

                for l in l_list:
                    for r in r_list:
                        nd = max(l["depth"], r["depth"]) + 1

                        if nd == target_depth:
                            a, b = sorted((j, k))
                            new_ops = l["ops"] | r["ops"]
                            new_ops.add(("pure", a, b))

                            res.append({
                                "ops": new_ops,
                                "depth": nd,
                            })

        else:
            # ax = a * x
            if target == 1:
                d = 1

                if d == target_depth:
                    res.append({
                        "ops": set(),
                        "depth": d,
                    })

                    memo[state_key] = res
                    return res

            # 이미 coefficient-bearing base로 사용 가능한 ax^target.
            if target in available_coeff_base:
                d = base_coeff_depth(target)

                if d == target_depth:
                    res.append({
                        "ops": set(),
                        "depth": d,
                    })

            # 이미 x^n이 있으면 a * x^n을 직접 수행한다.
            if target == n and target in available_pure:
                d = direct_coeff_depth(target)

                if d == target_depth:
                    res.append({
                        "ops": {("coeff_direct", 0, target)},
                        "depth": d,
                    })

            # 순수 x^target을 먼저 만든 뒤, 마지막에 a * x^target 수행.
            if target == n and target not in available_pure:
                pure_target_depth = min_depths[target][0]

                if pure_target_depth != float("inf") and pure_target_depth + 1 == target_depth:
                    pure_paths = find_paths(target, False)

                    for p in pure_paths:
                        if p["depth"] + 1 == target_depth:
                            res.append({
                                "ops": set(p["ops"]),
                                "depth": p["depth"] + 1,
                            })

            for j in range(1, target // 2 + 1):
                k = target - j

                # ax^j * x^k -> ax^target
                if max(min_depths[j][1], min_depths[k][0]) + 1 <= target_depth:
                    c_list = find_paths(j, True)
                    p_list = find_paths(k, False)

                    for c in c_list:
                        for p in p_list:
                            nd = max(c["depth"], p["depth"]) + 1

                            if nd == target_depth:
                                new_ops = c["ops"] | p["ops"]
                                new_ops.add(("coeff", j, k))

                                res.append({
                                    "ops": new_ops,
                                    "depth": nd,
                                })

                # ax^k * x^j -> ax^target
                if j != k:
                    if max(min_depths[k][1], min_depths[j][0]) + 1 <= target_depth:
                        c_list = find_paths(k, True)
                        p_list = find_paths(j, False)

                        for c in c_list:
                            for p in p_list:
                                nd = max(c["depth"], p["depth"]) + 1

                                if nd == target_depth:
                                    new_ops = c["ops"] | p["ops"]
                                    new_ops.add(("coeff", k, j))

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

    min_depth = min(c["depth"] for c in candidates)
    candidates = [
        c for c in candidates
        if c["depth"] == min_depth
    ]

    min_ops_count = min(len(c["ops"]) for c in candidates)

    kind_order = {
        "pure": 0,
        "coeff_direct": 1,
        "coeff": 2,
    }

    converted_candidates = []

    for cand in candidates:
        ops_set = cand["ops"]

        if len(ops_set) > min_ops_count:
            continue

        sorted_internal_route = sorted(
            list(ops_set),
            key=lambda op: (
                op[1] + op[2],
                kind_order.get(op[0], 99),
                op[1],
                op[2],
            ),
        )

        updated_made_powers = set(base_made_powers)

        for kind, a, b in sorted_internal_route:
            # 순수 power 생성 route만 made_powers를 증가시킨다.
            if kind == "pure":
                updated_made_powers.add(a + b)

        route_ops = [
            (kind, a, b)
            for kind, a, b in sorted_internal_route
        ]
        route = [(a, b) for _, a, b in route_ops]

        depth = cand["depth"]
        pmult = 1 if multA else 0

        # coeff_direct는 a * x^k 이므로 CMult에서 제외한다.
        cmult = sum(
            1
            for kind, _, _ in sorted_internal_route
            if kind != "coeff_direct"
        )

        converted_candidates.append({
            "route": route,
            "route_ops": route_ops,
            "depth": depth,
            "pmult": pmult,
            "cmult": cmult,
            "made_powers": updated_made_powers,
        })

    # 같은 외부 표현 route/depth/pmult/cmult를 갖는 후보 중,
    # made_powers가 strict superset인 후보가 있으면 dominated 후보를 제거한다.
    pruned_candidates = []

    for idx, item in enumerate(converted_candidates):
        dominated = False

        for jdx, other in enumerate(converted_candidates):
            if idx == jdx:
                continue

            same_external_route = (
                tuple(item["route_ops"]) == tuple(other["route_ops"])
                and item["depth"] == other["depth"]
                and item["pmult"] == other["pmult"]
                and item["cmult"] == other["cmult"]
            )

            if same_external_route and other["made_powers"] > item["made_powers"]:
                dominated = True
                break

        if not dominated:
            pruned_candidates.append(item)

    final_xi_list = []
    seen = set()

    for item in pruned_candidates:
        sig = (
            tuple(item["route_ops"]),
            tuple(sorted(item["made_powers"])),
            item["depth"],
            item["pmult"],
            item["cmult"],
        )

        if sig in seen:
            continue

        seen.add(sig)

        xi_obj = XI(multA, n)
        xi_obj.add_routes(
            item["route"],
            item["depth"],
            item["pmult"],
            item["cmult"],
            item["made_powers"],
            route_ops=item["route_ops"],
        )
        final_xi_list.append(xi_obj)

    return final_xi_list

# def solve_xn_routes(multA, n, made_powers=None) -> list[XI]:
#     """
#     기존에 생성한 x^i들을 활용하여 (a)x^n을 구성하는 최소의 경로들 탐색
#     made_powers: 기존에 생성한 x^i들.

#     Rules:
#         - made_powers tracks only pure powers x^k.
#         - coefficient routes do not add their target exponent to made_powers.
#         - if multA=True and x^n already exists, a * x^n is also considered.
#         - XI.add_routes(route, depth, pmult, cmult, made_powers)
#     """
#     from math import ceil, log2

#     BASE_COEFF_POWER = 2

#     def ceil_log2(v: int) -> int:
#         if v <= 1:
#             return 0
#         return ceil(log2(v))

#     def pure_depth(k: int) -> int:
#         return ceil_log2(k)

#     def direct_coeff_depth(k: int) -> int:
#         # a * x^k
#         return ceil_log2(k) + 1

#     def base_coeff_depth(k: int) -> int:
#         # ax^k as coefficient-bearing base depth
#         return ceil_log2(k + 1)

#     # ---------------------------------------------------------
#     # 0. Existing made powers
#     # ---------------------------------------------------------
#     input_made_powers = set(made_powers or set())

#     # Internally, 0 and 1 are always considered available.
#     base_made_powers = {0, 1} | input_made_powers

#     # 재사용이 가능한 x^i들의 집합.
#     available_pure = {
#         p for p in base_made_powers
#         if 1 <= p <= n
#     }

#     # 사용가능한 ax^i들의 집합.
#     available_coeff_base = set()
#     # x^2가 이미 만들어진 경우에는 ax^2를 사용할 수 있음.
#     if multA and BASE_COEFF_POWER in base_made_powers and BASE_COEFF_POWER <= n:
#         available_coeff_base.add(BASE_COEFF_POWER)

#     # 완성된 x^n이 이미 존재하고, multA=False인 경우에는 바로 반환.
#     if not multA and n in base_made_powers:
#         xi_obj = XI(multA, n)
#         xi_obj.add_routes(
#             [],
#             pure_depth(n),
#             0,
#             0,
#             set(base_made_powers),
#         )
#         return [xi_obj]

#     # ---------------------------------------------------------
#     # 1. DP: 목표로 하는 최적 depth를 연산.
#     # ---------------------------------------------------------
#     # min_depths[i] = (pure_depth: x^n 연산비용, coeff_depth: ax^n 연산비용)
#     min_depths = {}

#     for i in range(1, n + 1):
#         # Pure x^i
#         # 초기 depth는 이미 만들어진 경우 log2(i), 새로 만들어야하는 경우 inf.
#         d_pure = pure_depth(i) if i in available_pure else float("inf")

#         # x^i = x^j * x^k로 구성할 수 있는지 확인.
#         # j, k 모두 만들어진 경우 [max(depth(j), depth(k)) + 1] 로 x^i의 depth 갱신.
#         for j in range(1, i // 2 + 1):
#             k = i - j

#             if j in min_depths and k in min_depths:
#                 p_left = min_depths[j][0]
#                 p_right = min_depths[k][0]

#                 if p_left != float("inf") and p_right != float("inf"):
#                     p_val = max(p_left, p_right) + 1
#                     if p_val < d_pure:
#                         d_pure = p_val

#         # Coefficient ax^i
#         d_coeff = float("inf")

#         if multA:
            
#             if i == 1:
#                 d_coeff = min(d_coeff, 1)
                
#             # Existing coefficient-bearing base: ax^2
#             if i in available_coeff_base:
#                 d_coeff = min(d_coeff, base_coeff_depth(i))

#             # If target x^n already exists, consider a * x^n.
#             if i == n and i in available_pure:
#                 d_coeff = min(d_coeff, direct_coeff_depth(i))

#             # x^n이 없더라도 x^n을 만들 수 있다면 이후 multA를 수행할 수도 있다.
#             if i == n and i not in available_pure and d_pure != float("inf"):
#                 d_coeff = min(d_coeff, d_pure + 1)
        
#             for j in range(1, i // 2 + 1):
#                 k = i - j

#                 if j in min_depths and k in min_depths:
#                     # ax^j * x^k -> ax^i
#                     c_left = min_depths[j][1]
#                     p_right = min_depths[k][0]

#                     if c_left != float("inf") and p_right != float("inf"):
#                         c_val = max(c_left, p_right) + 1
#                         if c_val < d_coeff:
#                             d_coeff = c_val

#                     # ax^k * x^j -> ax^i
#                     if j != k:
#                         c_right = min_depths[k][1]
#                         p_left = min_depths[j][0]

#                         if c_right != float("inf") and p_left != float("inf"):
#                             c_val = max(c_right, p_left) + 1
#                             if c_val < d_coeff:
#                                 d_coeff = c_val

#         min_depths[i] = (d_pure, d_coeff)

#     # ---------------------------------------------------------
#     # 2. Recursive Search
#     # ---------------------------------------------------------
#     memo = {}

#     def find_paths(target, has_coeff):
#         state_key = (target, has_coeff)

#         if state_key in memo:
#             return memo[state_key]

#         target_depth = min_depths[target][1 if has_coeff else 0]

#         if target_depth == float("inf"):
#             memo[state_key] = []
#             return []

#         res = []

#         if not has_coeff:
#             # Existing pure power is reused directly.
#             if target in available_pure:
#                 d = pure_depth(target)

#                 if d == target_depth:
#                     res.append({
#                         "ops": set(),
#                         "depth": d,
#                     })

#                 memo[state_key] = res
#                 return res

#             # x^j * x^k -> x^target
#             for j in range(1, target // 2 + 1):
#                 k = target - j

#                 if max(min_depths[j][0], min_depths[k][0]) + 1 > target_depth:
#                     continue

#                 l_list = find_paths(j, False)
#                 r_list = find_paths(k, False)

#                 for l in l_list:
#                     for r in r_list:
#                         nd = max(l["depth"], r["depth"]) + 1

#                         if nd == target_depth:
#                             a, b = sorted((j, k))
#                             new_ops = l["ops"] | r["ops"]
#                             new_ops.add(("pure", a, b))

#                             res.append({
#                                 "ops": new_ops,
#                                 "depth": nd,
#                             })

#         else:
#             if target == 1:
#                 d = 1

#                 if d == target_depth:
#                     res.append({
#                         "ops": set(),
#                         "depth": d,
#                     })

#                     memo[state_key] = res
#                     return res
        
#             # Existing coefficient-bearing base: ax^2
#             if target in available_coeff_base:
#                 d = base_coeff_depth(target)

#                 if d == target_depth:
#                     res.append({
#                         "ops": set(),
#                         "depth": d,
#                     })

#             # Direct target route: a * x^n
#             if target == n and target in available_pure:
#                 d = direct_coeff_depth(target)

#                 if d == target_depth:
#                     res.append({
#                         "ops": {("coeff_direct", 0, target)},
#                         "depth": d,
#                     })

#             # x^target을 먼저 만든 뒤, 마지막에 plaintext multiplication으로 a*x^target 생성.
#             # route에는 pure-power 생성 경로만 남기고,
#             # pmult는 최종 XI 생성부의 pmult = 1 if multA else 0에서 반영한다.
#             if target == n and target not in available_pure:
#                 pure_target_depth = min_depths[target][0]

#                 if pure_target_depth != float("inf") and pure_target_depth + 1 == target_depth:
#                     pure_paths = find_paths(target, False)

#                     for p in pure_paths:
#                         if p["depth"] + 1 == target_depth:
#                             res.append({
#                                 "ops": set(p["ops"]),
#                                 "depth": p["depth"] + 1,
#                             })
                    
#             for j in range(1, target // 2 + 1):
#                 k = target - j

#                 # ax^j * x^k -> ax^target
#                 if max(min_depths[j][1], min_depths[k][0]) + 1 <= target_depth:
#                     c_list = find_paths(j, True)
#                     p_list = find_paths(k, False)

#                     for c in c_list:
#                         for p in p_list:
#                             nd = max(c["depth"], p["depth"]) + 1

#                             if nd == target_depth:
#                                 new_ops = c["ops"] | p["ops"]
#                                 new_ops.add(("coeff", j, k))

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
#                                     new_ops.add(("coeff", k, j))

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

#     min_depth = min(c["depth"] for c in candidates)
#     candidates = [
#         c for c in candidates
#         if c["depth"] == min_depth
#     ]

#     min_ops_count = min(len(c["ops"]) for c in candidates)

#     final_xi_list = []
#     filtered = {}

#     for cand in candidates:
#         ops_set = cand["ops"]

#         if len(ops_set) > min_ops_count:
#             continue

#         sorted_internal_route = sorted(
#             list(ops_set),
#             key=lambda op: (
#                 op[1] + op[2],
#                 kind_order.get(op[0], 99),
#                 op[1],
#                 op[2],
#             ),
#         )

#         updated_made_powers = set(base_made_powers)

#         for kind, a, b in sorted_internal_route:
#             if kind == "pure":
#                 updated_made_powers.add(a + b)

#         route = [
#             (a, b)
#             for _, a, b in sorted_internal_route
#         ]

#         depth = cand["depth"]
#         pmult = 1 if multA else 0

#         cmult = sum(
#             1
#             for kind, _, _ in sorted_internal_route
#             if kind != "coeff_direct"
#         )

#         # kind를 제거한 외부 표현 기준으로 같은 후보끼리 비교한다.
#         key = (
#             tuple(route),
#             depth,
#             pmult,
#             cmult,
#         )

#         prev = filtered.get(key)

#         # 같은 route/depth/pmult/cmult라면 made_powers가 더 많은 후보를 유지한다.
#         # 특히 solve_xn_routes(True, 2, {0, 1})에서
#         # {0, 1} 후보 대신 {0, 1, 2} 후보를 남긴다.
#         if prev is None:
#             filtered[key] = {
#                 "route": route,
#                 "depth": depth,
#                 "pmult": pmult,
#                 "cmult": cmult,
#                 "made_powers": updated_made_powers,
#             }
#         else:
#             prev_mp = prev["made_powers"]

#             if updated_made_powers > prev_mp:
#                 filtered[key] = {
#                     "route": route,
#                     "depth": depth,
#                     "pmult": pmult,
#                     "cmult": cmult,
#                     "made_powers": updated_made_powers,
#                 }
#             elif not (updated_made_powers < prev_mp):
#                 # 서로 strict subset 관계가 아니면 보수적으로 둘 다 살려야 하지만,
#                 # 현재 filtered는 key당 하나만 저장하므로 더 큰 집합 크기를 우선한다.
#                 if len(updated_made_powers) > len(prev_mp):
#                     filtered[key] = {
#                         "route": route,
#                         "depth": depth,
#                         "pmult": pmult,
#                         "cmult": cmult,
#                         "made_powers": updated_made_powers,
#                     }

#     seen = set()

#     for item in filtered.values():
#         sig = (
#             tuple(item["route"]),
#             tuple(sorted(item["made_powers"])),
#             item["depth"],
#             item["pmult"],
#             item["cmult"],
#         )

#         if sig in seen:
#             continue

#         seen.add(sig)

#         xi_obj = XI(multA, n)
#         xi_obj.add_routes(
#             item["route"],
#             item["depth"],
#             item["pmult"],
#             item["cmult"],
#             item["made_powers"],
#         )
#         final_xi_list.append(xi_obj)

#     return final_xi_list

#     final_xi_list = []
#     seen = set()

#     kind_order = {
#         "pure": 0,
#         "coeff_direct": 1,
#         "coeff": 2,
#     }

#     for cand in candidates:
#         ops_set = cand["ops"]

#         if len(ops_set) > min_ops_count:
#             continue

#         sorted_internal_route = sorted(
#             list(ops_set),
#             key=lambda op: (
#                 op[1] + op[2],
#                 kind_order.get(op[0], 99),
#                 op[1],
#                 op[2],
#             ),
#         )

#         updated_made_powers = set(base_made_powers)

#         for kind, a, b in sorted_internal_route:
#             if kind == "pure":
#                 updated_made_powers.add(a + b)

#         route = [
#             (a, b)
#             for _, a, b in sorted_internal_route
#         ]

#         depth = cand["depth"]
#         pmult = 1 if multA else 0

#         # (0, k)는 a * x^k 를 의미하므로 coefficient multiplication count에서 제외한다.
#         cmult = sum(
#             1
#             for kind, _, _ in sorted_internal_route
#             if kind != "coeff_direct"
#         )

#         sig = (
#             tuple(sorted_internal_route),
#             tuple(sorted(updated_made_powers)),
#             depth,
#             pmult,
#             cmult,
#         )

#         if sig in seen:
#             continue

#         seen.add(sig)

#         xi_obj = XI(multA, n)
#         xi_obj.add_routes(
#             route,
#             depth,
#             pmult,
#             cmult,
#             updated_made_powers,
#         )
#         final_xi_list.append(xi_obj)
    
#     # if(multA):
#     #     print(f"target: x^{n}, multA: {multA}, mp: {made_powers} -> {updated_made_powers}")
#     #     print(f"Found routes:")
#     #     for i, xi in enumerate(final_xi_list):
#     #         route = xi.route
#     #         xi.print_params()
#     #         print(f"{i}\t{' -> '.join([f"{a[0]}+{a[1]}" for a in route])}")
#     #     print()
    
#     return final_xi_list


def check_without_dcmp(made_powers: set[int], coeff: list[float]) -> bool:
    '''
    Check whether polynomial can be evaluated without cmult.
    '''
    for i, coef in enumerate(coeff):
        if coef != 0:
            if i not in made_powers:
                return False
    return True

def check_without_dcmp_v2(
    made_powers: set[int],
    coeff: list[float],
) -> tuple[
    list[tuple[int, int]] | None,
    list[dict] | None,
    Complexity | None,
    set[int] | None,
]:
    """
    추가 분해 없이 다항식을 평가하는 최적 terminal 계획을 생성한다.

    반환값:
        total_routes: 기존 코드와의 호환성을 위한 (lhs, rhs) 목록
        term_plans: 항별 평가계획
        complexity: terminal 평가 복잡도
        made_powers: 평가 후 사용 가능한 순수 power 집합
    """
    def ceil_log2(v: int) -> int:
        if v <= 1:
            return 0
        return ceil(log2(v))

    def make_term_plan(
        degree: int,
        coeff_type: str,
        multA: bool,
        route_ops,
    ) -> dict:
        return {
            "degree": int(degree),
            "coeff_type": str(coeff_type),
            "multA": bool(multA),
            "route_ops": [
                (str(kind), int(a), int(b))
                for kind, a, b in route_ops
            ],
        }

    class TermCandidate:
        def __init__(
            self,
            degree,
            coeff_type,
            multA,
            route_ops,
            depth,
            cmult,
            pmult,
            made_powers,
        ):
            self.term_plan = make_term_plan(
                degree,
                coeff_type,
                multA,
                route_ops,
            )
            self.depth = depth
            self.cmult = cmult
            self.pmult = pmult
            self.made_powers = set(made_powers)

    def plan_key(plan: dict):
        return (
            plan["degree"],
            plan["coeff_type"],
            int(plan["multA"]),
            tuple(plan["route_ops"]),
        )

    def candidate_key(cand: TermCandidate):
        return (
            cand.depth,
            cand.cmult,
            cand.pmult,
            plan_key(cand.term_plan),
            -len(cand.made_powers),
            tuple(sorted(cand.made_powers)),
        )

    def result_key(item):
        term_plans, comp, mp = item
        return (
            comp.depth,
            comp.cmult,
            comp.pmult,
            comp.add,
            tuple(plan_key(plan) for plan in term_plans),
            -len(mp),
            tuple(sorted(mp)),
        )

    base_made_powers = {0, 1} | set(made_powers or set())
    poly = Poly(coeff)

    required_powers = [
        i for i, coef in enumerate(poly.coeff)
        if i >= 1 and coef != 0
    ]

    add_count = len(poly.coeff_type) - poly.coeff_type.count("0") - 1

    memo: dict[
        tuple[int, tuple[int, ...]],
        tuple[list[dict], Complexity, set[int]] | None,
    ] = {}

    def copy_result(item):
        if item is None:
            return None
        term_plans, comp, mp = item
        copied_comp = Complexity()
        copied_comp.insert_value(
            comp.depth,
            comp.cmult,
            comp.pmult,
            comp.add,
        )
        copied_plans = [
            {
                "degree": plan["degree"],
                "coeff_type": plan["coeff_type"],
                "multA": plan["multA"],
                "route_ops": list(plan["route_ops"]),
            }
            for plan in term_plans
        ]
        return copied_plans, copied_comp, set(mp)

    def solve_from(
        pos: int,
        current_made_powers: set[int],
    ) -> tuple[list[dict], Complexity, set[int]] | None:
        state_key = (pos, tuple(sorted(current_made_powers)))

        if state_key in memo:
            return copy_result(memo[state_key])

        if pos >= len(required_powers):
            comp = Complexity()
            comp.insert_value(0, 0, 0, 0)
            res = ([], comp, set(current_made_powers))
            memo[state_key] = res
            return copy_result(res)

        target_power = required_powers[pos]
        coeff_type = poly.coeff_type[target_power]
        term_candidates: list[TermCandidate] = []

        if target_power in current_made_powers:
            term_candidates.append(
                TermCandidate(
                    degree=target_power,
                    coeff_type=coeff_type,
                    multA=False,
                    route_ops=[],
                    depth=ceil_log2(target_power)
                    + (1 if coeff_type == "F" else 0),
                    cmult=0,
                    pmult=1 if coeff_type == "F" else 0,
                    made_powers=current_made_powers,
                )
            )
        elif coeff_type == "I":
            for xi in solve_xn_routes(
                False,
                target_power,
                current_made_powers,
            ):
                term_candidates.append(
                    TermCandidate(
                        degree=target_power,
                        coeff_type=coeff_type,
                        multA=False,
                        route_ops=xi.route_ops,
                        depth=xi.depth,
                        cmult=xi.cmult,
                        pmult=xi.pmult,
                        made_powers=xi.made_powers,
                    )
                )
        elif coeff_type == "F":
            # 순수 x^i를 만든 뒤 마지막에 계수를 곱하는 경로
            for xi in solve_xn_routes(
                False,
                target_power,
                current_made_powers,
            ):
                term_candidates.append(
                    TermCandidate(
                        degree=target_power,
                        coeff_type=coeff_type,
                        multA=False,
                        route_ops=xi.route_ops,
                        depth=xi.depth + 1,
                        cmult=xi.cmult,
                        pmult=xi.pmult + 1,
                        made_powers=xi.made_powers,
                    )
                )

            # 계수를 먼저 결합한 coefficient-bearing 경로
            for xi in solve_xn_routes(
                True,
                target_power,
                current_made_powers,
            ):
                term_candidates.append(
                    TermCandidate(
                        degree=target_power,
                        coeff_type=coeff_type,
                        multA=True,
                        route_ops=xi.route_ops,
                        depth=xi.depth,
                        cmult=xi.cmult,
                        pmult=xi.pmult,
                        made_powers=xi.made_powers,
                    )
                )

        if not term_candidates:
            memo[state_key] = None
            return None

        term_candidates.sort(key=candidate_key)
        best_result = None

        for cand in term_candidates:
            next_made_powers = (
                set(current_made_powers)
                | set(cand.made_powers)
            )

            rest = solve_from(pos + 1, next_made_powers)
            if rest is None:
                continue

            rest_plans, rest_comp, rest_made_powers = rest

            comp = Complexity()
            comp.depth = max(cand.depth, rest_comp.depth)
            comp.cmult = cand.cmult + rest_comp.cmult
            comp.pmult = cand.pmult + rest_comp.pmult
            comp.add = 0

            term_plans = [cand.term_plan] + list(rest_plans)
            final_made_powers = (
                set(next_made_powers)
                | set(rest_made_powers)
            )

            item = (term_plans, comp, final_made_powers)

            if (
                best_result is None
                or result_key(item) < result_key(best_result)
            ):
                best_result = item

        memo[state_key] = best_result
        return copy_result(best_result)

    result = solve_from(0, set(base_made_powers))

    if result is None:
        return None, None, None, None

    term_plans, comp_res, final_made_powers = result

    if poly.coeff and poly.coeff[0] != 0:
        term_plans.insert(
            0,
            make_term_plan(
                degree=0,
                coeff_type=poly.coeff_type[0],
                multA=False,
                route_ops=[],
            ),
        )

    total_routes = [
        (a, b)
        for plan in term_plans
        for _, a, b in plan["route_ops"]
    ]

    comp_res.add = add_count
    return total_routes, term_plans, comp_res, final_made_powers

""" def check_without_dcmp_v2(made_powers: set[int], coeff: list[float],) -> tuple[list[tuple[int, int]] | None, Complexity | None, set[int] | None]:
    def ceil_log2(v: int) -> int:
        if v <= 1:
            return 0
        return ceil(log2(v))

    base_made_powers = {0, 1} | set(made_powers or set())

    # A. 필요한 차수 탐색
    missing_powers: list[int] = []
    for i, coef in enumerate(coeff):
        if coef != 0 and i not in base_made_powers:
            missing_powers.append(i)

    current_made_powers = set(base_made_powers)
    total_routes: list[tuple[int, int]] = []
    max_power_depth = 0

    # B. 오름차순으로 필요한 차수를 최소한의 연산으로 처리.
    for target_power in sorted(missing_powers):
        if target_power in current_made_powers:
            continue

        candidates = solve_xn_routes(False, target_power, current_made_powers)
        if not candidates:
            return None, None, None

        candidates.sort(
            key=lambda xi: (
                xi.depth,
                xi.cmult,
                xi.pmult,
                len(xi.route),
                tuple(xi.route),
            )
        )
        best_xi = candidates[0]

        # Store only effective new pure-power routes.
        for route in best_xi.route:
            a, b = route
            out_power = a + b
            if out_power not in current_made_powers:
                total_routes.append((a, b))
                current_made_powers.add(out_power)

        current_made_powers |= set(best_xi.made_powers)
        max_power_depth = max(max_power_depth, best_xi.depth)

    poly = Poly(coeff)

    comp_res = Complexity()
    val = 1 if (poly.coeff_type and poly.coeff_type[-1] == "F" and poly.deg >= 1) else 0

    # Same terminal evaluation rule as check_without_dcmp() branch in cal_polyEval.
    # eval_depth = ceil(log2(poly.deg) + val) if poly.deg > 0 else 0
    eval_depth = 0
    for i, ctype in enumerate(poly.coeff_type):
        if i == 0 or ctype == "0":
            continue
        add_pmult_depth = 1 if ctype == "F" else 0
        eval_depth = max(eval_depth, ceil_log2(i) + add_pmult_depth)
    
    comp_res.depth = max(max_power_depth, eval_depth)
    comp_res.cmult = len(total_routes)
    comp_res.pmult = poly.coeff_type[1:].count("F")
    comp_res.add = len(poly.coeff_type) - poly.coeff_type.count("0") - 1

    return total_routes, comp_res, current_made_powers
 """
# def check_without_dcmp_v2(made_powers: set[int], coeff: list[float]) -> tuple[list, Complexity]:
#     # A. 현재 필요한 x^i 검색.
#     empty_powers: set[int] = set()
#     for i, coef in enumerate(coeff):
#         if coef != 0:
#             if i not in made_powers:
#                 empty_powers.add(i)
            
#     # B. 낮은 차수부터 정렬.
#     for i in sorted(empty_powers):
#         # C. 하나씩 solve_xn_routes 수행.
        
#     # D. 최종 계산복잡도 연산.
#     # 최종 Decomp객체의 xi는 0차, multA=False로 추가, routes에 추가적으로 발생한 cmult경로를 추가하고 이에 따른 계산복잡도를 갱신한다.
#     # dcmp_p, dcmp_q는 None으로 유지한다.
#     # 
    
    
def attach(d1: Poly | XI | None, c1: Complexity, d2: Poly, c2: Complexity, attach_type: str) -> Complexity:
    '''
    Attach calculation complexity
    '''
    res = Complexity()
    if d2.coeff == []:
        return c1
    if attach_type == 'x': # Mult
        # d2 = constant
        if d2.deg == 0:
            res = c1
            add_mult = 1 if (d2.coeff_type[0] == "F" and d2.coeff[0] != 0) else 0
            res.depth += add_mult
            res.pmult += add_mult

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
        if d2.deg == 0 and d2.coeff[0] == 0.0:
            return c1
        res.depth = max(c1.depth, c2.depth)
        res.cmult = c1.cmult + c2.cmult
        res.pmult = c1.pmult + c2.pmult
        res.add = c1.add + c2.add + 1
        
    return res