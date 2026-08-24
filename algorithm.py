try:
    from .basic_class import Complexity, Poly, Decomp, XI, NULL_DECOMP
    from .poly_util import check_without_dcmp, check_without_dcmp_v2, solve_xn_routes, attach
except ImportError as e:
    from basic_class import Complexity, Poly, Decomp, XI, NULL_DECOMP
    from poly_util import check_without_dcmp, check_without_dcmp_v2, solve_xn_routes, attach
    
from math import log2, ceil, sqrt


def _ceil_log2(v: int) -> int:
    if v <= 1:
        return 0
    return ceil(log2(v))


def _best_xi(cands: list[XI], max_depth: int | None = None) -> XI | None:
    if max_depth is not None:
        cands = [xi for xi in cands if xi.depth <= max_depth]
    if not cands:
        return None

    cands.sort(
        key=lambda xi: (
            xi.depth,
            xi.cmult,
            xi.pmult,
            len(xi.route),
            tuple(xi.route),
            -len(xi.made_powers),
            tuple(sorted(xi.made_powers)),
        )
    )
    return cands[0]


def _copy_route_ops(route_ops) -> list[tuple[str, int, int]]:
    return [
        (str(kind), int(a), int(b))
        for kind, a, b in route_ops
    ]


def _make_term_plan(
    degree: int,
    coeff_type: str,
    multA: bool,
    route_ops,
) -> dict:
    return {
        "degree": int(degree),
        "coeff_type": str(coeff_type),
        "multA": bool(multA),
        "route_ops": _copy_route_ops(route_ops),
    }


def _append_constant_term_plan(poly: Poly, term_plans: list[dict]) -> None:
    if poly.coeff and poly.coeff[0] != 0:
        term_plans.insert(
            0,
            _make_term_plan(
                degree=0,
                coeff_type=poly.coeff_type[0],
                multA=False,
                route_ops=[],
            ),
        )


def _make_terminal_decomp(
    block: Poly,
    comp: Complexity,
    term_plans: list[dict],
    mp: set[int],
) -> Decomp:
    route_ops = [
        op
        for plan in term_plans
        for op in plan["route_ops"]
    ]
    routes = [(a, b) for _, a, b in route_ops]

    xi = XI(False, 0)
    xi.add_routes(
        routes,
        comp.depth,
        comp.pmult,
        comp.cmult,
        set(mp),
        route_ops=route_ops,
    )

    result = Decomp(block.coeff, comp, xi)
    result.eval_order = "terminal"
    result.term_plans = [
        {
            "degree": plan["degree"],
            "coeff_type": plan["coeff_type"],
            "multA": plan["multA"],
            "route_ops": list(plan["route_ops"]),
        }
        for plan in term_plans
    ]
    return result


def _eval_q_pure(
    poly_q: Poly,
    current_mp: set[int],
) -> tuple[Decomp, set[int]]:
    """
    q(x)는 critical path가 아니므로 각 항의 순수 x^i를 먼저 만든 뒤
    계수를 결합한다. 항별 route 종류를 term_plans에 보존한다.
    """
    mp = set(current_mp)
    term_plans: list[dict] = []
    comp = Complexity()

    for i, ctype in enumerate(poly_q.coeff_type):
        if i == 0 or ctype == "0":
            continue

        route_ops = []

        if i not in mp:
            xi = _best_xi(solve_xn_routes(False, i, mp))
            if xi is None:
                return NULL_DECOMP, mp

            route_ops = list(xi.route_ops)
            comp.cmult += xi.cmult
            mp |= set(xi.made_powers)

        add_pmult = 1 if ctype == "F" else 0
        comp.depth = max(
            comp.depth,
            _ceil_log2(i) + add_pmult,
        )
        comp.pmult += add_pmult

        term_plans.append(
            _make_term_plan(
                degree=i,
                coeff_type=ctype,
                multA=False,
                route_ops=route_ops,
            )
        )

    _append_constant_term_plan(poly_q, term_plans)
    comp.add = (
        len(poly_q.coeff_type)
        - poly_q.coeff_type.count("0")
        - 1
    )
    return _make_terminal_decomp(
        poly_q,
        comp,
        term_plans,
        mp,
    ), mp


def _eval_p_limited(
    poly_p: Poly,
    current_mp: set[int],
    max_px_depth: int,
) -> tuple[Decomp, set[int]]:
    """
    p(x)는 이후 x^k 또는 (a)x^k와 곱해지므로 각 항의 계수 포함 depth가
    max_px_depth를 넘지 않도록 구성한다. 선택된 항별 route 종류를 보존한다.
    """
    mp = set(current_mp)
    term_plans: list[dict] = []
    comp = Complexity()

    for i, ctype in enumerate(poly_p.coeff_type):
        if i == 0 or ctype == "0":
            continue

        if ctype == "I":
            route_ops = []

            if i not in mp:
                xi = _best_xi(
                    solve_xn_routes(False, i, mp),
                    max_px_depth,
                )
                if xi is None:
                    xi = _best_xi(
                        solve_xn_routes(False, i, mp)
                    )
                if xi is None:
                    return NULL_DECOMP, mp

                route_ops = list(xi.route_ops)
                comp.cmult += xi.cmult
                mp |= set(xi.made_powers)

            comp.depth = max(comp.depth, _ceil_log2(i))
            term_plans.append(
                _make_term_plan(
                    degree=i,
                    coeff_type=ctype,
                    multA=False,
                    route_ops=route_ops,
                )
            )
            continue

        use_existing_depth = _ceil_log2(i) + 1
        if i in mp and use_existing_depth <= max_px_depth:
            comp.depth = max(comp.depth, use_existing_depth)
            comp.pmult += 1
            term_plans.append(
                _make_term_plan(
                    degree=i,
                    coeff_type=ctype,
                    multA=False,
                    route_ops=[],
                )
            )
            continue

        xi = _best_xi(
            solve_xn_routes(True, i, mp),
            max_px_depth,
        )

        if xi is None:
            xi = _best_xi(solve_xn_routes(False, i, mp))
            if xi is None:
                return NULL_DECOMP, mp

            comp.cmult += xi.cmult
            comp.pmult += xi.pmult + 1
            comp.depth = max(comp.depth, xi.depth + 1)
            mp |= set(xi.made_powers)

            term_plans.append(
                _make_term_plan(
                    degree=i,
                    coeff_type=ctype,
                    multA=False,
                    route_ops=xi.route_ops,
                )
            )
        else:
            comp.cmult += xi.cmult
            comp.pmult += xi.pmult
            comp.depth = max(comp.depth, xi.depth)
            mp |= set(xi.made_powers)

            term_plans.append(
                _make_term_plan(
                    degree=i,
                    coeff_type=ctype,
                    multA=True,
                    route_ops=xi.route_ops,
                )
            )

    _append_constant_term_plan(poly_p, term_plans)
    comp.add = (
        len(poly_p.coeff_type)
        - poly_p.coeff_type.count("0")
        - 1
    )
    return _make_terminal_decomp(
        poly_p,
        comp,
        term_plans,
        mp,
    ), mp


def _eval_split_by_collected_powers(
    parent_poly: Poly,
    class_xi: XI,
    poly_p: Poly,
    poly_q: Poly,
) -> Decomp | None:
    """
    q(x)를 먼저 평가하여 생성된 power를 p(x)가 재사용하는 terminal 후보이다.
    """
    if poly_p.is_empty():
        return None

    if poly_q.is_empty():
        dcmp_q = None
        mp_after_q = set(class_xi.made_powers)
    else:
        dcmp_q, mp_after_q = _eval_q_pure(
            poly_q,
            class_xi.made_powers,
        )
        if dcmp_q is NULL_DECOMP:
            return None

    max_px_depth = (
        ceil(log2(poly_p.deg + 1))
        if poly_p.deg > 0
        else 0
    )
    dcmp_p, _ = _eval_p_limited(
        poly_p,
        mp_after_q,
        max_px_depth,
    )
    if dcmp_p is NULL_DECOMP:
        return None

    comp_i = Complexity(class_xi)
    comp_xp = attach(
        class_xi,
        comp_i,
        poly_p,
        dcmp_p.comp,
        'x',
    )
    comp_total = (
        comp_xp
        if dcmp_q is None
        else attach(None, comp_xp, poly_q, dcmp_q.comp, '+')
    )

    result = Decomp(parent_poly.coeff, comp_total)
    result.update(class_xi, dcmp_p, dcmp_q)
    result.eval_order = "q_then_p" if dcmp_q is not None else "p_only"
    return result

def cal_polyEval(poly: Poly | list[float], made_powers: set[int] | None = None, is_root: bool = True, consider_axi: bool = False) -> Decomp:
    if made_powers is None:
        made_powers = {0, 1}
    else:
        made_powers = set(made_powers)
        
    if type(poly) != Poly:
        poly = Poly(poly)
    results: list[Decomp] = []
    
    # 1. Check evaluation availability without decomposition. = All x^i is memoized.
    if check_without_dcmp(made_powers, poly.coeff):
        comp_res = Complexity()
        val = 1 if (poly.coeff_type[-1] == "F" and poly.deg >= 1) else 0
        comp_res.depth = ceil(log2(poly.deg) + val) if poly.deg > 0 else 0
        comp_res.cmult = 0
        comp_res.pmult = poly.coeff_type[1:].count("F")
        comp_res.add = len(poly.coeff_type) - poly.coeff_type.count("0") - 1

        term_plans = [
            _make_term_plan(i, ctype, False, [])
            for i, ctype in enumerate(poly.coeff_type)
            if ctype != "0"
        ]
        return _make_terminal_decomp(
            poly,
            comp_res,
            term_plans,
            set(made_powers),
        )
    
    # 2. Check evaluation availability without decomposition when not all x^i is memoized.
    # Rule: max_degree x^k must be made.
    if not is_root:
        new_routes, new_term_plans, new_comp, new_made_powers = check_without_dcmp_v2(
            made_powers,
            poly.coeff,
        )

        if (
            new_routes is not None
            and new_term_plans is not None
            and new_comp is not None
            and new_made_powers is not None
        ):
            val = 1 if poly.coeff_type[poly.deg] == "F" else 0
            target_depth = ceil(log2(poly.deg + val)) if poly.deg > 0 else 0
            target_cmult = poly.deg - 1 if poly.deg > 0 else 0

            if new_comp.depth <= target_depth and new_comp.cmult <= target_cmult:
                results.append(
                    _make_terminal_decomp(
                        poly,
                        new_comp,
                        new_term_plans,
                        new_made_powers,
                    )
                )
            
    # 3. Calculate complexity with decomposition.
    # best_i = [i for i in range(1, poly.deg+1) if poly.coeff_type[i] != '0']
    for xi in range(1, poly.deg+1):
        # Case 1: multA=True
        if poly.coeff_type[poly.deg] == "F" and consider_axi:
            XIs = solve_xn_routes(True, xi, made_powers)
            for class_xi in XIs:
                res = process_recursion(poly, class_xi, consider_axi)
                if res is not None:
                    results.extend(res)

        # Case 2: multA=False
        XIs = solve_xn_routes(False, xi, made_powers)
        for class_xi in XIs:
            res = process_recursion(poly, class_xi, consider_axi)
            if res is not None:
                results.extend(res)

    # Sort and return
    if len(results) == 0:
        return NULL_DECOMP
    
    results.sort()
    return results[0]


def process_recursion(poly: Poly, xi: XI, consider_axi: bool=True) -> list[Decomp] | None:
    max_deg = poly.deg
    val = 1 if poly.coeff_type[max_deg] == "F" else 0
    
    # Depth, CMult values are fixed by the degree of polynomial.
    # If the decomposition's complexity exceeds the target, we skip it.
    target_depth = ceil(log2(max_deg + val))
    target_cmult = max_deg - 1
    # target_cmult = ceil(sqrt(max_deg))+1
    
    poly_p, poly_q = poly.seperate(xi.n, xi.multA)
    comp_i = Complexity(xi)
    
    # Condition 1: poly_p, poly_q are all empty.
    if poly_p.is_empty() and poly_q.is_empty():
        if comp_i.depth <= target_depth and comp_i.cmult <= target_cmult:
            result = Decomp(poly.coeff, comp_i)
            result.update(xi, None, None)
            return [result]        
        return None
    
    # Condition 2: poly_q is empty
    if not poly_p.is_empty() and poly_q.is_empty():
        decomp_p: Decomp = cal_polyEval(poly_p, xi.made_powers, is_root=False, consider_axi=consider_axi)
        if decomp_p is NULL_DECOMP:
            return None
        comp_pi = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
        if comp_pi.depth <= target_depth and comp_pi.cmult <= target_cmult:
            result = Decomp(poly.coeff, comp_pi)
            result.update(xi, decomp_p, None)
            result.eval_order = "p_only"
            return [result]
        return None
    
    # Condition 3: poly_p is empty
    if poly_p.is_empty() and not poly_q.is_empty():
        decomp_q: Decomp = cal_polyEval(poly_q, xi.made_powers, is_root=False, consider_axi=consider_axi)
        if decomp_q is NULL_DECOMP:
            return None
        comp_iq = attach(xi, comp_i, poly_q, decomp_q.comp, '+')
        if comp_iq.depth <= target_depth and comp_iq.cmult <= target_cmult:
            result = Decomp(poly.coeff, comp_iq)
            result.update(xi, None, decomp_q)
            result.eval_order = "q_only"
            return [result]
        return None

    # Condition 4: all poly is not empty
    results = []
    if not poly_p.is_empty() and not poly_q.is_empty():
        # 0. p(x), q(x)에 필요한 차수를 모아 terminal 방식으로 한 번에 평가하는 후보를 추가한다.
        #    기존 p->q / q->p 재귀 후보와 경쟁시킨다.
        collected_result = _eval_split_by_collected_powers(poly, xi, poly_p, poly_q)
        if collected_result is not None:
            if collected_result.comp.depth == target_depth and collected_result.comp.cmult <= target_cmult:
                results.append(collected_result)

        # dcmp_p -> dcmp_q
        decomp_p = cal_polyEval(poly_p, xi.made_powers, is_root=False, consider_axi=consider_axi)
        decomp_q = cal_polyEval(poly_q, decomp_p.merge_mp() | xi.made_powers, is_root=False, consider_axi=consider_axi)
        
        if decomp_p is not NULL_DECOMP and decomp_q is not NULL_DECOMP:
            comp_pi = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
            comp_piq = attach(None, comp_pi, poly_q, decomp_q.comp, '+')

            if comp_piq.depth == target_depth and comp_piq.cmult <= target_cmult:
                result = Decomp(poly.coeff, comp_piq)
                result.update(xi, decomp_p, decomp_q)
                result.eval_order = "p_then_q"
                results.append(result)
                
        # dcmp_q -> dcmp_p
        decomp_q = cal_polyEval(poly_q, xi.made_powers, is_root=False, consider_axi=consider_axi)
        decomp_p = cal_polyEval(poly_p, decomp_q.merge_mp() | xi.made_powers, is_root=False, consider_axi=consider_axi)
        
        if decomp_p is not NULL_DECOMP and decomp_q is not NULL_DECOMP:
            comp_pi = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
            comp_piq = attach(None, comp_pi, poly_q, decomp_q.comp, '+')

            if comp_piq.depth == target_depth and comp_piq.cmult <= target_cmult:
                result = Decomp(poly.coeff, comp_piq)
                result.update(xi, decomp_p, decomp_q)
                result.eval_order = "q_then_p"
                results.append(result)
        
        if len(results) == 0:
            return None
        else:
            return results

    
# def cal_PSMethod(poly: Poly | list[float], made_powers: set[int] | None = None, is_root: bool = True) -> Decomp:
#     if made_powers is None:
#         made_powers = {0, 1}
#     else:
#         made_powers = set(made_powers)
        
#     if type(poly) != Poly:
#         poly = Poly(poly)   
        
#     # Calculate complexity with decomposition. k = ceil(n/2)
#     k = ceil(poly.deg/2)
#     XIs = solve_xn_routes(False, k, made_powers)
#     poly_p, poly_q = poly.seperate(k, False)
#     assert poly_p.deg >= poly_q.deg
    
#     # q(x) - a * x^i로 비효율적이게 구성해도 depth가 무조건 p(x) * x^k보다 작음.
#     # 따라서 q(x)의 모든 항은 solve_xn_routes에서 multA=False로 두고 생성.
#     # q(x)내부의 모든 원소(x^i들)를 최대한 효율적이게 만들고, made_powers 업데이트.
    
#     max_px_depth = ceil(log2(poly_p.deg+1)) # p(x) 내부 원소가 가질 수 있는 최대 깊이.
#     # p(x) - 작은 원소부터 판별.
#     # 경우의 수를 분리
#     # i가 made_powers에 존재, log2(i+1)이 max_px_depth 이하 -> pmult +1로 생성.
#     # i가 made_powers에 존재하지만 log2(i+1)이 max_px_depth 초과 -> solve_xn_routes로 max_px_depth를 넘지 않게 구성.
#     # i가 made_powers에 존재하지 않음 -> solve_xn_routes로 최대한 효율적으로 구성.
    
#     # 최종 계산복잡도 연산
#     return result


def cal_PSMethod(poly: Poly | list[float], made_powers: set[int] | None = None, is_root: bool = True) -> Decomp:
    if made_powers is None:
        made_powers = {0, 1}
    else:
        made_powers = set(made_powers)

    if type(poly) != Poly:
        poly = Poly(poly)

    if poly.deg == 0:
        return Decomp(poly.coeff, Complexity())

    # k = ceil(deg/2)로 한 번만 분해한다.
    k = ceil(poly.deg / 2)
    XIs = solve_xn_routes(False, k, made_powers)

    if not XIs:
        return NULL_DECOMP

    results: list[Decomp] = []

    for class_xi in XIs:
        poly_p, poly_q = poly.seperate(k, False)

        # PS split에서는 p(x)의 차수가 q(x)의 차수보다 낮지 않아야 한다.
        if poly_p.deg < poly_q.deg:
            continue

        result = _eval_split_by_collected_powers(poly, class_xi, poly_p, poly_q)
        if result is not None:
            results.append(result)

    if len(results) == 0:
        return NULL_DECOMP

    results.sort()
    return results[0]
