from .basic_class import Complexity, Poly, Decomp, XI, NULL_DECOMP
from .poly_util import check_without_dcmp, check_without_dcmp_v2_all, solve_xn_routes, solve_xn_routes_depth_limited, attach

from math import log2, ceil, sqrt


def _normalize_power_depths(power_cache) -> dict[int, int]:
    '''
    기능:
        set 또는 power-depth dictionary를 OPD가 사용하는 실제 depth 상태로 변환한다.
    입력:
        power_cache: 순수 power 집합, {power: depth} dictionary 또는 None.
    출력:
        x^0, x^1을 포함하는 power-depth dictionary.
    '''
    if isinstance(power_cache, dict):
        result = {
            int(power): int(power_depth)
            for power, power_depth in power_cache.items()
        }
    else:
        result = {
            int(power): ceil(log2(power)) if power > 1 else 0
            for power in set(power_cache or set())
        }
    result.setdefault(0, 0)
    result.setdefault(1, 0)
    return result


def _merge_power_depths(*power_caches) -> dict[int, int]:
    '''
    기능:
        여러 power-depth cache를 병합하고 같은 power에는 가장 작은 depth를 유지한다.
    입력:
        power_caches: 병합할 power-depth dictionary 목록.
    출력:
        병합된 power-depth dictionary.
    '''
    result = {0: 0, 1: 0}
    for power_cache in power_caches:
        for power, power_depth in power_cache.items():
            result[power] = min(
                result.get(power, power_depth),
                power_depth,
            )
    return result


def _power_state_key(power_depths: dict[int, int]) -> frozenset:
    '''
    기능:
        power별 실제 depth 상태를 cache 비교용 hashable 값으로 변환한다.
    입력:
        power_depths: {power: actual depth} dictionary.
    출력:
        (power, depth) 쌍으로 구성된 frozenset.
    '''
    return frozenset(
        (int(power), int(power_depth))
        for power, power_depth in power_depths.items()
    )


def _make_cache_key(poly: Poly, consider_axi: bool, assumption1: bool, depth_budget: int) -> tuple:
    '''
    기능:
        OPD 탐색상태의 다항식 식별자와 탐색옵션을 캐시 키로 변환한다.
    입력:
        poly: 현재 탐색할 Poly 객체.
        consider_axi: ax^i 후보 탐색 여부.
        assumption1: Assumption 1 적용 여부.
        depth_budget: 현재 탐색에 허용되는 최대 multiplicative depth.
    출력:
        타입 입력은 계수타입, 수치 입력은 실제 계수와 탐색옵션으로 구성된 tuple.
    '''
    poly_key = (
        tuple(poly.coeff_type)
        if poly.is_type_poly
        else tuple(poly.coeff)
    )

    return (
        poly_key,
        bool(poly.is_type_poly),
        bool(consider_axi),
        bool(assumption1),
        int(depth_budget),
    )


def _make_decomp_signature(dcmp: Decomp | None) -> tuple | None:
    '''
    기능:
        중복 캐시 후보를 판별할 수 있도록 decomposition 트리를 tuple로 변환한다.
    입력:
        dcmp: signature를 생성할 Decomp 객체 또는 None.
    출력:
        계수타입, factor, 평가순서, 자식구조를 포함한 tuple 또는 None.
    '''
    if dcmp is None:
        return None

    return (
        tuple(dcmp.coeff_type),
        bool(dcmp.xi.multA),
        int(dcmp.xi.n),
        tuple(dcmp.xi.route_ops),
        tuple(sorted(dcmp.xi.power_depths.items())),
        str(dcmp.eval_order),
        tuple(
            (
                int(plan["degree"]),
                str(plan["coeff_type"]),
                bool(plan["multA"]),
                tuple(plan["route_ops"]),
            )
            for plan in dcmp.term_plans
        ),
        _make_decomp_signature(dcmp.dcmp_p),
        _make_decomp_signature(dcmp.dcmp_q),
    )


def _find_cached_decomps(cache: dict, key: tuple, power_depths: dict[int, int]) -> list[Decomp]:
    '''
    기능:
        동일한 다항식 타입과 입력 power cache에 저장된 최적 decomposition을 찾는다.
    입력:
        cache: SearchResult tuple 목록을 저장한 dictionary.
        key: 다항식 타입과 탐색옵션으로 구성된 캐시 키.
        power_depths: 현재 입력 pure power의 실제 depth cache.
    출력:
        일치하는 모든 Decomp 후보 목록.
    '''
    power_key = _power_state_key(power_depths)
    candidates = [
        dcmp for powers, dcmp, _ in cache.get(key, [])
        if powers == power_key
    ]

    candidates.sort()
    return candidates


# Store the decomposition candidates in the cache as (power cache, decomposition, complexity) tuples.
def _store_cached_decomps(cache: dict, key: tuple, power_depths: dict[int, int], results: list[Decomp]) -> None:
    '''
    기능:
        탐색 후보들을 (power cache, decomposition, complexity) tuple로 저장한다.
    입력:
        cache: SearchResult tuple 목록을 저장할 dictionary.
        key: 다항식 타입과 탐색옵션으로 구성된 캐시 키.
        power_depths: 후보 탐색에 사용된 입력 pure power의 실제 depth cache.
        results: 저장할 Decomp 후보 목록.
    출력:
        없음. cache dictionary를 직접 갱신한다.
    '''
    power_key = _power_state_key(power_depths)
    stored_results = cache.setdefault(key, [])
    stored_signatures = {
        (powers, _make_decomp_signature(dcmp), comp.return_params())
        for powers, dcmp, comp in stored_results
    }

    for result in results:
        signature = (
            power_key,
            _make_decomp_signature(result),
            result.comp.return_params(),
        )

        if signature in stored_signatures:
            continue

        stored_results.append((power_key, result, result.comp))
        stored_signatures.add(signature)


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


def _make_term_plan(degree: int, coeff_type: str, multA: bool, route_ops) -> dict:
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
    power_depths: set[int] | dict[int, int],
) -> Decomp:
    '''
    기능:
        항별 terminal 평가계획과 출력 power-depth cache를 leaf Decomp로 구성한다.
    입력:
        block: terminal로 평가한 다항식.
        comp: terminal 평가의 계산복잡도.
        term_plans: 항별 power 및 계수 결합 계획.
        power_depths: 평가 후 사용 가능한 순수 power의 실제 depth cache.
    출력:
        terminal 평가계획을 저장한 Decomp 객체.
    '''
    power_depths = _normalize_power_depths(power_depths)
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
        set(power_depths),
        route_ops=route_ops,
        power_depths=power_depths,
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


def _eval_split_by_collected_powers(parent_poly: Poly, class_xi: XI, poly_p: Poly, poly_q: Poly, eval_order: str = "q_then_p",) -> Decomp | None:
    '''
    기능:
        p(x), q(x) 중 지정된 다항식을 먼저 terminal 방식으로 평가하고,
        이때 생성된 순수 power를 나중 다항식이 재사용하는 후보를 만든다.
        XI와 곱해지는 p(x)는 평가순서와 관계없이 제한된 깊이로 평가한다.
    입력:
        parent_poly: 분리 전 부모 다항식.
        class_xi: 부모와 quotient를 결합할 x^i 또는 ax^i factor 정보.
        poly_p: 현재 분리의 quotient 다항식.
        poly_q: 현재 분리의 remainder 다항식.
        eval_order: "q_then_p" 또는 "p_then_q" terminal 평가순서.
    출력:
        먼저 생성된 순수 power를 재사용하여 구성한 Decomp 객체.
        입력이 유효하지 않거나 평가할 수 없으면 None.
    '''
    if poly_p.is_empty():
        return None

    if eval_order not in ("q_then_p", "p_then_q"):
        return None

    max_px_depth = (
        ceil(log2(poly_p.deg + 1))
        if poly_p.deg > 0
        else 0
    )

    if eval_order == "q_then_p":
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

        dcmp_p, _ = _eval_p_limited(
            poly_p,
            mp_after_q,
            max_px_depth,
        )
        if dcmp_p is NULL_DECOMP:
            return None
    else:
        dcmp_p, mp_after_p = _eval_p_limited(
            poly_p,
            class_xi.made_powers,
            max_px_depth,
        )
        if dcmp_p is NULL_DECOMP:
            return None

        if poly_q.is_empty():
            dcmp_q = None
        else:
            dcmp_q, _ = _eval_q_pure(
                poly_q,
                mp_after_p,
            )
            if dcmp_q is NULL_DECOMP:
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
    result.eval_order = eval_order if dcmp_q is not None else "p_only"
    return result

'''
cal_polyEval
 ├─ log-depth budget과 입력 power별 실제 depth로 탐색상태 구성
 ├─ 최종 budget을 만족하는 모든 terminal 계획 추가
 ├─ i = 1, ..., deg에 대해
 │   ├─ 중간 power의 최소 depth를 강제하지 않는 모든 x^i 생성 DAG
 │   ├─ 최고차항이 F이면 동일 조건의 모든 ax^i 생성 DAG
 │   └─ f = (x^i 또는 ax^i)p + q로 분해
 │       ├─ p → power cache 확장 → q의 모든 후보 조합
 │       └─ q → power cache 확장 → p의 모든 후보 조합
 ├─ CMult 상한 없이 최종 depth 이하 후보만 유지
 ├─ 구조가 정확히 같은 후보만 제거하고 전체를 search cache에 저장
 └─ cal_polyEval wrapper가 CKKS complexity 최소 후보 반환
'''
def _get_depth_limit(poly: Poly) -> int:
    '''
    기능:
        다항식 차수와 최고차항 타입으로 log-depth 상한을 계산한다.
    입력:
        poly: depth 상한을 계산할 Poly 객체.
    출력:
        최고차항이 F이면 한 단계를 반영한 최대 multiplicative depth.
    '''
    if poly.deg <= 0:
        return 0
    val = 1 if poly.coeff_type[poly.deg] == "F" else 0
    return ceil(log2(poly.deg + val))


def _dedupe_decomps(results: list[Decomp]) -> list[Decomp]:
    '''
    기능:
        구조와 complexity가 완전히 같은 decomposition 후보만 제거한다.
    입력:
        results: 중복을 포함할 수 있는 Decomp 후보 목록.
    출력:
        최초 등장 순서를 유지한 고유 Decomp 후보 목록.
    '''
    unique_results = []
    seen = set()
    for result in results:
        signature = (
            _make_decomp_signature(result),
            result.comp.return_params(),
            tuple(sorted(result.merge_power_depths().items())),
        )
        if signature in seen:
            continue
        seen.add(signature)
        unique_results.append(result)
    return unique_results


def cal_polyEval_candidates(poly: Poly | list[float] | list[str], made_powers: set[int] | dict[int, int] | None = None, depth_budget: int | None = None, is_root: bool = True, consider_axi: bool = True, assumption1: bool = True, initial_cache: dict | None = None, intermediate_cache: dict | None = None) -> list[Decomp]:
    '''
    기능:
        중간 power의 개별 최소 깊이를 강제하지 않고 최종 다항식 depth 상한을
        만족하는 모든 비중복 decomposition 및 하위 power 재사용 조합을 반환한다.
    입력:
        poly: 수치 계수, 0/I/F 계수타입 목록 또는 Poly 객체.
        made_powers: 현재 순수 x^i 집합 또는 각 power의 실제 depth dictionary.
        depth_budget: 현재 다항식 평가에 허용되는 최대 multiplicative depth.
        is_root: 초기 다항식 탐색 여부.
        consider_axi: True이면 ax^i factor 후보를 탐색한다.
        assumption1: ax^i 정규화에 Assumption 1을 적용할지 여부.
        initial_cache: 입력 power cache가 {0,1}인 후보 cache.
        intermediate_cache: 그 외 입력 power cache의 후보 cache.
    출력:
        depth_budget 이하의 모든 고유 Decomp 후보 목록.
    '''
    power_depths = _normalize_power_depths(made_powers)
    made_powers = set(power_depths)

    if not isinstance(poly, Poly):
        poly = Poly(poly)

    if depth_budget is None:
        depth_budget = _get_depth_limit(poly)

    if initial_cache is None:
        initial_cache = {}
    if intermediate_cache is None:
        intermediate_cache = {}

    cache = (
        initial_cache
        if power_depths == {0: 0, 1: 0}
        else intermediate_cache
    )
    cache_key = _make_cache_key(
        poly,
        consider_axi,
        assumption1,
        depth_budget,
    )
    cached_results = _find_cached_decomps(
        cache,
        cache_key,
        power_depths,
    )
    if len(cached_results) > 0:
        return cached_results

    results = []

    if check_without_dcmp(made_powers, poly.coeff):
        comp_res = Complexity()
        comp_res.depth = max(
            (
                power_depths[degree]
                + (1 if coeff_type == "F" else 0)
            )
            for degree, coeff_type in enumerate(poly.coeff_type)
            if degree >= 1 and coeff_type != "0"
        ) if poly.deg > 0 else 0
        comp_res.cmult = 0
        comp_res.pmult = poly.coeff_type[1:].count("F")
        comp_res.add = max(
            len(poly.coeff_type) - poly.coeff_type.count("0") - 1,
            0,
        )
        if comp_res.depth <= depth_budget:
            term_plans = [
                _make_term_plan(i, ctype, False, [])
                for i, ctype in enumerate(poly.coeff_type)
                if ctype != "0"
            ]
            results.append(
                _make_terminal_decomp(
                    poly,
                    comp_res,
                    term_plans,
                    power_depths,
                )
            )
    else:
        for _, term_plans, comp_res, new_power_depths in check_without_dcmp_v2_all(
            power_depths,
            poly.coeff_type if poly.is_type_poly else poly.coeff,
            depth_budget,
        ):
            results.append(
                _make_terminal_decomp(
                    poly,
                    comp_res,
                    term_plans,
                    new_power_depths,
                )
            )

    for factor_power in range(1, poly.deg):
        factor_types = [False]
        if consider_axi and poly.coeff_type[poly.deg] == "F":
            factor_types.insert(0, True)

        for multA in factor_types:
            factor_routes = solve_xn_routes_depth_limited(
                multA,
                factor_power,
                power_depths,
                depth_budget,
            )
            for class_xi in factor_routes:
                case_results = process_recursion(
                    poly,
                    class_xi,
                    depth_budget,
                    consider_axi,
                    assumption1,
                    initial_cache,
                    intermediate_cache,
                )
                results.extend(case_results)

    results = _dedupe_decomps(results)
    results.sort()
    if len(results) > 0:
        _store_cached_decomps(
            cache,
            cache_key,
            power_depths,
            results,
        )
    return results


def cal_polyEval(poly: Poly | list[float] | list[str], made_powers: set[int] | dict[int, int] | None = None, is_root: bool = True, consider_axi: bool = True, assumption1: bool = True, initial_cache: dict | None = None, intermediate_cache: dict | None = None, depth_budget: int | None = None) -> Decomp:
    '''
    기능:
        depth 제한을 만족하는 모든 decomposition을 탐색한 뒤 CKKS complexity가
        가장 작은 후보를 반환한다.
    입력:
        poly: 수치 계수, 0/I/F 계수타입 목록 또는 Poly 객체.
        made_powers: 현재 순수 x^i 집합 또는 각 power의 실제 depth dictionary.
        is_root: 초기 다항식 탐색 여부.
        consider_axi: True이면 ax^i factor 후보를 탐색한다.
        assumption1: ax^i 정규화에 Assumption 1을 적용할지 여부.
        initial_cache: 입력 power cache가 {0,1}인 후보 cache.
        intermediate_cache: 그 외 입력 power cache의 후보 cache.
        depth_budget: 허용할 최대 multiplicative depth.
    출력:
        전체 depth-feasible 후보 중 최소 complexity Decomp 또는 NULL_DECOMP.
    '''
    candidates = cal_polyEval_candidates(
        poly,
        made_powers,
        depth_budget,
        is_root,
        consider_axi,
        assumption1,
        initial_cache,
        intermediate_cache,
    )
    return candidates[0] if len(candidates) > 0 else NULL_DECOMP


def process_recursion(poly: Poly, xi: XI, depth_budget: int, consider_axi: bool=True, assumption1: bool=True, initial_cache: dict | None = None, intermediate_cache: dict | None = None) -> list[Decomp]:
    '''
    기능:
        하나의 x^i 또는 ax^i factor에 대한 모든 quotient 타입과 depth-feasible
        하위 decomposition 조합을 탐색한다.
    입력:
        poly: 분리 전 부모 다항식.
        xi: 현재 factor의 power 생성정보.
        depth_budget: 부모 다항식에 허용되는 최대 multiplicative depth.
        consider_axi: 하위 재귀에서 ax^i 후보를 탐색할지 여부.
        assumption1: ax^i 정규화에 적용할 Assumption 1 여부.
        initial_cache: 입력 power cache가 {0,1}인 후보 cache.
        intermediate_cache: 그 외 입력 power cache의 후보 cache.
    출력:
        현재 factor에서 생성된 모든 고유 Decomp 후보 목록.
    '''
    results = []
    for poly_p, poly_q in poly.seperate_cases(
        xi.n,
        xi.multA,
        assumption1,
    ):
        results.extend(
            _process_recursion_case(
                poly,
                xi,
                poly_p,
                poly_q,
                depth_budget,
                consider_axi,
                assumption1,
                initial_cache,
                intermediate_cache,
            )
        )
    return _dedupe_decomps(results)


def _process_recursion_case(poly: Poly, xi: XI, poly_p: Poly, poly_q: Poly, depth_budget: int, consider_axi: bool=True, assumption1: bool=True, initial_cache: dict | None = None, intermediate_cache: dict | None = None) -> list[Decomp]:
    '''
    기능:
        확정된 quotient/remainder 타입에서 p->q 및 q->p의 모든 후보 조합을
        만들고 depth_budget을 만족하는 결과만 반환한다.
    입력:
        poly: 분리 전 부모 다항식.
        xi: 부모와 quotient를 결합할 x^i 또는 ax^i factor.
        poly_p: quotient 다항식.
        poly_q: remainder 다항식.
        depth_budget: 부모 다항식에 허용되는 최대 multiplicative depth.
        consider_axi: 하위 재귀에서 ax^i 후보를 탐색할지 여부.
        assumption1: 하위 ax^i 정규화에 적용할 Assumption 1 여부.
        initial_cache: 입력 power cache가 {0,1}인 후보 cache.
        intermediate_cache: 그 외 입력 power cache의 후보 cache.
    출력:
        depth 조건을 만족하는 모든 고유 Decomp 후보 목록.
    '''
    comp_i = Complexity(xi)
    if comp_i.depth > depth_budget:
        return []

    if poly_p.is_empty() and poly_q.is_empty():
        result = Decomp(poly.coeff, comp_i)
        result.update(xi, None, None)
        return [result]

    if not poly_p.is_empty() and poly_q.is_empty():
        results = []
        for decomp_p in cal_polyEval_candidates(
            poly_p,
            xi.power_depths,
            depth_budget,
            False,
            consider_axi,
            assumption1,
            initial_cache,
            intermediate_cache,
        ):
            comp_total = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
            if comp_total.depth <= depth_budget:
                result = Decomp(poly.coeff, comp_total)
                result.update(xi, decomp_p, None)
                result.eval_order = "p_only"
                results.append(result)
        return _dedupe_decomps(results)

    if poly_p.is_empty() and not poly_q.is_empty():
        results = []
        for decomp_q in cal_polyEval_candidates(
            poly_q,
            xi.power_depths,
            depth_budget,
            False,
            consider_axi,
            assumption1,
            initial_cache,
            intermediate_cache,
        ):
            comp_total = attach(xi, comp_i, poly_q, decomp_q.comp, '+')
            if comp_total.depth <= depth_budget:
                result = Decomp(poly.coeff, comp_total)
                result.update(xi, None, decomp_q)
                result.eval_order = "q_only"
                results.append(result)
        return _dedupe_decomps(results)

    results = []

    for decomp_p in cal_polyEval_candidates(
        poly_p,
        xi.power_depths,
        depth_budget,
        False,
        consider_axi,
        assumption1,
        initial_cache,
        intermediate_cache,
    ):
        powers_after_p = _merge_power_depths(
            xi.power_depths,
            decomp_p.merge_power_depths(),
        )
        for decomp_q in cal_polyEval_candidates(
            poly_q,
            powers_after_p,
            depth_budget,
            False,
            consider_axi,
            assumption1,
            initial_cache,
            intermediate_cache,
        ):
            comp_xp = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
            comp_total = attach(None, comp_xp, poly_q, decomp_q.comp, '+')
            if comp_total.depth <= depth_budget:
                result = Decomp(poly.coeff, comp_total)
                result.update(xi, decomp_p, decomp_q)
                result.eval_order = "p_then_q"
                results.append(result)

    for decomp_q in cal_polyEval_candidates(
        poly_q,
        xi.power_depths,
        depth_budget,
        False,
        consider_axi,
        assumption1,
        initial_cache,
        intermediate_cache,
    ):
        powers_after_q = _merge_power_depths(
            xi.power_depths,
            decomp_q.merge_power_depths(),
        )
        for decomp_p in cal_polyEval_candidates(
            poly_p,
            powers_after_q,
            depth_budget,
            False,
            consider_axi,
            assumption1,
            initial_cache,
            intermediate_cache,
        ):
            comp_xp = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
            comp_total = attach(None, comp_xp, poly_q, decomp_q.comp, '+')
            if comp_total.depth <= depth_budget:
                result = Decomp(poly.coeff, comp_total)
                result.update(xi, decomp_p, decomp_q)
                result.eval_order = "q_then_p"
                results.append(result)

    return _dedupe_decomps(results)


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
