# save.py
import json, os
from math import ceil, log2
from .basic_class import Poly, Decomp, XI, Complexity, NULL_DECOMP
from .algorithm import cal_polyEval

# 다항식 계수타입을 키로 반환
def get_poly_type_key(poly_type: list, assumption1: bool=True) -> str:
    '''
    기능:
        다항식 계수타입을 캐시 키로 변환한다.
    입력:
        poly_type: 수치 계수 목록 또는 0/I/F 계수타입 목록.
        assumption1: 기존 호출부와의 호환성을 위해 유지하는 Assumption 1 설정.
    출력:
        0/I/F 계수타입만 연결한 문자열 캐시 키.
    '''
    if type(poly_type[0]) is not str:
        poly_type = Poly(poly_type).coeff_type
    return ''.join(poly_type)


def serialize_route_ops(route_ops) -> list[dict]:
    return [
        {
            "kind": str(kind),
            "lhs": int(lhs),
            "rhs": int(rhs),
        }
        for kind, lhs, rhs in route_ops
    ]


def deserialize_route_ops(raw_ops) -> list[tuple[str, int, int]]:
    result = []
    for index, op in enumerate(raw_ops or []):
        if not isinstance(op, dict):
            raise ValueError(
                f"route_ops[{index}] must be an object"
            )
        kind = op.get("kind")
        if kind not in {"pure", "coeff", "coeff_direct"}:
            raise ValueError(
                f"route_ops[{index}] has invalid kind: {kind}"
            )
        result.append(
            (
                kind,
                int(op.get("lhs", 0)),
                int(op.get("rhs", 0)),
            )
        )
    return result


def serialize_term_plan(plan: dict) -> dict:
    return {
        "degree": int(plan["degree"]),
        "coeff_type": str(plan["coeff_type"]),
        "multA": bool(plan["multA"]),
        "route_ops": serialize_route_ops(plan["route_ops"]),
    }


def deserialize_term_plan(plan: dict) -> dict:
    if not isinstance(plan, dict):
        raise ValueError("term plan must be an object")
    return {
        "degree": int(plan["degree"]),
        "coeff_type": str(plan["coeff_type"]),
        "multA": bool(plan.get("multA", False)),
        "route_ops": deserialize_route_ops(
            plan.get("route_ops", [])
        ),
    }


def serialize_complexity(comp: Complexity) -> dict:
    return {
        "depth": int(comp.depth),
        "cmult": int(comp.cmult),
        "pmult": int(comp.pmult),
        "add": int(comp.add),
    }


def deserialize_complexity(raw: dict | None) -> Complexity:
    comp = Complexity()
    if not raw:
        return comp
    comp.insert_value(
        int(raw.get("depth", 0)),
        int(raw.get("cmult", 0)),
        int(raw.get("pmult", 0)),
        int(raw.get("add", 0)),
    )
    return comp


def normalize_power_depths(power_cache) -> dict[int, int]:
    '''
    기능:
        set, (power, depth) 쌍 또는 dictionary를 실제 power-depth 상태로 변환한다.
    입력:
        power_cache: 순수 power cache 표현.
    출력:
        x^0, x^1을 포함하는 {power: actual depth} dictionary.
    '''
    if isinstance(power_cache, dict):
        result = {
            int(power): int(power_depth)
            for power, power_depth in power_cache.items()
        }
    else:
        values = list(power_cache or [])
        if values and isinstance(values[0], (list, tuple)):
            result = {}
            for power, power_depth in values:
                power = int(power)
                power_depth = int(power_depth)
                result[power] = min(
                    result.get(power, power_depth),
                    power_depth,
                )
        else:
            result = {
                int(power): ceil(log2(power)) if power > 1 else 0
                for power in values
            }
    result.setdefault(0, 0)
    result.setdefault(1, 0)
    return result


def serialize_power_depths(power_depths: dict[int, int]) -> list[list[int]]:
    '''
    기능:
        power-depth dictionary를 JSON 직렬화 가능한 정렬된 pair 목록으로 변환한다.
    입력:
        power_depths: {power: actual depth} dictionary.
    출력:
        [[power, depth], ...] 형식의 정렬된 목록.
    '''
    return [
        [int(power), int(power_depth)]
        for power, power_depth in sorted(power_depths.items())
    ]


def serialize_xi(xi: XI | None) -> dict | None:
    if xi is None:
        return None

    return {
        "n": int(xi.n),
        "multA": bool(xi.multA),
        # 구버전 reader 확인용으로 유지한다.
        "route": [
            [int(lhs), int(rhs)]
            for lhs, rhs in xi.route
        ],
        "route_ops": serialize_route_ops(xi.route_ops),
        "power_depths": serialize_power_depths(xi.power_depths),
        "depth": int(xi.depth),
        "cmult": int(xi.cmult),
        "pmult": int(xi.pmult),
    }


def deserialize_xi(
    raw: dict | None,
    incoming_powers: set[int] | dict[int, int],
) -> XI:
    if raw is None:
        xi = XI()
        xi.add_routes(
            [],
            0,
            0,
            0,
            set(normalize_power_depths(incoming_powers)),
            route_ops=[],
            power_depths=normalize_power_depths(incoming_powers),
        )
        return xi

    route_ops = deserialize_route_ops(raw.get("route_ops", []))

    # 구형 cache 호환: route_ops가 없으면 기존 route를 pure로 간주한다.
    if not route_ops:
        route_ops = [
            ("pure", int(lhs), int(rhs))
            for lhs, rhs in raw.get("route", [])
        ]

    routes = [(lhs, rhs) for _, lhs, rhs in route_ops]
    available_depths = normalize_power_depths(incoming_powers)

    available_coeff = {1}
    if 2 in available_depths:
        available_coeff.add(2)

    for index, (kind, lhs, rhs) in enumerate(route_ops):
        if kind == "pure":
            if lhs not in available_depths or rhs not in available_depths:
                raise ValueError(
                    f"route_ops[{index}] uses unavailable pure power"
                )
            target = lhs + rhs
            generated_depth = max(
                available_depths[lhs],
                available_depths[rhs],
            ) + 1
            available_depths[target] = min(
                available_depths.get(target, generated_depth),
                generated_depth,
            )
        elif kind == "coeff":
            if lhs not in available_coeff or rhs not in available_depths:
                raise ValueError(
                    f"route_ops[{index}] uses unavailable coefficient power"
                )
            available_coeff.add(lhs + rhs)
        elif kind == "coeff_direct":
            if rhs not in available_depths:
                raise ValueError(
                    f"route_ops[{index}] uses unavailable x^{rhs}"
                )
            available_coeff.add(rhs)

    stored_power_depths = raw.get("power_depths", [])
    if stored_power_depths:
        available_depths = normalize_power_depths(stored_power_depths)
    elif raw.get("made_powers"):
        available_depths = normalize_power_depths(raw.get("made_powers"))

    xi = XI(
        bool(raw.get("multA", False)),
        int(raw.get("n", 0)),
    )
    xi.add_routes(
        routes,
        int(raw.get("depth", 0)),
        int(raw.get("pmult", 0)),
        int(raw.get("cmult", 0)),
        set(available_depths),
        route_ops=route_ops,
        power_depths=available_depths,
    )
    return xi


def serialize_decomp(dcmp: Decomp | None) -> dict | None:
    if dcmp is None:
        return None

    is_leaf = (
        dcmp.dcmp_p is None
        and dcmp.dcmp_q is None
    )

    return {
        "is_leaf": is_leaf,
        "coeff_type": list(dcmp.coeff_type),
        "eval_order": str(dcmp.eval_order),
        "complexity": serialize_complexity(dcmp.comp),
        "xi": serialize_xi(dcmp.xi),
        "term_plans": [
            serialize_term_plan(plan)
            for plan in dcmp.term_plans
        ],
        "p": (
            None
            if is_leaf
            else serialize_decomp(dcmp.dcmp_p)
        ),
        "q": (
            None
            if is_leaf
            else serialize_decomp(dcmp.dcmp_q)
        ),
    }


def reconstruct_decomp(
    plan: dict,
    poly: Poly | list | tuple,
    incoming_powers: set[int] | dict[int, int] | None = None,
    depth: int = 0,
    symbolic: bool = False,
    validate_type: bool = False,
) -> Decomp:
    '''
    기능:
        직렬화된 decomposition을 탐색 없이 Python Decomp 트리로 복원한다.
    입력:
        plan: serialize_decomp으로 생성된 decomposition dictionary.
        poly: 현재 노드에 적용할 수치 또는 타입 다항식.
        incoming_powers: 현재 노드에서 사용 가능한 pure power 집합 또는 depth dictionary.
        depth: 오류 메시지에 사용할 현재 재귀 깊이.
        symbolic: True이면 자식 plan의 저장 계수타입으로 다항식을 복원한다.
        validate_type: True이면 실제 정규화 타입과 저장 타입의 일치를 검증한다.
    출력:
        복원된 Decomp 객체. 타입이 일치하지 않으면 ValueError.
    '''
    if not isinstance(plan, dict):
        raise ValueError(
            f"Invalid plan at depth {depth}: not an object"
        )

    if not isinstance(poly, Poly):
        poly = Poly(poly)

    stored_coeff_type = plan.get("coeff_type")
    if stored_coeff_type is not None:
        stored_coeff_type = [str(ctype) for ctype in stored_coeff_type]

        if validate_type and poly.coeff_type != stored_coeff_type:
            raise ValueError(
                f"Coefficient type mismatch at depth {depth}: "
                f"{poly.coeff_type} != {stored_coeff_type}"
            )

    incoming_power_depths = normalize_power_depths(incoming_powers)

    xi = deserialize_xi(
        plan.get("xi"),
        incoming_power_depths,
    )
    comp = deserialize_complexity(plan.get("complexity"))

    result = Decomp(poly.coeff, comp, xi)
    if stored_coeff_type is not None:
        result.coeff_type = list(stored_coeff_type)
    result.eval_order = plan.get("eval_order", "terminal")
    result.term_plans = [
        deserialize_term_plan(item)
        for item in plan.get("term_plans", [])
    ]

    if plan.get("is_leaf", False):
        return result

    if xi.n == 0 and not xi.multA:
        raise ValueError(
            f"Invalid non-leaf xi at depth {depth}: n == 0"
        )

    poly_p, poly_q = poly.seperate(xi.n, xi.multA)
    current_power_depths = dict(xi.power_depths)
    dcmp_p = None
    dcmp_q = None

    order = result.eval_order

    if order in {"q_then_p", "q_only"}:
        if not poly_q.is_empty():
            q_plan = plan.get("q")
            if q_plan is None:
                raise ValueError(
                    f"Missing q plan at depth {depth}"
                )
            dcmp_q = reconstruct_decomp(
                q_plan,
                Poly(q_plan["coeff_type"]) if symbolic and q_plan.get("coeff_type") is not None else poly_q,
                current_power_depths,
                depth + 1,
                symbolic,
                validate_type,
            )
            current_power_depths = normalize_power_depths(
                list(current_power_depths.items())
                + list(dcmp_q.merge_power_depths().items())
            )

        if order != "q_only" and not poly_p.is_empty():
            p_plan = plan.get("p")
            if p_plan is None:
                raise ValueError(
                    f"Missing p plan at depth {depth}"
                )
            dcmp_p = reconstruct_decomp(
                p_plan,
                Poly(p_plan["coeff_type"]) if symbolic and p_plan.get("coeff_type") is not None else poly_p,
                current_power_depths,
                depth + 1,
                symbolic,
                validate_type,
            )
    else:
        if not poly_p.is_empty():
            p_plan = plan.get("p")
            if p_plan is None:
                raise ValueError(
                    f"Missing p plan at depth {depth}"
                )
            dcmp_p = reconstruct_decomp(
                p_plan,
                Poly(p_plan["coeff_type"]) if symbolic and p_plan.get("coeff_type") is not None else poly_p,
                current_power_depths,
                depth + 1,
                symbolic,
                validate_type,
            )
            current_power_depths = normalize_power_depths(
                list(current_power_depths.items())
                + list(dcmp_p.merge_power_depths().items())
            )

        if order != "p_only" and not poly_q.is_empty():
            q_plan = plan.get("q")
            if q_plan is None:
                raise ValueError(
                    f"Missing q plan at depth {depth}"
                )
            dcmp_q = reconstruct_decomp(
                q_plan,
                Poly(q_plan["coeff_type"]) if symbolic and q_plan.get("coeff_type") is not None else poly_q,
                current_power_depths,
                depth + 1,
                symbolic,
                validate_type,
            )

    result.update(xi, dcmp_p, dcmp_q)
    result.eval_order = order
    result.term_plans = [
        deserialize_term_plan(item)
        for item in plan.get("term_plans", [])
    ]
    result.made_powers = result.merge_mp()
    result.power_depths = result.merge_power_depths()
    return result


def serialize_search_key(key: tuple) -> dict:
    '''
    기능:
        algorithm.py의 hashable OPD cache key를 JSON dictionary로 변환한다.
    입력:
        key: (poly key, is_type_poly, consider_axi, assumption1, depth_budget) tuple.
    출력:
        JSON 직렬화가 가능한 cache key dictionary.
    '''
    poly_key, is_type_poly, consider_axi, assumption1, depth_budget = key
    return {
        "poly_key": list(poly_key),
        "is_type_poly": bool(is_type_poly),
        "consider_axi": bool(consider_axi),
        "assumption1": bool(assumption1),
        "depth_budget": int(depth_budget),
    }


def deserialize_search_key(raw: dict) -> tuple:
    '''
    기능:
        JSON cache key dictionary를 algorithm.py의 hashable tuple로 복원한다.
    입력:
        raw: serialize_search_key로 생성된 cache key dictionary.
    출력:
        (poly key, is_type_poly, consider_axi, assumption1, depth_budget) tuple.
    '''
    if not isinstance(raw, dict):
        raise ValueError("Search cache key must be an object.")

    is_type_poly = bool(raw.get("is_type_poly", False))
    poly_key = tuple(
        str(value) if is_type_poly else float(value)
        for value in raw.get("poly_key", [])
    )
    degree = max(len(poly_key) - 1, 0)
    leading_is_float = (
        degree > 0
        and is_type_poly
        and str(poly_key[-1]) == "F"
    )
    default_depth_budget = (
        ceil(log2(degree + (1 if leading_is_float else 0)))
        if degree > 0
        else 0
    )
    return (
        poly_key,
        is_type_poly,
        bool(raw.get("consider_axi", True)),
        bool(raw.get("assumption1", True)),
        int(raw.get("depth_budget", default_depth_budget)),
    )


def serialize_search_cache(cache: dict) -> list[dict]:
    '''
    기능:
        OPD SearchResult cache를 JSON 직렬화가 가능한 entry 목록으로 변환한다.
    입력:
        cache: key별 (power cache, decomposition, complexity) tuple 목록.
    출력:
        cache key와 SearchResult 목록을 포함한 JSON entry 목록.
    '''
    entries = []

    for key, results in cache.items():
        entries.append({
            "key": serialize_search_key(key),
            "results": [
                {
                    "input_power_depths": [
                        [int(power), int(power_depth)]
                        for power, power_depth in sorted(powers)
                    ],
                    "complexity": serialize_complexity(comp),
                    "decomposition": serialize_decomp(dcmp),
                }
                for powers, dcmp, comp in results
            ],
        })

    return entries


def deserialize_search_cache(entries: list[dict]) -> dict:
    '''
    기능:
        JSON entry 목록을 OPD SearchResult cache dictionary로 복원한다.
    입력:
        entries: serialize_search_cache로 생성된 cache entry 목록.
    출력:
        key별 (power cache, decomposition, complexity) tuple 목록.
    '''
    cache = {}

    for index, entry in enumerate(entries or []):
        if not isinstance(entry, dict):
            raise ValueError(f"Search cache entry {index} must be an object.")

        key = deserialize_search_key(entry.get("key", {}))
        poly_key, is_type_poly, _, _, _ = key
        poly = Poly(list(poly_key))
        results = []

        for result_index, raw_result in enumerate(entry.get("results", [])):
            if not isinstance(raw_result, dict):
                raise ValueError(
                    f"Search result {index}:{result_index} must be an object."
                )

            power_depths = normalize_power_depths(
                raw_result.get("input_power_depths", [])
                or raw_result.get("input_powers", [])
            )
            powers = frozenset(power_depths.items())
            plan = raw_result.get("decomposition")
            dcmp = reconstruct_decomp(
                plan,
                poly,
                power_depths,
                symbolic=is_type_poly,
            )
            comp = deserialize_complexity(raw_result.get("complexity"))
            dcmp.comp = comp
            results.append((powers, dcmp, comp))

        cache[key] = results

    return cache


def save_opd_caches(initial_cache: dict, intermediate_cache: dict, filename: str) -> None:
    '''
    기능:
        초기/중간 OPD SearchResult cache를 하나의 JSON 파일로 저장한다.
    입력:
        initial_cache: 입력 pure power가 {0,1}인 SearchResult cache.
        intermediate_cache: 그 외 입력 pure power의 SearchResult cache.
        filename: 저장할 JSON 파일 경로.
    출력:
        없음. 지정한 경로에 search cache JSON 파일을 생성한다.
    '''
    data = {
        "initial": serialize_search_cache(initial_cache),
        "intermediate": serialize_search_cache(intermediate_cache),
    }

    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, sort_keys=True)


def load_opd_caches(filename: str) -> tuple[dict, dict]:
    '''
    기능:
        JSON 파일에서 초기/중간 OPD SearchResult cache를 복원한다.
    입력:
        filename: save_opd_caches로 저장한 JSON 파일 경로.
    출력:
        (initial_cache, intermediate_cache) dictionary tuple.
        파일이 없으면 두 개의 빈 dictionary를 반환한다.
    '''
    if not os.path.exists(filename):
        return {}, {}

    with open(filename, 'r', encoding='utf-8') as f:
        data = json.load(f)

    initial_cache = deserialize_search_cache(data.get("initial", []))
    intermediate_cache = deserialize_search_cache(data.get("intermediate", []))
    return initial_cache, intermediate_cache


def select_decomp_for_coeff(poly: Poly | list[float], initial_cache: dict, consider_axi: bool=True, assumption1: bool=False) -> Decomp:
    '''
    기능:
        실제 계수의 ax^i 정규화 타입과 일치하는 cache 후보 중 최적 decomposition을 선택한다.
    입력:
        poly: 실제 수치 계수를 가진 Poly 객체 또는 계수 목록.
        initial_cache: 타입 입력으로 생성한 초기 SearchResult cache.
        consider_axi: cache 생성 시 사용한 ax^i 후보 탐색 여부.
        assumption1: cache 생성 시 사용한 Assumption 1 여부.
    출력:
        실제 정규화 타입과 일치하는 최소 complexity Decomp.
        일치 후보가 없으면 NULL_DECOMP.
    '''
    if not isinstance(poly, Poly):
        poly = Poly(poly)

    key = (
        tuple(poly.coeff_type),
        True,
        bool(consider_axi),
        bool(assumption1),
        ceil(
            log2(
                poly.deg
                + (1 if poly.deg > 0 and poly.coeff_type[poly.deg] == "F" else 0)
            )
        ) if poly.deg > 0 else 0,
    )
    valid_results = []

    for powers, dcmp, _ in initial_cache.get(key, []):
        if powers != frozenset({(0, 0), (1, 0)}):
            continue

        try:
            result = reconstruct_decomp(
                serialize_decomp(dcmp),
                poly,
                {0: 0, 1: 0},
                validate_type=True,
            )
        except ValueError:
            continue

        valid_results.append(result)

    if len(valid_results) == 0:
        return NULL_DECOMP

    valid_results.sort()
    return valid_results[0]


# 전역 변수 혹은 클래스 멤버로 캐시 로드
def load_cache(CACHE_FILE: str):
    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, 'r') as f:
            cache = json.load(f)

        # 기존 A0| / A1| 키를 계수타입만 사용하는 신규 키로 변환한다.
        cache = {
            key[3:] if key.startswith(("A0|", "A1|")) else key: value
            for key, value in cache.items()
        }
    return cache


def format_cache_json(cache: dict) -> str:
    '''
    기능:
        최상위 cache와 decomposition root 및 p/q 노드는 여러 줄로 유지하고,
        term_plans는 plan별로 줄바꿈하며 그 외 dictionary와 list는 한 줄로
        압축한 JSON 문자열을 생성한다.
    입력:
        cache: 계수타입 문자열을 decomposition dictionary에 대응시킨 cache.
    출력:
        정렬된 키와 선택적 줄바꿈 형식을 적용한 JSON 문자열.
    '''
    def render(value, level: int, parent_key: str | None = None, expand: bool = False) -> str:
        '''
        기능:
            JSON 값을 현재 들여쓰기와 부모 키에 따라 재귀적으로 문자열화한다.
        입력:
            value: 직렬화할 JSON 호환 값.
            level: 현재 들여쓰기 깊이.
            parent_key: 현재 값에 대응하는 부모 dictionary의 키.
            expand: 현재 dictionary를 여러 줄로 출력할지 여부.
        출력:
            지정된 줄바꿈 규칙을 적용한 JSON 문자열.
        '''
        if isinstance(value, dict):
            if len(value) == 0:
                return "{}"

            items = []
            for key in sorted(value):
                child = value[key]
                child_expand = (
                    (
                        isinstance(child, dict)
                        and key in {"p", "q"}
                    )
                    or (
                        isinstance(child, list)
                        and key == "term_plans"
                    )
                )
                rendered = render(
                    child,
                    level + 1,
                    str(key),
                    child_expand,
                )
                items.append(
                    f"{json.dumps(str(key), ensure_ascii=False)}: {rendered}"
                )

            if not expand:
                return "{" + ", ".join(items) + "}"

            indent = "  " * level
            child_indent = "  " * (level + 1)
            return (
                "{\n"
                + ",\n".join(child_indent + item for item in items)
                + "\n"
                + indent
                + "}"
            )

        if isinstance(value, list):
            if len(value) == 0:
                return "[]"

            items = [
                render(item, level + 1)
                for item in value
            ]

            if not expand:
                return "[" + ", ".join(items) + "]"

            indent = "  " * level
            child_indent = "  " * (level + 1)
            return (
                "[\n"
                + ",\n".join(child_indent + item for item in items)
                + "\n"
                + indent
                + "]"
            )

        return json.dumps(value, ensure_ascii=False)

    if len(cache) == 0:
        return "{}"

    root_items = []
    for key in sorted(cache):
        root_items.append(
            "  "
            + json.dumps(str(key), ensure_ascii=False)
            + ": "
            + render(cache[key], 1, str(key), isinstance(cache[key], dict))
        )

    return "{\n" + ",\n".join(root_items) + "\n}"


def save_cache(cache: dict, filename: str):
    parent_dir = os.path.dirname(os.path.abspath(filename))
    os.makedirs(parent_dir, exist_ok=True)

    normalized_cache = {
        key[3:] if key.startswith(("A0|", "A1|")) else key: value
        for key, value in cache.items()
    }

    with open(filename, 'w', encoding='utf-8') as f:
        f.write(format_cache_json(normalized_cache))
        f.write('\n')

def cal_polyEval_cached(poly: Poly, cache: dict, assumption1: bool=True) -> Decomp:
    '''
    기능:
        계수타입과 Assumption 1 설정에 맞는 OPD 결과를 캐시에서 복원하거나 계산한다.
    입력:
        poly: 평가할 Poly 객체.
        cache: 직렬화된 decomposition을 저장한 dictionary.
        assumption1: ax^i 정규화에 Assumption 1을 적용할지 여부.
    출력:
        캐시에서 복원하거나 새로 계산한 최적 Decomp 객체.
    '''
    # 1. 키 생성
    key = get_poly_type_key(poly.coeff_type, assumption1)

    # 2. 캐시 히트 (Cache Hit)
    if key in cache:
        # 저장된 구조를 불러와서 현재 계수(poly)에 적용
        plan = cache[key]
        return reconstruct_decomp(plan, poly)

    # 3. 캐시 미스 (Cache Miss) -> 기존 알고리즘 실행
    # 기존 로직 수행 (cal_polyEval 내부 로직)
    result: Decomp = cal_polyEval(poly, assumption1=assumption1)

    # 4. 결과 저장
    if result:
        cache[key] = serialize_decomp(result)
        # 주기적으로 save_cache() 호출 필요 (또는 프로그램 종료 시)

    return result
