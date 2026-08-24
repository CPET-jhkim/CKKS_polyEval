# save.py
import json, os
import numpy as np
from basic_class import Poly, Decomp, XI, Complexity
from algorithm import attach, cal_polyEval

# 다항식 계수타입을 키로 반환
def get_poly_type_key(poly_type: list) -> str:
    if type(poly_type[0]) is not str:
        poly_type = Poly(poly_type).coeff_type
    return "".join(poly_type)




CACHE_SCHEMA_VERSION = 2


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
        "made_powers": sorted(
            int(power) for power in xi.made_powers
        ),
        "depth": int(xi.depth),
        "cmult": int(xi.cmult),
        "pmult": int(xi.pmult),
    }


def deserialize_xi(
    raw: dict | None,
    incoming_powers: set[int],
) -> XI:
    if raw is None:
        xi = XI()
        xi.add_routes(
            [],
            0,
            0,
            0,
            set(incoming_powers),
            route_ops=[],
        )
        return xi

    route_ops = deserialize_route_ops(raw.get("route_ops", []))

    # schema v1 호환: route_ops가 없으면 기존 route를 pure로 간주한다.
    if not route_ops:
        route_ops = [
            ("pure", int(lhs), int(rhs))
            for lhs, rhs in raw.get("route", [])
        ]

    routes = [(lhs, rhs) for _, lhs, rhs in route_ops]
    available = set(incoming_powers)

    for index, (kind, lhs, rhs) in enumerate(route_ops):
        if lhs not in available:
            raise ValueError(
                f"route_ops[{index}] uses unavailable x^{lhs}"
            )
        if rhs not in available:
            raise ValueError(
                f"route_ops[{index}] uses unavailable x^{rhs}"
            )
        if kind == "pure":
            available.add(lhs + rhs)

    stored_made_powers = {
        int(power)
        for power in raw.get("made_powers", [])
    }
    if stored_made_powers:
        available |= stored_made_powers

    xi = XI(
        bool(raw.get("multA", False)),
        int(raw.get("n", 0)),
    )
    xi.add_routes(
        routes,
        int(raw.get("depth", 0)),
        int(raw.get("pmult", 0)),
        int(raw.get("cmult", 0)),
        available,
        route_ops=route_ops,
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
        "schema_version": CACHE_SCHEMA_VERSION,
        "is_leaf": is_leaf,
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
    poly: Poly | list | np.ndarray,
    incoming_powers: set[int] | None = None,
    depth: int = 0,
) -> Decomp:
    """schema v2 캐시를 탐색 없이 Python Decomp로 복원한다."""
    if not isinstance(plan, dict):
        raise ValueError(
            f"Invalid plan at depth {depth}: not an object"
        )

    if not isinstance(poly, Poly):
        poly = Poly(poly)

    if incoming_powers is None:
        incoming_powers = {0, 1}
    else:
        incoming_powers = set(incoming_powers)

    xi = deserialize_xi(
        plan.get("xi"),
        incoming_powers,
    )
    comp = deserialize_complexity(plan.get("complexity"))

    result = Decomp(poly.coeff, comp, xi)
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
    current_powers = set(xi.made_powers)
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
                poly_q,
                current_powers,
                depth + 1,
            )
            current_powers |= dcmp_q.merge_mp()

        if order != "q_only" and not poly_p.is_empty():
            p_plan = plan.get("p")
            if p_plan is None:
                raise ValueError(
                    f"Missing p plan at depth {depth}"
                )
            dcmp_p = reconstruct_decomp(
                p_plan,
                poly_p,
                current_powers,
                depth + 1,
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
                poly_p,
                current_powers,
                depth + 1,
            )
            current_powers |= dcmp_p.merge_mp()

        if order != "p_only" and not poly_q.is_empty():
            q_plan = plan.get("q")
            if q_plan is None:
                raise ValueError(
                    f"Missing q plan at depth {depth}"
                )
            dcmp_q = reconstruct_decomp(
                q_plan,
                poly_q,
                current_powers,
                depth + 1,
            )

    result.update(xi, dcmp_p, dcmp_q)
    result.eval_order = order
    result.term_plans = [
        deserialize_term_plan(item)
        for item in plan.get("term_plans", [])
    ]
    result.made_powers = result.merge_mp()
    return result

# 전역 변수 혹은 클래스 멤버로 캐시 로드
def load_cache(CACHE_FILE: str):
    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, 'r') as f:
            cache = json.load(f)
    return cache

def save_cache(cache: dict, filename: str):
    with open(filename, 'w') as f:
        json.dump(cache, f, indent=2, sort_keys=True)

def cal_polyEval_cached(poly: Poly, cache: dict) -> Decomp:
    # 1. 키 생성
    key = get_poly_type_key(poly.coeff_type)
    
    # 2. 캐시 히트 (Cache Hit)
    if key in cache:
        # 저장된 구조를 불러와서 현재 계수(poly)에 적용
        plan = cache[key]
        return reconstruct_decomp(plan, poly)
    
    # 3. 캐시 미스 (Cache Miss) -> 기존 알고리즘 실행
    # 기존 로직 수행 (cal_polyEval 내부 로직)
    result: Decomp = cal_polyEval(poly) 
    
    # 4. 결과 저장
    if result:
        cache[key] = serialize_decomp(result)
        # 주기적으로 save_cache() 호출 필요 (또는 프로그램 종료 시)

    return result


if __name__ == '__main__':
    from poly_util import make_all_polys
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from tqdm import tqdm
    
    # nohup /home/doodle/.virtualenvs/remez/bin/python /home/doodle/python/py_approx_test/polyEval/save.py > /dev/null 2>&1 &
    # n~m차까지 모든 다항식 타입에 대하여 저장
    load_filename = 'data/decomp_cache.json'
    save_filename = 'data/decomp_cache.json'
    dcmp_cache = load_cache(load_filename)
    for deg in range(1, 9):
        polys = make_all_polys(deg)
        total_polys = len(polys)
        
        print(f"Degree {deg}: Processing {total_polys} polys...")
        
        try:
            future_to_key = {}
            # 에러 발생 시 즉시 중단을 위해 context manager 내에서 엄격하게 관리
            with ProcessPoolExecutor(max_workers=30) as executor:
                # 1. 태스크 제출
                for coeff in polys:
                    key = get_poly_type_key(coeff)
                    future = executor.submit(cal_polyEval, coeff)
                    future_to_key[future] = key
                
                # 2. 결과 수집 및 에러 감시
                with tqdm(total=total_polys, desc=f"Deg {deg:02d}", unit="poly", ascii=True) as pbar:
                    for future in as_completed(future_to_key):
                        key = future_to_key[future]
                        # future.result() 호출 시 해당 프로세스에서 발생한 예외가 raise됨
                        result = future.result() 
                        
                        if type(result) == Decomp:
                            data = serialize_decomp(result)
                            dcmp_cache[key] = data
                        
                        pbar.update(1)
            
            # 정상 종료 시 저장
            save_cache(dcmp_cache, save_filename)
            print(f"deg={deg} save completed.\n" + "-"*30)
            
        except Exception as e:
            # 개별 태스크 혹은 전체 프로세스 에러 시 catch
            print(f"\n[CRITICAL] Error occurred at Degree {deg}, Key {key}: {e}")
            # 루프 전체를 종료하기 위해 예외를 다시 던지거나 sys.exit() 사용
            raise e
        
