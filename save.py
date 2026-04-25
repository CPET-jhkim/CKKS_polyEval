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

# 분해식 데이터를 dict형으로 직렬화
def serialize_decomp(dcmp: Decomp) -> dict | None:
    if dcmp is None:
        return None
    
    # 더 이상 분해되지 않는 말단 노드(Base case) 처리
    # xi.n이 0이거나 coeff가 비어있을 경우 등 로직에 맞게 조정
    if dcmp.xi.n == 0 and dcmp.dcmp_p is None:
        return {
            "is_leaf": True,
            "xi": None
        }

    return {
        "is_leaf": False,
        "xi": {
            "n": dcmp.xi.n,
            "multA": dcmp.xi.multA,
            "route": dcmp.xi.route  # XI 복원을 위해 경로 저장
        },
        "p": serialize_decomp(dcmp.dcmp_p),
        "q": serialize_decomp(dcmp.dcmp_q)
    }

# 저장한 데이터로 탐색 없이 복원
def reconstruct_decomp(plan: dict, poly: Poly | list | np.ndarray, depth: int = 0) -> Decomp:
    indent = "  " * depth
    
    # 0. 입력 데이터 타입 방어 코드
    try:
        if plan is None:
            raise ValueError("Plan is None")
            
        # poly가 list나 ndarray면 Poly 객체로 변환
        if not isinstance(poly, Poly):
            poly = Poly(poly)
    except Exception as e:
        print(f"{indent}[Error] Input Validation Failed at depth {depth}")
        print(f"{indent}Poly type: {type(poly)}")
        print(f"{indent}Plan: {plan}")
        raise e

    try:
        # --- 디버깅용 로그 (필요 없으면 주석 처리) ---
        # print(f"{indent}> Reconstruct Start: deg={poly.deg}, coeff={poly.coeff}")

        # 1. Base Case: 더 이상 분해되지 않는 경우 (Leaf)
        if plan.get("is_leaf", False):
            # 기본 Complexity 생성 (Leaf 노드용)
            # 주의: Leaf 노드일 때도 기본 연산 비용이 있다면 계산 로직 필요
            # 여기서는 기본값으로 처리
            res = Decomp(poly.coeff, Complexity(), XI()) 
            res.made_powers = {0, 1}
            return res

        # 2. 구조 정보 복원 (XI)
        plan_xi = plan.get("xi")
        if not plan_xi:
            raise ValueError(f"Missing 'xi' in plan at depth {depth}")

        xi = XI(plan_xi.get("multA", False), plan_xi.get("n", 1))
        
        # made_powers 복원
        made_powers = {0, 1}
        raw_routes = plan_xi.get("route", [])
        routes = [tuple(r) for r in raw_routes]            
        for op in routes:
            # op가 리스트인지 튜플인지 확인하여 안전하게 합계 계산
            # val = sum(op) if isinstance(op, (list, tuple)) else op
            made_powers.add(sum(op))
        
        # XI 객체에 route 정보 주입
        if hasattr(xi, 'add_routes'):
            xi.add_routes(routes, made_powers)
        else:
            xi.route = routes # fallback

        # 3. 다항식 분리
        # seperate 메서드가 실패하지 않도록 try 감싸기
        try:
            poly_p, poly_q = poly.seperate(xi.n, xi.multA)
        except Exception as sep_err:
            print(f"{indent}[Error] Separation failed: xi.n={xi.n}, multA={xi.multA}")
            raise sep_err

        # 4. 재귀적 복원
        dcmp_p = None
        dcmp_q = None

        # P 복원 (P는 보통 존재해야 함)
        if plan.get("p"):
            dcmp_p = reconstruct_decomp(plan["p"], poly_p, depth + 1)
        else:
            # Plan에는 없는데 Poly P가 비어있지 않다면 문제
            if poly_p.coeff and poly_p.deg >= 0:
                print(f"{indent}[Warning] Plan for P is missing but Poly P exists: {poly_p.coeff}")

        # Q 복원 (Q는 없을 수 있음)
        if poly_q.coeff and plan.get("q"):
            dcmp_q = reconstruct_decomp(plan["q"], poly_q, depth + 1)
        
        # 5. Complexity 병합
        comp_i = Complexity()
        comp_i.insert_value(xi.depth, xi.add_count, xi.pmult, 0)

        # p_comp 안전하게 가져오기
        p_comp = dcmp_p.comp if dcmp_p else Complexity()
        
        # attach 함수 호출
        # (attach 함수가 None을 처리하지 못할 경우를 대비해 p_comp 전달)
        comp_pi = attach(xi, comp_i, poly_p, p_comp, 'x')
        
        comp_piq = comp_pi
        if dcmp_q:
            comp_piq = attach(None, comp_pi, poly_q, dcmp_q.comp, '+')

        # 6. 최종 결과 생성
        res = Decomp(poly.coeff, comp_piq)
        res.update(xi, dcmp_p, dcmp_q)
        
        # made_powers 병합
        res.made_powers = res.merge_mp()
        
        return res

    except Exception as e:
        # 에러 발생 시 상세 문맥 출력
        print(f"\n{'='*20} CRITICAL ERROR at Depth {depth} {'='*20}")
        print(f"{indent}Error Type: {type(e).__name__}")
        print(f"{indent}Message: {e}")
        print(f"{indent}Current Poly Coeff: {poly.coeff}")
        print(f"{indent}Current XI Plan: {plan.get('xi', 'Missing')}")
        
        # Traceback 출력 (선택 사항)
        # traceback.print_exc() 
        
        # 상위 호출 스택으로 에러 전파
        raise e


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
    save_filename = 'data/decomp_cache0408.json'
    dcmp_cache = load_cache(load_filename)
    for deg in range(9, 30):
        polys = make_all_polys(deg)
        total_polys = len(polys)
        
        print(f"Degree {deg}: Processing {total_polys} polys...")
        
        try:
            future_to_key = {}
            # 에러 발생 시 즉시 중단을 위해 context manager 내에서 엄격하게 관리
            with ProcessPoolExecutor(max_workers=1) as executor:
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
        
