# util.py
# 유틸리티 함수들

import random
from itertools import product
from math import log2, ceil

def make_all_polys(max_deg: int) -> list[list[float]]:
    """
    max_deg에 대해 가능한 '모든 계수 타입의 조합'을 구성하고,
    각 조합에 대해 실제 값은 랜덤하게 생성하여 반환합니다.
    
    [조합 규칙]
    - 0 ~ (max_deg-1) 차수: [0, 정수, 소수] 3가지 경우의 수
    - max_deg 차수(최고차항): [정수, 소수] 2가지 경우의 수
    
    [값 생성 규칙]
    - 0: 0.0 고정
    - 정수: 1~10 사이 랜덤 정수 (float 변환)
    - 소수: 1.1 ~ 9.9 사이 (소수점 첫째 자리가 0이 되지 않도록 강제)
    """
    
    # 1. 값을 생성하는 익명 함수(lambda) 정의
    gen_zero = lambda: 0.0
    gen_int  = lambda: float(random.randint(1, 10))
    
    # [수정됨] 소수 생성 로직 변경
    # 정수부(1~9) + 소수부(0.1~0.9)를 더하여 무조건 X.1 ~ X.9 형태가 되도록 함
    # 예: 1 + 0.1 = 1.1, 9 + 0.9 = 9.9 (10.0 등은 생성되지 않음)
    gen_float = lambda: float(random.randint(1, 9)) + random.randint(1, 9) / 10.0

    # 2. 각 차수별 가능한 생성기(Generator) 목록 정의
    # 중간 차수용 후보군 (3개)
    generators_middle = [gen_zero, gen_int, gen_float]
    # 최고 차수용 후보군 (2개)
    generators_highest = [gen_int, gen_float]

    # 3. itertools.product를 위한 인자 리스트 구성
    # [중간, 중간, ..., 중간, 최고] 형태의 리스트 생성
    iterables = [generators_middle] * max_deg + [generators_highest]
    
    # 4. 모든 생성기 조합(Skeleton) 생성
    generator_combinations = product(*iterables)

    # 5. 각 조합을 순회하며 실제 랜덤 값을 생성하여 리스트로 변환
    all_polys = []
    for combo in generator_combinations:
        # combo는 (gen_zero, gen_int, ...) 같은 함수들의 튜플입니다.
        # 각 함수를 실행(func())하여 실제 랜덤 값을 뽑아냅니다.
        poly_values = [func() for func in combo]
        all_polys.append(poly_values)
        
    return all_polys
        
def is_multA_required(n: int) -> bool:
    """
    n이 주어졌을 때 a를 중간에 결합해야 하는지(multA=True) 판별합니다.
    
    Returns:
        True: n이 2의 거듭제곱이 아님 (예: 3, 5, 6, 7...) -> 중간 결합 필수
        False: n이 2의 거듭제곱임 (예: 1, 2, 4, 8...) -> 상관없음(나중에 결합 가능)
    """
    if n < 1:
        return False
        
    # (n & (n-1)) == 0 이면 2의 거듭제곱입니다.
    # 따라서 0이 아니면(!=) 2의 거듭제곱이 아니므로 True를 반환합니다.
    return (n & (n - 1)) != 0


def solve_xn_operation(multA: bool, n: int) -> list[tuple[set[int], dict[int, tuple[list[int], tuple[int, int]]]]]:
    """
    n을 만들기 위한 최적 연산 과정을 구하고, 
    각 과정에서 사용된 pure_x_set의 요소들(i)에 대해 개별적인 생성 비용(chain, cost)을 매핑하여 반환합니다.
    """
    
    # ---------------------------------------------------------
    # 1. Helper: 특정 k값(Pure)을 만드는 최적 경로/비용 계산 함수
    # ---------------------------------------------------------
    def get_pure_k_info(k: int) -> tuple[list[int], tuple[int, int]]:
        if k == 1: return ([1], (0, 0))
        
        # DP로 k까지의 최소 Depth 계산
        depth_map = {1: 0}
        for i in range(2, k + 1):
            min_d = float('inf')
            for j in range(1, i // 2 + 1):
                d = max(depth_map[j], depth_map[i-j]) + 1
                if d < min_d: min_d = d
            depth_map[i] = min_d
            
        # 최적 경로 역추적 (Min Cmult 우선)
        best_res = None
        
        def backtrack(target, path, cmult):
            nonlocal best_res
            if best_res and cmult >= best_res[1][1]: return
            
            if target == 1:
                res_chain = sorted(path + [1])
                res_cost = (depth_map[k], cmult)
                if best_res is None or cmult < best_res[1][1]:
                    best_res = (res_chain, res_cost)
                return

            # 분할 탐색
            for j in range(1, target // 2 + 1):
                rem = target - j
                # Critical Path 유지 조건
                if max(depth_map[j], depth_map[rem]) + 1 == depth_map[target]:
                    backtrack(rem, path + [target], cmult + 1)
                    
        backtrack(k, [], 0)
        return best_res if best_res else ([1, k], (depth_map[k], 1))

    # ---------------------------------------------------------
    # 2. Main Logic: 전체 n에 대한 Global 최적 경로 탐색
    # ---------------------------------------------------------
    
    # DP 초기화
    min_depths = {1: (0, 0)} # (pure_depth, coeff_depth)
    for i in range(2, n + 1):
        d_pure, d_coeff = float('inf'), float('inf')
        for j in range(1, i // 2 + 1):
            k = i - j
            # Pure
            p_val = max(min_depths[j][0], min_depths[k][0]) + 1
            if p_val < d_pure: d_pure = p_val
            # Coeff (multA=True)
            if multA:
                c_val = max(min_depths[j][1], min_depths[k][0]) + 1
                if c_val < d_coeff: d_coeff = c_val
                c_val2 = max(min_depths[k][1], min_depths[j][0]) + 1
                if c_val2 < d_coeff: d_coeff = c_val2
                
        min_depths[i] = (d_pure, d_coeff)

    target_opt_depth = min_depths[n][1] if multA else min_depths[n][0]
    
    def find_global_paths(target, has_coeff):
        if target == 1:
            return [{'pure': {1}, 'depth': 0, 'cmult': 0}]
            
        res = []
        for j in range(1, target // 2 + 1):
            k = target - j
            
            # Case 1: Pure 생성
            if not has_coeff:
                if max(min_depths[j][0], min_depths[k][0]) + 1 <= min_depths[target][0]:
                    l_list = find_global_paths(j, False)
                    r_list = find_global_paths(k, False)
                    for l in l_list:
                        for r in r_list:
                            nd = max(l['depth'], r['depth']) + 1
                            if nd == min_depths[target][0]:
                                res.append({
                                    'pure': l['pure'] | r['pure'] | {target},
                                    'depth': nd,
                                    'cmult': l['cmult'] + r['cmult'] + 1
                                })
            
            # Case 2: Coeff 생성 (has_coeff True)
            else:
                # ax^j * x^k
                if max(min_depths[j][1], min_depths[k][0]) + 1 <= target_opt_depth:
                    c_list = find_global_paths(j, True)
                    p_list = find_global_paths(k, False)
                    for c in c_list:
                        for p in p_list:
                            nd = max(c['depth'], p['depth']) + 1
                            if nd == target_opt_depth:
                                res.append({
                                    'pure': c['pure'] | p['pure'],
                                    'depth': nd,
                                    'cmult': c['cmult'] + p['cmult'] + 1
                                })
                
                # ax^k * x^j (j!=k)
                if j != k:
                    if max(min_depths[k][1], min_depths[j][0]) + 1 <= target_opt_depth:
                        c_list = find_global_paths(k, True)
                        p_list = find_global_paths(j, False)
                        for c in c_list:
                            for p in p_list:
                                nd = max(c['depth'], p['depth']) + 1
                                if nd == target_opt_depth:
                                    res.append({
                                        'pure': c['pure'] | p['pure'],
                                        'depth': nd,
                                        'cmult': c['cmult'] + p['cmult'] + 1
                                    })
        return res

    candidates = find_global_paths(n, multA)
    
    # ---------------------------------------------------------
    # 3. Result Processing
    # ---------------------------------------------------------
    final_output = []
    seen_sigs = set()
    
    if not candidates: return []
    min_total_cmult = min(c['cmult'] for c in candidates)
    
    for cand in candidates:
        if cand['cmult'] > min_total_cmult: continue
        
        pure_set = cand['pure']
        pure_sorted = sorted(list(pure_set))
        
        sig = tuple(pure_sorted)
        if sig in seen_sigs: continue
        seen_sigs.add(sig)
        
        details_map = {}
        for val in pure_sorted:
            details_map[val] = get_pure_k_info(val)
            
        final_output.append((pure_set, details_map))
        
    return final_output  
