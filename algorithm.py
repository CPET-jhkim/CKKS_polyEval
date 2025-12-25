# algorithm.py
from complexity import Complexity, attach, compare
from polynomial import Poly
from math import log2, ceil, sqrt
from print import print_step
from util import *
from collections import deque

def calculate(poly: Poly, made_powers: set[int]) -> tuple[Complexity, set[int], tuple[bool, int, tuple, tuple]]:
    max_deg = poly.deg
    ct = poly.coeff_type
    # =============================
    '''
    0. 분해가 필요없는 기초다항식 처리
        빈 다항식의 경우 복잡도 0,0,0,0을 반환.
    '''
    ## 0-1) 빈 다항식 or 0차식 처리
    if max_deg <= 0:
        return Complexity(), made_powers, (False, 0, (), ())
    
    ## 0-2) 1차식 처리
    elif max_deg == 1:
        comp_res = Complexity()
        comp_res.depth = 0 if ct[-1] == "I" else 1
        comp_res.cmult = 0
        comp_res.pmult = 0 if ct[-1] == "I" else 1
        comp_res.add = 1 if ct[0] != "0" else 0
        return comp_res, made_powers, (False, 0, (), ())
    
    # ==============================
    '''
    1. 분해하지 않는 경우
        x^0, x^1, ..., x^n 항을 모두 따로 계산할 때의 복잡도 측정.
        depth:  계수가 정수이면 log2(n)이며, 소수일 경우 log2(n+1).
        cmult:  x^0, x^1, ..., x^n을 각각 만들기 위한 최소 비용 알고리즘을 통해 측정.
        pmult:  (0, 정수)가 아닌 계수의 개수
        add:    덧셈 횟수(= 계수가 0이 아닌 총 항의 개수)
    '''
    comp_temp = Complexity()
    
    ## 1-1) depth
    if poly.coeff_type[-1] == "I":
        comp_temp.depth = ceil(log2(max_deg))
    else:
        comp_temp.depth = ceil(log2(max_deg+1))
    
    ## 1-2) cmult
    final_powers, ops, add_count = build_addition_plan_for_poly(poly.coeff, made_powers)
    comp_temp.cmult = add_count
    
    ## 1-3) pmult
    for i in range(1, max_deg+1):
        if ct[i] == "F":
            comp_temp.pmult += 1

    ## 1-4) add
    rp = {i for i, c in enumerate(poly.coeff) if c != 0}
    comp_temp.add = len(rp) - 1
    
    results = [(comp_temp, final_powers, (False, 0, (), ()))]
    
    '''
    2. f(x) = (x^i)*p(x) + q(x)로 분해하는 경우
        A. (x^i)의 계산복잡도 측정(made_powers 업데이트)
        B. p(x), q(x) 중 최고차수가 가장 낮은 것부터 복잡도 측정.
            이 때 (x^i)를 구성할 수 있는 made_powers 각각에 대하여 측정.
        C. 모든 항의 복잡도 측정 후 결합
            (x^i) * p(x)의 경우 곱셈 -> depth(큰쪽) + 1, cmult(합) + 1
            (x^i)p(x) + q(x)의 경우 덧셈 -> add + 1        
    '''
    for i in range(max_deg-1, 0, -1):
        mp = made_powers.copy()
        
        '''
        ax^i로 분해하는 조건
        1. 계수 a가 소수
        2. x^i , p(x) 중 x^i에 붙이는 것이 효율적인 경우
        '''
        multA = ct[max_deg] == "F" and check_axi_optimize(i, max_deg-i)
            
        xi_routes = cal_xi_routes(i, mp, multA)
        poly_p, poly_q = poly.seperate(i, multA)
        
        for depth, mp2, ops_list in xi_routes:
            add_count = len(ops_list)
            # x^i의 계산복잡도 정리
            comp_i = Complexity()
            pm = 1 if multA else 0
            comp_i.insert_value(depth, add_count, pm, 0)
            
            # 차수가 작은 다항식부터 연산
            if poly_p.deg < poly_q.deg:
                comp_p, mp3, decomp_p = calculate(poly_p, mp2)
                comp_q, mp4, decomp_q = calculate(poly_q, mp3)
            else:
                comp_q, mp3, decomp_q = calculate(poly_q, mp2)
                comp_p, mp4, decomp_p = calculate(poly_p, mp3)
                
            # 계산복잡도 결합
            comp_pi = attach(comp_i, comp_p, 'x')
            if poly_q.coeff != []:
                comp_piq = attach(comp_pi, comp_q, '+')
            else:
                comp_piq = comp_pi
                
            results.append((comp_piq, mp, (multA, i, decomp_p, decomp_q)))

    '''
    3. 최적의 결과 비교
        모든 분해식 결과는 (x^i 차수, 총 계산복잡도, 만들어진 x^i 차수들) tuple로 구성.
        총 계산복잡도 기준으로 depth, cmult, pmult, add 순으로 가장 낮은 값을 검색.
    '''
    best = results[0]
    for c in results[1:]:
        if compare(best[0], c[0]) == 2:
            best = c    

    return best

# # 다항식에서 0이 아닌 값의 index를 반환하는 함수
# def required_powers(poly: list[float]) -> set[int]:
#     return {i for i, c in enumerate(poly) if c != 0}


# x^i를 구성하기 위한 모든 경우의 수를 탐색.
# BFS방식으로 구현.
def cal_xi_routes(target_i: int, made_powers: set[int], multA: bool = False) -> list[tuple[int, set[int], list]]:
    """
    x^i를 구성할 수 있는 최적 경로 탐색 (BFS)
    :param multA: True일 경우 ax^i 연산을 고려하여, 마지막 단계에서 a가 곱해지지 않은 항의 집합을 반환
    :return: list of (depth, result_set, operations_list)
    """
    # 이미 목표 값이 존재하는 경우
    if target_i in made_powers:
        depth = 1 if multA else 0
        return [(depth, made_powers, [])]

    # {숫자: 깊이} (기존 숫자는 깊이 0으로 가정)
    initial_depths = {num: 0 for num in made_powers}
    queue = deque([(initial_depths, [])])
    
    visited_states = set()
    visited_states.add(tuple(sorted(initial_depths.items())))
    
    candidates = [] 
    min_ops_len = float('inf') 
    
    while queue:
        curr_depths, history = queue.popleft()
        
        # 가지치기: 현재 경로가 이미 찾은 최단 경로보다 길면 중단
        if len(history) > min_ops_len:
            continue
            
        curr_elements = list(curr_depths.keys())
        
        # 목표 값을 찾았는지 확인 (BFS 특성상 큐에서 꺼낼 때 확인해야 함)
        if target_i in curr_depths:
            target_depth = curr_depths[target_i]
            current_len = len(history)
            
            if current_len <= min_ops_len:
                min_ops_len = current_len
                
                # --- [수정된 부분] multA 로직 분기 ---
                if multA and history:
                    # 마지막 연산이 u + v = target_i 일 때,
                    # a가 u에 붙거나 v에 붙는 두 가지 경우를 분리
                    u, v, _ = history[-1]
                    
                    # 1. v에 a가 붙은 경우 -> 순수 집합은 {u}
                    candidates.append((target_depth, {u}, history))
                    
                    # 2. u != v 일 때, u에 a가 붙은 경우 -> 순수 집합은 {v}
                    if u != v:
                        candidates.append((target_depth, {v}, history))
                else:
                    # 기본 동작: 현재까지 만든 모든 수 집합 반환
                    candidates.append((target_depth, set(curr_elements), history))
                # -----------------------------------
                
            continue
            
        # 다음 단계 탐색 (Addition Chain 생성)
        next_moves = []
        n = len(curr_elements)
        for j in range(n):
            for k in range(j, n):
                val_a, val_b = curr_elements[j], curr_elements[k]
                new_val = val_a + val_b
                
                # target_i보다 작거나 같고, 아직 만들지 않은 수라면 추가
                if new_val <= target_i and new_val not in curr_depths:
                    new_depth = max(curr_depths[val_a], curr_depths[val_b]) + 1
                    next_moves.append((new_val, new_depth, val_a, val_b))
        
        for val, d, a, b in next_moves:
            new_depths = curr_depths.copy()
            new_depths[val] = d
            state_key = tuple(sorted(new_depths.items()))
            
            if state_key not in visited_states:
                visited_states.add(state_key)
                new_hist = history + [(a, b, val)]
                queue.append((new_depths, new_hist))
    
    if not candidates:
        return []

    # 필터링: 1. Depth 최소, 2. Count(연산 횟수) 최소
    min_depth = min(c[0] for c in candidates)
    depth_filtered = [c for c in candidates if c[0] == min_depth]
    
    min_count = min(len(c[2]) for c in depth_filtered)
    final_results = [c for c in depth_filtered if len(c[2]) == min_count]
    
    return final_results

def build_addition_plan_for_poly(poly: list, made_powers: set[int]):
    """
    다항식 계산을 위한 최적 덧셈 사슬 계획 수립 (재귀적 탐색)
    """
    # 0, 1차항은 기본 제공으로 간주하고 그 이상 차수만 계산
    required = sorted([i for i, c in enumerate(poly) if c != 0 and i > 1])
    
    # 초기 집합
    base_powers = {0, 1} | set(made_powers)
    
    if not required:
        return base_powers, [], 0
    
    # 전역 최적해 저장용
    global_best = {
        'ops': None,
        'count': float('inf')
    }
    
    def recursive_solver(target_idx, current_set, accumulated_ops):
        # Base Case: 모든 타겟 해결
        if target_idx == len(required):
            if len(accumulated_ops) < global_best['count']:
                global_best['count'] = len(accumulated_ops)
                global_best['ops'] = accumulated_ops
            return

        # Pruning: 이미 최적해보다 나빠졌으면 중단
        if len(accumulated_ops) >= global_best['count']:
            return

        target = required[target_idx]
        
        # 이미 있다면 다음으로
        if target in current_set:
            recursive_solver(target_idx + 1, current_set, accumulated_ops)
            return

        # BFS를 통해 현재 상태에서 target을 만드는 후보 경로들 탐색
        candidates = cal_xi_routes(target, current_set)
        
        # Branching
        for _, new_set, ops in candidates:
            recursive_solver(target_idx + 1, new_set, accumulated_ops + ops)

    # 재귀 실행
    recursive_solver(0, base_powers, [])
    
    final_ops = global_best['ops'] if global_best['ops'] is not None else []
    final_count = global_best['count'] if global_best['count'] != float('inf') else 0
    
    # 연산 순서대로 집합 복원
    final_powers = set(base_powers)
    for a, b, c in final_ops:
        final_powers.add(c)
        
    return final_powers, final_ops, final_count