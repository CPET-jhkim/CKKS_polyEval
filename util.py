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


# ax^i검증 -> a를 x^i에 붙이는 것이 효율적인가??
# def check_axi_optimize(i: int, j: int) -> tuple[bool, bool]:
#     '''
#     판단기준
#     x^i의 깊이 증가성
#     x^j의 깊이 증가성
    
#     반환
#     x^i에 a결합 여부, x^j에 a결합 여부
#     '''
#     di = ceil(log2(i)) < ceil(log2(i+1))
#     dj = ceil(log2(j)) < ceil(log2(j+1))
    
#     if dj:
        
    
def check_axi_optimize(i: int, j: int) -> bool:
    """
    ax^i * p(x) (p(x)의 대표 차수 x^j) 연산 시,
    a를 x^i 쪽에 붙이는 것이 더 효율적인지 판별.
    
    판단 기준:
    1. '중간 삽입(최적화)'이 가능한 쪽을 우선 (비트 1이 2개 이상)
    2. 조건이 같다면, 차수가 더 낮은 쪽을 우선 (깊이가 얕음)
    
    :param i: x^i의 지수
    :param j: p(x)의 대표 지수 (일반적으로 분할된 나머지 차수)
    :return: True (a를 x^i에 붙임), False (a를 p(x)쪽 x^j에 붙임)
    """
    # 0. 예외 처리
    if i <= 0: return False # i가 없으면 j에 붙여야 함
    if j <= 0: return True  # j가 없으면 i에 붙여야 함

    # 1. 최적화 가능 여부 (Popcount > 1) 확인
    # (n & (n-1)) != 0 이면 2의 거듭제곱이 아님 -> 최적화 가능
    opt_i = (i & (i - 1)) != 0
    opt_j = (j & (j - 1)) != 0

    # Case 1: 한쪽만 최적화 가능한 경우 -> 가능한 쪽 선택
    if opt_i and not opt_j:
        return True
    if not opt_i and opt_j:
        return False
    
    # Case 2: 둘 다 최적화 가능하거나, 둘 다 불가능한 경우 (Tie-breaking)
    # -> 차수가 낮은 쪽에 붙이는 것이 유리함.
    # -> 문제의 i=1, j=2 케이스: 둘 다 최적화 불가능(False) -> 1 < 2 이므로 True 반환
    return i < j