# basic_class.py
from .print import poly_to_str, print_poly, print_poly_type
from math import log2, ceil
from functools import total_ordering
from itertools import product

@total_ordering
class Complexity:
    '''
    Calculation complexity class
    '''
    def __init__(self, xi = None):
        if xi is not None:
            self.depth = xi.depth
            self.cmult = xi.cmult
            self.pmult = xi.pmult
            self.add = 0
        else:
            self.depth = 0
            self.cmult = 0
            self.pmult = 0
            self.add = 0
    
    def insert_value(self, depth, cmult, pmult, add):
        self.depth = depth
        self.cmult = cmult
        self.pmult = pmult
        self.add = add

    # Comparison - smaller means higher complexity.
    def __lt__(self, other):
        return (self.depth, self.cmult, self.pmult, self.add) < (other.depth, other.cmult, other.pmult, other.add)
    
    def __eq__(self, other):
        return (self.depth, self.cmult, self.pmult, self.add) == (other.depth, other.cmult, other.pmult, other.add)
    
    def print_params(self):
        print(f"{'DCPA:':<8}{self.depth}|{self.cmult}|{self.pmult}|{self.add}")
    
    def return_params(self):
        return f"{self.depth}|{self.cmult}|{self.pmult}|{self.add}"
        
class Poly:
    '''
    Polynomial info class
    coeff: coefficient info
    complexity: calculation complexity info.
    mp: created x^i info.
    '''
    def __init__(self, coeff: list[float] | list[str]):
        self.is_type_poly = (
            len(coeff) > 0
            and all(type(cf) == str for cf in coeff)
        )

        if self.is_type_poly:
            if not all(cf in {"0", "I", "F"} for cf in coeff):
                raise ValueError("Coefficient type must be one of 0, I, F.")

            self.coeff_type = [str(cf) for cf in coeff]
            self.coeff = [
                0.0 if cf == "0" else float(i + 1) if cf == "I" else float(i + 1) + 0.5
                for i, cf in enumerate(coeff)
            ]
        else:
            self.coeff = [float(cf) for cf in coeff]
            self.coeff_type: list[str] = []
            self.check_type()

        self.deg = max(len(coeff) - 1, 0)
        self.complexity = Complexity()
        # self.mp: set[int] = set([0, 1])
        # self.ops_list = None
    
    # Check each coefficient's type.
    # 0: 0, I: integer, F: float
    def check_type(self):
        for c in self.coeff:
            if c == 0:
                self.coeff_type.append("0")
            elif c.is_integer():
                self.coeff_type.append("I")
            else:
                self.coeff_type.append("F")
                
    # Return two Poly instance by seperating polynomial.
    def seperate(self, i: int, multA=False) -> tuple["Poly", "Poly"]:
        if self.is_type_poly:
            return self.seperate_cases(i, multA, True)[0]

        def trim(coeff: list[float]) -> list[float]:
            while coeff and coeff[-1] == 0:
                coeff.pop()
            return coeff
        
        coeff_p = trim(self.coeff[i:])
        coeff_q = trim(self.coeff[:i])
        
        if multA and coeff_p:
            leading_coeff = coeff_p[-1]
            if leading_coeff != 0:
                coeff_p = [c / leading_coeff for c in coeff_p]
                
        return Poly(coeff_p), Poly(coeff_q)

    def seperate_cases(self, i: int, multA=False, assumption1: bool=True) -> list[tuple["Poly", "Poly"]]:
        '''
        기능:
            다항식을 x^i 기준으로 분리하고, 타입 다항식의 ax^i 분리에서
            발생 가능한 quotient 계수 타입을 생성한다.
        입력:
            i: 분리에 사용할 x의 차수.
            multA: True이면 quotient의 최고차항 계수 a를 ax^i에 포함한다.
            assumption1: True이면 a로 나눈 모든 비영 하위 계수를 F로 둔다.
                         False이면 각 비영 하위 계수에 I/F를 모두 허용한다.
        출력:
            가능한 (quotient Poly, remainder Poly) 튜플의 목록.
            수치 다항식은 실제 나눗셈 결과 하나만 반환한다.
        '''
        def trim_type(coeff_type: list[str]) -> list[str]:
            '''
            기능:
                계수타입 목록의 후행 0을 제거한다.
            입력:
                coeff_type: 낮은 차수부터 정렬된 0/I/F 계수타입 목록.
            출력:
                후행 0이 제거된 계수타입 목록.
            '''
            while coeff_type and coeff_type[-1] == "0":
                coeff_type.pop()
            return coeff_type

        if not self.is_type_poly:
            poly_p, poly_q = self.seperate(i, multA)

            if multA and assumption1 and not poly_p.is_empty():
                poly_p.coeff_type = [
                    "0" if ctype == "0" else "I" if j == poly_p.deg else "F"
                    for j, ctype in enumerate(poly_p.coeff_type)
                ]

            return [(poly_p, poly_q)]

        coeff_type_p = trim_type(list(self.coeff_type[i:]))
        coeff_type_q = trim_type(list(self.coeff_type[:i]))
        poly_q = Poly(coeff_type_q)

        if not multA or not coeff_type_p:
            return [(Poly(coeff_type_p), poly_q)]

        if coeff_type_p[-1] != "F":
            return [(Poly(coeff_type_p), poly_q)]

        variable_indices = [
            j for j, ctype in enumerate(coeff_type_p[:-1])
            if ctype != "0"
        ]

        if assumption1:
            normalized_type = [
                "0" if ctype == "0" else "F"
                for ctype in coeff_type_p
            ]
            normalized_type[-1] = "I"
            return [(Poly(normalized_type), poly_q)]

        results = []
        seen = set()

        for type_case in product(["I", "F"], repeat=len(variable_indices)):
            normalized_type = list(coeff_type_p)

            for j, ctype in zip(variable_indices, type_case):
                normalized_type[j] = ctype

            normalized_type[-1] = "I"
            signature = tuple(normalized_type)

            if signature in seen:
                continue

            seen.add(signature)
            results.append((Poly(normalized_type), poly_q))

        return results
    
    # print
    def print(self, type="poly"):
        if type == "poly":
            print_poly(self.coeff)
        elif type == "type":
            print_poly_type(self.coeff_type)

    # isempty
    def is_empty(self) -> bool:
        return len(self.coeff) == 0
  
class XI:
    def __init__(self, multA: bool=False, n: int=0):
        self.multA = multA
        self.n = n
        self.made_powers = {0, 1}
        self.power_depths = {0: 0, 1: 0}
        # self.add_count = 0
        self.depth = 0
        self.pmult = 0
        self.cmult = 0
        self.route = []
        # (kind, lhs, rhs). kind in {"pure", "coeff", "coeff_direct"}
        self.route_ops = []
        # val = 1 if multA else 0
        # try:
        #     self.depth = ceil(log2(n + val))
        # except Exception as e:
        #     self.depth = 0
        # self.pmult = val
        
    def add_routes(self, route, depth, pmult, cmult, made_powers, route_ops=None, power_depths=None):
        '''
        기능:
            power 생성 route와 계산복잡도 및 출력 power-depth cache를 XI에 저장한다.
        입력:
            route: 기존 호환용 (lhs, rhs) 곱셈 목록.
            depth: 목표 x^n 또는 ax^n의 실제 multiplicative depth.
            pmult: plaintext multiplication 수.
            cmult: ciphertext multiplication 수.
            made_powers: 출력 cache에 존재하는 순수 power 차수 집합.
            route_ops: kind를 포함한 상세 곱셈 목록.
            power_depths: 각 순수 power의 실제 생성 depth dictionary.
        출력:
            없음. 현재 XI 객체를 직접 갱신한다.
        '''
        self.route = [tuple(r) for r in route]
        self.depth = depth
        self.pmult = pmult
        self.cmult = cmult
        if power_depths is None:
            self.power_depths = {
                int(power): ceil(log2(power)) if power > 1 else 0
                for power in ({0, 1} | set(made_powers))
            }
        else:
            self.power_depths = {
                int(power): int(power_depth)
                for power, power_depth in power_depths.items()
            }
            self.power_depths.setdefault(0, 0)
            self.power_depths.setdefault(1, 0)
        self.made_powers = set(self.power_depths)

        if route_ops is None:
            # 구버전 호출과의 호환성. 명시 정보가 없으면 순수 power route로 간주한다.
            self.route_ops = [("pure", int(a), int(b)) for a, b in self.route]
        else:
            self.route_ops = [
                (str(kind), int(a), int(b))
                for kind, a, b in route_ops
            ]
    
    def print_params(self):
        print(f"XI(n={self.n}, multA={self.multA}, DCP={self.depth}/{self.cmult}/{self.pmult}, route={self.route})")


class Decomp:
    def __init__(self, coeff: list[float], comp: Complexity, xi=None):
        self.coeff = coeff
        self.coeff_type = Poly(coeff).coeff_type if coeff else []
        self.comp = comp
        self.xi = xi if xi is not None else XI()
        self.dcmp_p = None
        self.dcmp_q = None
        self.made_powers = set(self.xi.made_powers)
        self.power_depths = dict(self.xi.power_depths)
        # terminal, p_only, q_only, p_then_q, q_then_p
        self.eval_order = "terminal"
        # 말단 다항식의 항별 평가계획. 각 항은 degree/coeff_type/multA/route_ops를 가진다.
        self.term_plans = []

    def update(self, xi: XI, dcmp_p, dcmp_q):
        self.xi = xi
        self.dcmp_p = dcmp_p
        self.dcmp_q = dcmp_q
        self.made_powers = set(xi.made_powers)
        self.power_depths = dict(xi.power_depths)
        
    def __lt__(self, other):
        # return (self.comp, int(self.xi.multA), self.check_depth()) < (other.comp, int(other.xi.multA), other.check_depth())
        return (self.comp, int(self.xi.multA), self.xi.n) < (other.comp, int(other.xi.multA), other.xi.n)

    def __eq__(self, other):
        # return (self.comp, self.xi.multA, self.check_depth()) == (other.comp, other.xi.multA, other.check_depth())
        return (self.comp, self.xi.multA) == (other.comp, other.xi.multA)
    
    def is_empty(self) -> bool:
        return self.coeff == []
    
    def restore_dcmp(self) -> str:               
        # (a)x^i
        res = ""
        if self.xi.n != 0:
            if self.xi.multA:
                res += f"{self.coeff[-1]}"
            res += "x"
            if self.xi.n >= 2:
                res += f"^{self.xi.n}"
        else:
            return poly_to_str(self.coeff)
        
        # p(x)
        if self.dcmp_p is not None:
            res += f"({self.dcmp_p.restore_dcmp()})"
            
        # q(x)
        if self.dcmp_q is not None:
            res += f" + ({self.dcmp_q.restore_dcmp()})"

        return res


    def merge_mp(self):
        mp2 = self.dcmp_p.merge_mp() if self.dcmp_p is not None else {0, 1}
        mp3 = self.dcmp_q.merge_mp() if self.dcmp_q is not None else {0, 1}
        return self.xi.made_powers | mp2 | mp3

    def merge_power_depths(self):
        '''
        기능:
            현재 decomposition 전체에서 생성되거나 재사용된 순수 power의 실제
            depth를 병합한다. 같은 power가 여러 경로에 있으면 더 작은 depth를 유지한다.
        입력:
            없음.
        출력:
            power 차수를 key, 실제 생성 depth를 value로 갖는 dictionary.
        '''
        result = dict(self.xi.power_depths)
        for child in (self.dcmp_p, self.dcmp_q):
            if child is None:
                continue
            for power, power_depth in child.merge_power_depths().items():
                result[power] = min(
                    result.get(power, power_depth),
                    power_depth,
                )
        result.setdefault(0, 0)
        result.setdefault(1, 0)
        return result
    
    def merge_route(self):
        routes = set(self.xi.route) if self.xi else set()

        for child in (self.dcmp_p, self.dcmp_q):
            if child is not None:
                routes |= child.merge_route()

        return routes
        
    def check_floats(self) -> int:
        res = 0
        for c in self.coeff:
            if c == 0:
                continue
            elif c.is_integer():
                continue
            else:
                res += 1
        if self.dcmp_p is not None:
            res += self.dcmp_p.check_floats()
        if self.dcmp_q is not None:
            res += self.dcmp_q.check_floats()

        return res
    
    # def check_depth(self) -> int:
    #     res = 1
    #     if self.xi.n != 0:
    #         res += ceil(log2(self.xi.n+1))

    #     pDepth = 0
    #     qDepth = 0
    #     if self.dcmp_p is not None:
    #         pDepth = self.dcmp_p.check_depth()
    #     if self.dcmp_q is not None:
    #         qDepth = self.dcmp_q.check_depth()
            
    #     return res+max(pDepth, qDepth)
    def check_depth(self) -> int:
        res = 1
        dp = 0
        dq = 0
        if self.dcmp_p:
            dp = self.dcmp_p.check_depth()
        if self.dcmp_q:
            dq = self.dcmp_q.check_depth()
        return res + max(dp, dq)
            
            
NULL_DECOMP = Decomp([], Complexity(), XI())
    
    
        
    
