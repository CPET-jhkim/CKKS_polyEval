# class.py
from print import *
from math import log2, ceil
# from util import *

class Complexity:
    '''
    연산복잡도 클래스
    '''
    def __init__(self):
        self.depth = 0
        self.cmult = 0
        self.pmult = 0
        self.add = 0
    
    def insert_value(self, depth, cmult, pmult, add):
        self.depth = depth
        self.cmult = cmult
        self.pmult = pmult
        self.add = add
    
    # 크기 비교 - 작으면 연산복잡도가 더 높음.
    def __lt__(self, other):
        return (self.depth, self.cmult, self.pmult, self.add) < (other.depth, other.cmult, other.pmult, other.add)
    
    def __eq__(self, other):
        return (self.depth, self.cmult, self.pmult, self.add) == (other.depth, other.cmult, other.pmult, other.add)
    
    def print_params(self):
        print(f"{'DCPA:':<8}{self.depth}|{self.cmult}|{self.pmult}|{self.add}")
        
class Poly:
    '''
    다항식 정보 클래스
    coeff: 계수 정보
    complexity: 연산복잡도 정보
    mp: 생성한 x^i차수 정보
    '''
    def __init__(self, coeff: list[float]):
        self.coeff = coeff
        self.deg = len(coeff) - 1
        self.coeff_type: list[str] = []
        self.check_type()
        self.complexity = Complexity()
        # self.mp: set[int] = set([0, 1])
        # self.ops_list = None
    
    # 다항식 각 계수의 타입 검사.
    # 0: 0, I: 정수, F: 소수
    def check_type(self):
        for c in self.coeff:
            if c == 0:
                self.coeff_type.append("0")
            elif c.is_integer():
                self.coeff_type.append("I")
            else:
                self.coeff_type.append("F")
                
    # 다항식 분해 -> 2개의 poly 클래스 반환.
    def seperate(self, i: int, multA=False) -> tuple["Poly", "Poly"]:
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
    
    # 데이터 출력용
    def print(self, type="poly"):
        if type == "poly":
            print_poly(self.coeff)
        elif type == "type":
            print_poly_type(self.coeff_type)

  
class XI:
    def __init__(self, multA: bool=False, n: int=0):
        self.multA = multA
        self.n = n
        self.made_powers = {0, 1}
        self.route = []
        val = 1 if multA else 0
        try:
            self.depth = ceil(log2(n + val))
        except Exception as e:
            self.depth = 0
        self.pmult = val
        
    def add_routes(self, route, made_powers):
        self.route = route
        self.add_count = len(route)
        self.made_powers = made_powers
    
    def print_params(self):
        print(f"XI(n={self.n}, multA={self.multA}, count={self.add_count}, route={self.route})")


class Decomp:
    def __init__(self, coeff: list[float], comp: Complexity, xi=XI()):
        self.coeff = coeff
        self.comp = comp
        self.xi = xi
        self.dcmp_p = None
        self.dcmp_q = None
        self.made_powers = self.xi.made_powers
        
    def update(self, xi: XI, dcmp_p: "Decomp", dcmp_q: "Decomp"):
        self.xi = xi
        self.dcmp_p = dcmp_p
        self.dcmp_q = dcmp_q
        self.made_powers |= xi.made_powers
        
    # 크기 비교
    def __lt__(self, other):
        return (self.comp, int(self.xi.multA), self.check_depth()) < (other.comp, int(other.xi.multA), other.check_depth())

    def __eq__(self, other):
        return (self.comp, self.xi.multA, self.check_depth()) == (other.comp, other.xi.multA, other.check_depth())
        
    def restore_dcmp(self) -> str:
        ###########        
        # poly가 비어있는 경우
        if self.coeff == []:
            return ""
        # i=0
        if self.xi.n == 0:
            return pp(self.coeff)
        else:
            coeff_p = self.dcmp_p.coeff
            coeff_q = self.dcmp_q.coeff
            if self.xi.multA and coeff_p:
                leading_coeff = coeff_p[-1]
                if leading_coeff != 0:
                    coeff_p = [c / leading_coeff for c in coeff_p]
            coeff = f"{self.coeff[-1]}" if self.xi.multA else ""
            if self.xi.n == 1:
                res = f"{coeff}x"
            else:
                res = f"{coeff}x^{self.xi.n}"
            
            q = f" + ({self.dcmp_q.restore_dcmp()})" if coeff_q != [] else ""
            return f"({res})[ {self.dcmp_p.restore_dcmp()} ]{q}" 

    def merge_mp(self):
        mp2 = self.dcmp_p.merge_mp() if self.dcmp_p is not None else {0, 1}
        mp3 = self.dcmp_q.merge_mp() if self.dcmp_q is not None else {0, 1}
        return self.xi.made_powers | mp2 | mp3
    
    def merge_route(self):
        return set(self.xi.route) | set(self.dcmp_p.xi.route) | set(self.dcmp_q.xi.route)
        
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
    
    def check_depth(self) -> int:
        res = 1
        if self.xi.n != 0:
            res += ceil(log2(self.xi.n+1))

        pDepth = 0
        qDepth = 0
        if self.dcmp_p is not None:
            pDepth = self.dcmp_p.check_depth()
        if self.dcmp_q is not None:
            qDepth = self.dcmp_q.check_depth()
            
        return res+max(pDepth, qDepth)
    

