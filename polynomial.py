# polynomial.py
from complexity import Complexity
from collections import deque
from print import *
from util import *

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
        self.mp: set[int] = set([0, 1])
        self.ops_list = None
    
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
            

class Decomp:
    def __init__(self, coeff: list[float], comp: Complexity, mp: set[int]):
        self.coeff = coeff
        self.multA = False
        self.compA = 1
        self.i = 0
        self.comp = comp
        self.mp = mp
        self.dcmp_p = None
        self.dcmp_q = None

    def update(self, multA: bool, i: int, dcmp_p: "Decomp", dcmp_q: "Decomp"):
        self.multA = multA
        self.compA = 1-multA
        self.i = i
        self.dcmp_p = dcmp_p
        self.dcmp_q = dcmp_q
    
    # 크기 비교
    def __lt__(self, other):
        return (self.comp, int(self.multA), self.check_depth()) < (other.comp, int(other.multA), other.check_depth())

    def __eq__(self, other):
        return (self.comp, self.compA, self.check_depth()) == (other.comp, other.compA, other.check_depth())
        
    def restore_dcmp(self) -> str:
        ###########        
        # poly가 비어있는 경우
        if self.coeff == []:
            return ""
        # i=0
        if self.i == 0:
            return pp(self.coeff)
        else:
            coeff_p = self.dcmp_p.coeff
            coeff_q = self.dcmp_q.coeff
            if self.multA and coeff_p:
                leading_coeff = coeff_p[-1]
                if leading_coeff != 0:
                    coeff_p = [c / leading_coeff for c in coeff_p]
            coeff = f"{self.coeff[-1]}" if self.multA else ""
            if self.i == 1:
                res = f"{coeff}x"
            else:
                res = f"{coeff}x^{self.i}"
            
            q = f" + ({self.dcmp_q.restore_dcmp()})" if coeff_q != [] else ""
            return f"({res})[ {self.dcmp_p.restore_dcmp()} ]{q}" 

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
        if self.i != 0:
            res += ceil(log2(self.i+1))

        pDepth = 0
        qDepth = 0
        if self.dcmp_p is not None:
            pDepth = self.dcmp_p.check_depth()
        if self.dcmp_q is not None:
            qDepth = self.dcmp_q.check_depth()
            
        return res+max(pDepth, qDepth)
