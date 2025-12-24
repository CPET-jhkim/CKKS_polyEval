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
    '''
    def __init__(self, coeff: list[float]):
        self.coeff = coeff
        self.deg = len(coeff) - 1
        self.coeff_type: list[str] = []
        self.check_type()
        self.complexity = Complexity()

    
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
    def seperate(self, i: int, divide=False) -> tuple["Poly", "Poly"]:
        def trim(coeff: list[float]) -> list[float]:
            while coeff and coeff[-1] == 0:
                coeff.pop()
            return coeff
        
        coeff_p = trim(self.coeff[i:])
        coeff_q = trim(self.coeff[:i])
        
        if divide and coeff_p:
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