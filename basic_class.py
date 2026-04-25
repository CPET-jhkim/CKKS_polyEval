# basic_class.py
try:
    from print import poly_to_str, print_poly, print_poly_type
except:
    from .print import poly_to_str, print_poly, print_poly_type
from math import log2, ceil

class Complexity:
    '''
    Calculation complexity class
    '''
    def __init__(self, xi = None):
        if xi is not None:
            self.depth = xi.depth
            self.cmult = xi.add_count
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
        
class Poly:
    '''
    Polynomial info class
    coeff: coefficient info
    complexity: calculation complexity info.
    mp: created x^i info.
    '''
    def __init__(self, coeff: list[float]):
        self.coeff = [float(cf) for cf in coeff]
        self.deg = max(len(coeff) - 1, 0)
        self.coeff_type: list[str] = []
        self.check_type()
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
        self.add_count = 0
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
        
    def update(self, xi: XI, dcmp_p, dcmp_q):
        self.xi = xi
        self.dcmp_p = dcmp_p
        self.dcmp_q = dcmp_q
        self.made_powers |= xi.made_powers
        
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
    
    
        
    
