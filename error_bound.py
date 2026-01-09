from math import sqrt, pow, log2
from polynomial import Decomp

class EB:
    def __init__(self, sigma: float, N: int, h: int, s: int):
        self.Bc = B_clean(sigma, N, h)
        self.Bs = B_scale(N, h)
        self.scale = pow(2, s)
        self.B2 = 0.0
        self.B4 = 0.0
        self.B6 = 0.0

# Error bound of coeff, x
def cal_bound(eb: EB, x: float, coeff: list, dcmp: Decomp) -> float:
    '''
    Error bound calculation
    1. Calculate all error bound of x^i from made_powers
    2. Attach ciphertext and calculate error bound from dcmp (Bottom-Up)
    '''
    max_deg = len(coeff) - 1
    
    # 1
    xi_eb = cal_xi_bound(eb, x, dcmp.mp)
    
    # 2
    
    
def cal_xi_bound(eb: EB, x: float, mp: set[int]) -> dict[int, float]:
    """
    덧셈 조합에 따라 x^i의 error bound 연산.
    """
    mp_sorted: list[int] = sorted(list(mp))
    res: dict[int, float] = {}

    for i in mp_sorted:
        if i <= 1:
            res[i] = 0.0
            continue
                    
        best_pair: tuple[int, int]
        
        start_index = 0
        candidates = [n for n in mp_sorted if n <= i / 2]
        
        found = False
        for x in reversed(candidates):
            y = i - x
            if y in mp:
                best_pair = (x, y)
                found = True
                break
        if found and best_pair:
            err1, err2 = res[best_pair[0]], res[best_pair[1]]
            x1, x2 = pow(x, best_pair[0]), pow(x, best_pair[1])
            val = eb_attach(eb, x1, x2, err1, err2, 'x')
            res[i] = val
        else:
            res[i] = -99

    return res

def cal_dcmp_bound(dcmp: Decomp, eb: EB, x: float, xi_eb: dict[int, float]) -> float:
    '''
    Error bound of (x^i)p(x) + q(x).
    '''
    # 0. 상수다항식
    if len(dcmp.coeff) == 1:
        return 0.0        
    
    # 1. x^i
    pt_xi = pow(x, dcmp.i)
    eb_xi = xi_eb[dcmp.i]
    
    # 2. p(x)
    pt_px = evalP(x, dcmp.dcmp_p.coeff)
    eb_px = cal_dcmp_bound(dcmp.dcmp_p, eb, x, xi_eb)
    
    # 3. q(x)
    pt_qx = evalP(x, dcmp.dcmp_q.coeff)
    eb_qx = cal_dcmp_bound(dcmp.dcmp_q, eb, x, xi_eb)
    
    # 4. attach
    pt_xip, eb_xip = eb_attach(eb, pt_xi, pt_px, eb_xi, eb_px, 'x')
    pt_xipq, eb_xipq = eb_attach(eb, pt_xip, pt_qx, eb_xip, eb_qx, '+')
    
    return eb_xipq
    

def evalP(x: float, coeff: list[float]) -> float:
    return sum(c*pow(x, i) for i, c in enumerate(coeff))

def eb_attach(eb: EB, x1: float, x2: float, err1: float, err2: float, op: str) -> tuple[float, float]:
    bs = eb.Bs
    if op == 'x': # 곱셈
        return x1*x2, abs((((x1+err1) * (x2+err2)) / eb.scale + bs) - (x1*x2))
    elif op == '+':
        return x1+x2, abs(err1 + err2)
    else:
        return -99, -99


    
# Bclean
def B_clean(sigma: float, N: int, h: int) -> float:
    return (8 * sqrt(2.0) * sigma * N) + (6 * sigma * sqrt(N)) + (16 * sigma * sqrt(h * N))

# Bscale
def B_scale(N: int, h: int) -> float:
    return sqrt(N / 3.0) * (3.0 + 8.0 * sqrt(h))


