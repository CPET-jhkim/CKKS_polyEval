try:
    from math import sqrt, log2
    import numpy as np
    from dataclasses import dataclass
    from collections import deque
    from basic_class import Decomp
except:
    from ..basic_class import Decomp
    
import sys
sys.path.append("./build")
from EBModule import *


""" class ErrBound:
    def __init__(self, sigma: float, N: int, h: int, s: int):
        self.Bc = B_clean(sigma, N, h)
        self.Bs = B_scale(N, h)
        self.scale = pow(2, s)

def B_clean(sigma: float, N: int, h: int) -> float:
    # (8*sqrt(2)*sigma*N) + (6*sigma*sqrt(N)) + (16*sigma*sqrt(h*N))
    term1 = 8 * sqrt(2) * sigma * N
    term2 = 6 * sigma * sqrt(N)
    term3 = 16 * sigma * sqrt(h * N)
    return term1 + term2 + term3
    # return ((8*sqrt(2))*sigma*N) + (6*sigma*sqrt(N)) + (16*sigma*sqrt(h*N))

def B_scale(N: int, h: int) -> float:
    # sqrt(N/3) * (3 + 8*sqrt(h))
    return sqrt(N / 3) * (3 + 8 * sqrt(h))
    # return sqrt(3*N) + (8*sqrt((h*N)/3)) """

""" @dataclass
class Step:
    op: str
    key1: str
    key2: str
    save_key: str
    
class EvalStep:
    def __init__(self, dcmp: Decomp | None = None):
        self.term_count = 0
        self.coeff_count = 0
        self.powers: dict[str, str] = {}
        self.coeffs: dict[str, str] = {}
        self.terms: dict[str, str] = {}
        self.eval_step: list[Step] = []

        if dcmp is not None:
            self.add_value('P', 1, "x")
            
            # 모든 생성된 route와 power 가져오기
            made_powers = dcmp.merge_mp()
            routes = list(dcmp.merge_route())
            # 작은 차수부터 먼저 연산될 수 있도록 합(target power) 기준으로 정렬
            routes.sort(key=lambda r: r[0] + r[1])
            
            # x^i 구성
            for route in routes:
                target_power = route[0] + route[1]
                if target_power in made_powers:
                    key1 = f"P{route[0]}"
                    key2 = f"P{route[1]}"
                    self.add_step(key1, key2, 'x')
                    
            self.serialize_dcmp(dcmp)

    def serialize_dcmp(self, dcmp: Decomp) -> str:
        xi = dcmp.xi
        res_key = ""

        # 1. Leaf polynomial
        if dcmp.dcmp_p is None and dcmp.dcmp_q is None:
            for i in range(len(dcmp.coeff) - 1, -1, -1):
                if i == 0 and dcmp.coeff[0] != 0:
                    coeff_key = self.add_value('C', self.coeff_count, f"{dcmp.coeff[0]}")
                    self.coeff_count += 1
                    res_key = coeff_key if res_key == "" else self.add_step(res_key, coeff_key, '+')
                
                if i != 0 and dcmp.coeff[i] != 0:
                    if dcmp.coeff[i] != 1:
                        coeff_key = self.add_value('C', self.coeff_count, f"{dcmp.coeff[i]}")
                        self.coeff_count += 1
                        xi_key = f"P{i}"
                        axi_key = self.add_step(coeff_key, xi_key, 'x')
                        res_key = axi_key if res_key == "" else self.add_step(res_key, axi_key, '+')
                    else:
                        xi_key = f"P{i}"
                        res_key = xi_key if res_key == "" else self.add_step(res_key, xi_key, '+')
            return res_key

        # 2. multA 여부에 따른 연산 순서 등록
        made_powers = dcmp.merge_mp()
        
        if xi.multA:
            xi_coeff = dcmp.coeff[-1] 
            coeff_key = self.add_value('C', self.coeff_count, f"{xi_coeff}")
            self.coeff_count += 1
            
            if xi.n in made_powers:
                # ax^n인 경우
                xi_key = f"P{xi.n}"
                axi_key = self.add_step(coeff_key, xi_key, 'x')
            else:
                # ax^i * x^j인 경우
                xi_key = f"P{xi.route[-1][0]}"
                axi_key = self.add_step(coeff_key, xi_key, 'x')
                
                xj_key = f"P{xi.route[-1][1]}"
                axi_key = self.add_step(axi_key, xj_key, 'x')
        else:
            axi_key = f"P{xi.n}"

        axipx_key = axi_key
        if dcmp.dcmp_p is not None:
            px_key = self.serialize_dcmp(dcmp.dcmp_p)
            axipx_key = self.add_step(axi_key, px_key, 'x')

        res_key = axipx_key
        if dcmp.dcmp_q is not None:
            qx_key = self.serialize_dcmp(dcmp.dcmp_q)
            res_key = self.add_step(axipx_key, qx_key, '+')

        return res_key

    def get_value(self, key: str) -> str:
        prefix = key[0]
        if prefix == 'P':
            return self.powers[key]
        elif prefix == 'C':
            return self.coeffs[key]
        elif prefix == 'T':
            return self.terms[key]
        else:
            raise ValueError(f"EvalStep.get_value: invalid key prefix: {key}")

    def add_value(self, key_header: str, key_number: int, data: str) -> str:
        key = f"{key_header}{key_number}"
        
        if key_header == 'C':
            self.coeffs[key] = data
            self.eval_step.append(Step('o', key, data, ""))
        elif key_header == 'P':
            self.powers[key] = data
            if key_number == 1:
                self.eval_step.append(Step('o', key, data, ""))
        elif key_header == 'T':
            self.terms[key] = data
            
        return key

    def add_step(self, key1: str, key2: str, op: str) -> str:
        if op == '+':
            key_header = 'T'
            term_number = self.term_count
            self.term_count += 1
            data = f"({self.get_value(key1)} + {self.get_value(key2)})"
            
        elif key1[0] == 'P' and key2[0] == 'P':
            key_header = 'P'
            term_number = int(key1[1:]) + int(key2[1:]) 
            data = f"x^{term_number}"
            
        elif key1[0] == 'C' and key2[0] == 'P':
            key_header = 'T'
            term_number = self.term_count
            self.term_count += 1
            data = self.get_value(key1) + self.get_value(key2)
            
        elif key1[0] == 'P' and key2[0] == 'C':
            key_header = 'T'
            term_number = self.term_count
            self.term_count += 1
            data = self.get_value(key1) + self.get_value(key2)
            
        elif key1[0] == 'P' and key2[0] == 'T':
            key_header = 'T'
            term_number = self.term_count
            self.term_count += 1
            data = self.get_value(key1) + self.get_value(key2)
            
        elif key1[0] == 'T' and key2[0] == 'T':
            key_header = 'T'
            term_number = self.term_count
            self.term_count += 1
            data = self.get_value(key1) + self.get_value(key2)
            
        else:
            raise ValueError(f"EvalStep.add_step: invalid key data {key1} {key2} {op}")
        
        save_key = self.add_value(key_header, term_number, data)
        self.eval_step.append(Step(op, key1, key2, save_key))
        return save_key

    def print_step(self):
        count = 1
        for s in self.eval_step:
            self.print_line(s, count)
            count += 1

    def print_line(self, s: Step, count: int):
        if s.op == '+':
            op_string = "ADD"
        elif s.op == 'x':
            op_string = "MUL"
        elif s.op == 'o':
            op_string = "RES"
            
        print(f"{count}\t{op_string}\t{s.key1}\t{s.key2}\t{s.save_key}") """

""" class FakeCiphertext:
    def __init__(self, x: float | np.ndarray, size=None, err=0.0, scale=1.0, ct_high=None, ct_low=None):
        self.pt: np.ndarray = np.array([], dtype=np.float64)
        self.ct_high: np.ndarray = np.array([], dtype=np.float64)
        self.ct_low: np.ndarray = np.array([], dtype=np.float64)
        self.err: float = 0.0
        self.scale: float = 1.0

        # 생성자
        if x is None:
            pass
        
        # x가 단일 스칼라 값으로 제시된 경우
        elif isinstance(x, float) and size is not None:
            self.pt = np.full(size, x * scale, dtype=np.float64)
            self.set_error(err)
            self.scale = scale
            
        # x가 벡터 형식으로 주어진 경우
        elif isinstance(x, np.ndarray):
            self.pt = x
            # x_arr = np.array(x, dtype=np.float64)
            # self.pt = x_arr * scale
            self.scale = scale
            
            if ct_high is not None and ct_low is not None:
                self.ct_high = np.array(ct_high, dtype=np.float64)
                self.ct_low = np.array(ct_low, dtype=np.float64)
            else:
                self.set_error(err)

    def set_error(self, err):
        self.err = err
        self.ct_high: np.ndarray = self.pt + err
        self.ct_low: np.ndarray = self.pt - err

def add(fct1: FakeCiphertext, fct2: FakeCiphertext) -> FakeCiphertext:
        assert fct1.scale == fct2.scale
        
        pt = fct1.pt + fct2.pt
        ct_high = fct1.ct_high + fct2.ct_high
        ct_low = fct1.ct_low + fct2.ct_low
        
        return FakeCiphertext(x=pt, scale=fct1.scale, err=0.0, ct_high=ct_high, ct_low=ct_low)

def add_pt(fct: FakeCiphertext, scalar: float) -> FakeCiphertext:
    scalar *= fct.scale
    pt = fct.pt + scalar
    ct_high = fct.ct_high + scalar
    ct_low = fct.ct_low + scalar

    return FakeCiphertext(x=pt, err=0.0, scale=fct.scale, ct_high=ct_high, ct_low=ct_low)

def mult(fct1: FakeCiphertext, fct2: FakeCiphertext, eb: ErrBound) -> FakeCiphertext:
    pt = (fct1.pt * fct2.pt)
    err = eb.Bs
    scale = (fct1.scale * fct2.scale) / eb.scale
    
    bound1 = fct1.ct_high * fct2.ct_high
    # bound2 = fct1.ct_high * fct2.ct_low
    # bound3 = fct1.ct_low * fct2.ct_high
    bound4 = fct1.ct_low * fct2.ct_low

    bounds = np.stack([bound1, bound4])
    ct_high = np.max(bounds, axis=0)
    ct_low = np.min(bounds, axis=0)

    ct_high /= eb.scale
    ct_low /= eb.scale
    pt /= eb.scale

    ct_high = ct_high + err
    ct_low = ct_low - err

    return FakeCiphertext(x=pt, scale=scale, err=err, ct_high=ct_high, ct_low=ct_low)

def square(fct: FakeCiphertext, eb: ErrBound) -> FakeCiphertext:
    pt = fct.pt * fct.pt
    err = eb.Bs
    scale = pow(fct.scale, 2) / eb.scale
    
    bound1 = fct.ct_high * fct.ct_high
    bound4 = fct.ct_low * fct.ct_low

    bounds = np.stack([bound1, pt, bound4])
    ct_high = np.max(bounds, axis=0)
    ct_low = np.min(bounds, axis=0)

    ct_high /= eb.scale
    ct_low /= eb.scale
    pt /= eb.scale

    ct_high = ct_high + err
    ct_low = ct_low - err

    return FakeCiphertext(x=pt, scale=scale, err=err, ct_high=ct_high, ct_low=ct_low)

def mult_pt(fct: FakeCiphertext, scalar: float, eb: ErrBound, rescale=True) -> FakeCiphertext:
    err = 0.0
    scale = fct.scale
    if rescale:
        scalar *= fct.scale
        err = eb.Bs
    
    pt = fct.pt * scalar
    ct_high = fct.ct_high * scalar
    ct_low = fct.ct_low * scalar
    
    if scalar < 0:
        ct_high, ct_low = ct_low, ct_high
    
    if rescale:
        pt /= eb.scale
        ct_high /= eb.scale
        ct_high += err
        ct_low /= eb.scale
        ct_low -= err
        
    return FakeCiphertext(x=pt, err=err, scale=scale, ct_high=ct_high, ct_low=ct_low)

 """
""" def cal_bound(eb: ErrBound, x: float | list[float], dcmp: Decomp) -> float | list[float]:
    if type(dcmp) == bool:
        return float(99)
    
    try:
        # 1) Calculate every x^i error bound.
        x_arr = np.atleast_1d(x)
        zero_eb = np.zeros_like(x_arr)
        xi_eb = {0: zero_eb, 1: zero_eb}
        for xi1, xi2 in dcmp.merge_route():
            xi_eb[xi1 + xi2] = eb_attach(eb, x_arr, x_arr, xi_eb[xi1], xi_eb[xi2], 'x')[1]
    except Exception as e:
        print(e)

    try:
        if len(dcmp.coeff) == 7:
            pass
        # 2) Calculate f(x)'s error bound recursively.
        res = cal_dcmp_bound(dcmp, eb, x, xi_eb)
        return res.tolist()
    except Exception as e:
        if np.isscalar(x):
            return 99.0
        return [99.0] * len(x) """

def get_value(key: str, powers, coeffs, terms) -> Ciphertext:
    key_header = key[0]
    if key_header == 'P':
        return powers[key]
    elif key_header == 'C':
        return coeffs[key]
    elif key_header == 'T':
        return terms[key]
    else:
        print(f"get_value error! key={key}")
        return Ciphertext()

# EvalStep 클래스의 각 연산순서에 따라 fakeCiphertext를 연산하여 에러 상한/하한을 연산.
def evaluate_polynomial_dcmp(es: EvalStep, x: Ciphertext, eb: ErrBound) -> Ciphertext:
    powers: dict[str, Ciphertext] = {}
    coeffs: dict[str, Ciphertext] = {}
    terms: dict[str, Ciphertext] = {}
    
    size = len(x.pt)
    powers['P1'] = x
    
    for step in es.eval_step:
        if step.op == 'o':
            if step.key1[0] != 'C':
                if step.key1[0] == 'P':
                    continue
                raise ValueError(f"FRES: invalid key data {step.key1} {step.key2}")

            fct_coeff = Ciphertext(x=float(step.key2), size=size, err=eb.Bc, scale=eb.scale)
            coeffs[step.key1] = fct_coeff
            continue
            
        # 2. 덧셈 연산
        elif step.op == '+':            
            v1 = get_value(step.key1, powers, coeffs, terms)
            v2 = get_value(step.key2, powers, coeffs, terms)
            terms[step.save_key] = ct_add(v1, v2)
            
        # 3. 곱셈 연산
        elif step.op == 'x':
            v1 = get_value(step.key1, powers, coeffs, terms)
            v2 = get_value(step.key2, powers, coeffs, terms)
            key1_header = step.key1[0]
            key2_header = step.key2[0]
            
            if key1_header == 'P' and key2_header == 'P':
                if step.key1 == step.key2:
                    powers[step.save_key] = ct_square(v1, eb)
                else:
                    powers[step.save_key] = ct_mult(v1, v2, eb)
            else:
                terms[step.save_key] = ct_mult(v1, v2, eb)
            
    # 값 반환
    last_key = es.eval_step[-1].save_key
    return terms[last_key]

def evaluate_polynomial_cleanse(x: Ciphertext, eb: ErrBound) -> Ciphertext:
    # x^2 구성
    x2 = ct_square(x, eb)
    # compare_fct(x2)
    
    # -2x + 3
    tmp: Ciphertext = ct_mult_pt(-2, x, eb, False)
    tmp = ct_add_pt(tmp, 3)
    # compare_fct(tmp)
    
    res = ct_mult(x2, tmp, eb)
    # compare_fct(res)
    
    return res


def compare_fct(x: Ciphertext):
    pt = x.pt / x.scale
    ct1 = x.ct_high / x.scale
    ct2 = x.ct_low / x.scale
    
    print(f"LBound < {bit_diff(ct1, pt)} < pt < {bit_diff(ct2, pt)} < HBound")
    

def bit_diff(bound: np.ndarray, data: np.ndarray) -> str:
    bound = np.asarray(bound, dtype=float)
    data = np.asarray(data, dtype=float)

    if bound.shape != data.shape:
        raise ValueError("bound와 data의 shape가 같아야 합니다.")

    diff = np.abs(bound - data)
    min_idx = np.unravel_index(np.argmin(diff), diff.shape)

    min_diff = diff[min_idx]
    if min_diff == 0.0:
        return "0 bits"

    sign = '-' if bound[min_idx] > data[min_idx] else '+'
    log2_diff = abs(log2(min_diff))

    return f"{sign}{log2_diff:.4f} bits"

""" def cal_dcmp_bound(dcmp: Decomp, eb: ErrBound, x: float | list[float], xi_eb: dict[int, float]) -> float | np.ndarray:
    x = np.atleast_1d(x)
    
    # 0) max_deg = 0
    if len(dcmp.coeff) == 1:
        return float(0)

    # 1) x^i
    n = int(dcmp.xi.n)
    pt_xi = float(1) if n == 0 else x ** n
    eb_xi = xi_eb[n]

    # 2) p(x)
    if dcmp.dcmp_p is not None:
        poly_p = dcmp.dcmp_p.coeff
        pt_px = evalP(x, poly_p)
        eb_px = cal_dcmp_bound(dcmp.dcmp_p, eb, x, xi_eb)
    else:
        pt_px, eb_px = float(0), float(0)

    # 3) q(x)
    if dcmp.dcmp_q is not None:
        poly_q = dcmp.dcmp_q.coeff
        pt_qx = evalP(x, poly_q)
        eb_qx = cal_dcmp_bound(dcmp.dcmp_q, eb, x, xi_eb)
    else:
        pt_qx, eb_qx = float(0), float(0)

    # 4) attach
    pt_xip, eb_xip = eb_attach(eb, pt_xi, pt_px, eb_xi, eb_px, 'x')
    pt_xipq, eb_xipq = eb_attach(eb, pt_xip, pt_qx, eb_xip, eb_qx, '+')

    return eb_xipq / eb.scale """

""" def evalP(x: np.ndarray, coeff: list) -> np.ndarray:
    res = float(0)
    for c in reversed(coeff):
        res = res * x + c
    return res """

""" def eb_attach(eb: ErrBound, x1: float | list, x2: float | list, err1: float | list, err2: float | list, op: str) -> tuple:
    x1, x2, err1, err2 = map(np.atleast_1d, [x1, x2, err1, err2])
    
    bs = eb.Bs
    scale = eb.scale

    dx1 = scale * x1
    dx2 = scale * x2

    if op == 'x':  # 곱셈
        val = x1 * x2
        approx = ((dx1 + err1) * (dx2 + err2)) / scale + bs
        target = scale * val
        err = abs(approx - target)
        return val, err

    elif op == '+':
        val = x1 + x2
        err = abs(err1 + err2)
        return val, err

    else:
        return [-99], [-99] """

