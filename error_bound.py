from decimal import Decimal, getcontext
from basic_class import Decomp

getcontext().prec = 128

# ---------- helpers ----------
def D(v) -> Decimal:
    """안전한 Decimal 변환: float는 str 경유(이진 float 오차를 그대로 끌고오지 않게)"""
    if isinstance(v, Decimal):
        return v
    if isinstance(v, int):
        return Decimal(v)
    # float, str 등
    return Decimal(str(v))

def sqrtD(v) -> Decimal:
    return D(v).sqrt()


# ---------- EB ----------
class EB:
    def __init__(self, sigma: float, N: int, h: int, s: int):
        # 전부 Decimal로 보관
        self.Bc = B_clean(sigma, N, h)   # Decimal
        self.Bs = B_scale(N, h)          # Decimal
        self.scale = Decimal(2) ** int(s)

        self.B2 = Decimal(0)
        self.B4 = Decimal(0)
        self.B6 = Decimal(0)


# ---------- bounds ----------
def cal_bound(eb: EB, x: float, dcmp: Decomp) -> Decimal:
    """
    Error bound calculation (Decimal 버전)
    """
    xd = D(x)

    # 1) 모든 x^i 조합의 error bound 연산
    xi_eb: dict[int, Decimal] = {0: Decimal(0), 1: Decimal(0)}
    for xi1, xi2 in dcmp.merge_route():
        # x^xi1, x^xi2의 에러를 곱셈으로 합성
        xi_eb[xi1 + xi2] = eb_attach(eb, xd, xd, xi_eb[xi1], xi_eb[xi2], 'x')[1]

    # 2) 실제 분해식의 error bound 연산(재귀)
    res = cal_dcmp_bound(dcmp, eb, xd, xi_eb)
    return res


def cal_dcmp_bound(dcmp: Decomp, eb: EB, x: Decimal, xi_eb: dict[int, Decimal]) -> Decimal:
    """
    Error bound of (x^i)p(x) + q(x). (Decimal 버전)
    """
    # 0) 상수 다항식
    if len(dcmp.coeff) == 1:
        return Decimal(0)

    # 1) x^i
    n = int(dcmp.xi.n)
    pt_xi = x ** n
    eb_xi = xi_eb[n]

    # 2) p(x)
    if dcmp.dcmp_p is not None:
        poly_p = dcmp.dcmp_p.coeff
        pt_px = evalP(x, poly_p)
        eb_px = cal_dcmp_bound(dcmp.dcmp_p, eb, x, xi_eb)
    else:
        pt_px, eb_px = Decimal(0), Decimal(0)

    # 3) q(x)
    if dcmp.dcmp_q is not None:
        poly_q = dcmp.dcmp_q.coeff
        pt_qx = evalP(x, poly_q)
        eb_qx = cal_dcmp_bound(dcmp.dcmp_q, eb, x, xi_eb)
    else:
        pt_qx, eb_qx = Decimal(0), Decimal(0)

    # 4) attach
    pt_xip, eb_xip = eb_attach(eb, pt_xi, pt_px, eb_xi, eb_px, 'x')
    pt_xipq, eb_xipq = eb_attach(eb, pt_xip, pt_qx, eb_xip, eb_qx, '+')

    # 원 코드 유지: 마지막에 scale로 정규화
    return eb_xipq / eb.scale


def evalP(x: Decimal, coeff: list) -> Decimal:
    """
    다항식 평가: sum c_i * x^i
    - coeff 원소가 float일 수 있으므로 Decimal로 변환
    - Horner로 구현(정밀도/속도 모두 유리)
    """
    # coeff[0] + coeff[1] x + ... 형태를 Horner로
    # 높은 차수부터 내려오기
    res = Decimal(0)
    for c in reversed(coeff):
        res = res * x + D(c)
    return res


def eb_attach(eb: EB, x1, x2, err1, err2, op: str) -> tuple[Decimal, Decimal]:
    """
    모든 입력/출력을 Decimal로 처리
    """
    bs = D(eb.Bs)
    scale = eb.scale

    x1d = D(x1)
    x2d = D(x2)
    e1d = D(err1)
    e2d = D(err2)

    dx1 = scale * x1d
    dx2 = scale * x2d

    if op == 'x':  # 곱셈
        val = x1d * x2d
        approx = ((dx1 + e1d) * (dx2 + e2d)) / scale + bs
        target = scale * val
        err = abs(approx - target)
        return val, err

    elif op == '+':
        val = x1d + x2d
        err = abs(e1d + e2d)
        return val, err

    else:
        return Decimal(-99), Decimal(-99)


# ---------- B bounds (Decimal) ----------
def B_clean(sigma: float, N: int, h: int) -> Decimal:
    sd = D(sigma)
    Nd = D(N)
    hd = D(h)

    # (8*sqrt(2)*sigma*N) + (6*sigma*sqrt(N)) + (16*sigma*sqrt(h*N))
    term1 = D(8) * sqrtD(2) * sd * Nd
    term2 = D(6) * sd * sqrtD(Nd)
    term3 = D(16) * sd * sqrtD(hd * Nd)
    return term1 + term2 + term3


def B_scale(N: int, h: int) -> Decimal:
    Nd = D(N)
    hd = D(h)

    # sqrt(N/3) * (3 + 8*sqrt(h))
    return sqrtD(Nd / D(3)) * (D(3) + D(8) * sqrtD(hd))
