from algorithm import cal_polyEval
from contextlib import redirect_stdout
from basic_class import Poly
from print import decomp_poly
from util import make_all_polys
from deprecated.error_bound import ErrBound, EvalStep, evaluate_polynomial_dcmp

if __name__ == '__main__':
    filename = "deg.txt"
    poly = [0.1, 0, 0.2, 0, 0.3, 0, 0.4]
    result = cal_polyEval(poly)
    print(f"{'dcmp:':<8}{result.restore_dcmp()}")
    a = EvalStep(result)
    a.print_step()
    
    fct = evaluate_polynomial_dcmp(es=a, x=[0, 1, 2, 3, 4], eb=ErrBound(sigma=3.1, N=pow(2, 17), h=192, s=50))
    print(fct.ct_high)
    print(fct.pt)
    print(fct.ct_low)
    
    # with open(filename, 'w', encoding='utf-8') as f:
    #     with redirect_stdout(f):
    #         for deg in range(16, 17):
    #             print('#' * 20)
    #             print(f"Degree: {deg}")
    #             print('#' * 20)

    #             polys = make_all_polys(deg)
    #             for coeff in polys:
    #                 result = cal_polyEval(coeff)
    #                 # if type(result) == bool:
    #                 poly = Poly(coeff)
    #                 poly.print("poly")
    #                 poly.print("type")
    #                 result = cal_polyEval(coeff)
    #                 print(f"{'dcmp:':<8}{result.restore_dcmp()}")
    #                 poly.print("type")
    #                 result.comp.print_params()

    #                 print('-'*50)