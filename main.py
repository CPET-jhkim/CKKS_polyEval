# main.py
from algorithm import calculate
from contextlib import redirect_stdout
from print import decomp_poly
from util import *
from basic_class import Poly
from error_bound import EB, cal_bound, eb_attach

if __name__ == '__main__': 
    coeff = make_type_poly(["F", "0", "F", "F", "F", "0", "F", "0", "F"])
    eb = EB(sigma=3.1, N=pow(2,17), h=192, s=50) 
    x = 1024
    c = Poly(coeff)
    c.print()
    result = calculate(c, {0, 1})    

    # print("Error bound test")
    # print(f"x: {x}")
    # print(f"Error bound of poly: {cal_bound(eb, x, result):.10f}")
    
    # print("-"*30)
    # xi_eb: dict = {0: 0.0, 1: 0.0}
    # route = [(1, 1), (1, 2), (2, 2), (1, 4), (3, 3), (3, 4), (4, 4)]
    # for xi1, xi2 in route:
    #     xi_eb[xi1+xi2] = eb_attach(eb, x, x, xi_eb[xi1], xi_eb[xi2], 'x')[1]
    #     print(f"x^{xi1+xi2} eb: {xi_eb[xi1+xi2]:.10f}")
        
    # filename = "deg.txt"
    # for deg in range(3, 10):
    #     filename = f"deg{deg}.txt"
    #     with open(filename, 'w', encoding='utf-8') as f:
    #         with redirect_stdout(f):
    #             print('#' * 20)
    #             print(f"Degree: {deg}")
    #             print('#' * 20)

    #             polys = make_all_polys(deg)
    #             for coeff in polys:
    #                 poly = Poly(coeff)
    #                 poly.print("poly")
    #                 poly.print("type")
    #                 result = calculate(poly, made_powers)
                    
    #                 print('-'*50)
