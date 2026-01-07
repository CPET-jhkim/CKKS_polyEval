# main.py
from algorithm import calculate
from contextlib import redirect_stdout
from print import decomp_poly
from util import make_all_polys
from polynomial import Poly

if __name__ == '__main__':
    made_powers: set[int] = {0, 1}
    # coeff = [2, 0, 0, 0.6, 0, 0, 1]
    # poly = Poly(coeff)
    # poly.print("poly")
    # result = calculate(poly, made_powers)
    # print(f"dcmp:\t{decomp_poly(coeff, result[-1])}")
    # poly.print("type")
    # result[0].print_params()
    
    filename = "deg.txt"
    
    with open(filename, 'w', encoding='utf-8') as f:
        with redirect_stdout(f):
            for deg in range(5, 6):
                print('#' * 20)
                print(f"Degree: {deg}")
                print('#' * 20)

                polys = make_all_polys(deg)
                for coeff in polys:
                    poly = Poly(coeff)
                    poly.print("poly")
                    result = calculate(poly, made_powers)
                    print(f"{'dcmp:':<8}{result[-1].restore_dcmp()}")
                    poly.print("type")
                    result[0].print_params()
                    
                    print('-'*50)
