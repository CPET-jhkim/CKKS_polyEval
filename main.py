# main.py
from algorithm import calculate
from contextlib import redirect_stdout
from print import print_poly, print_poly_type, decomp_poly
from util import make_all_polys
from polynomial import Poly

if __name__ == '__main__':
    made_powers: set[int] = {0, 1}
    coeff = [8, 2.6, 1.2, 1]
    poly = Poly(coeff)
    res = calculate(poly, made_powers)
    res[0].print_params()
    dd = decomp_poly(coeff, res[-1])
    print(dd)
    
    filename = "deg.txt"
    
    with open(filename, 'w', encoding='utf-8') as f:
        with redirect_stdout(f):
            for deg in range(2, 6):
                print('#' * 20)
                print(f"Degree: {deg}")
                print('#' * 20)

                polys = make_all_polys(deg)
                for coeff in polys:
                    poly = Poly(coeff)
                    poly.print("poly")
                    print()
                    poly.print("type")
                    result = calculate(poly, made_powers)
                    print(f"분해식:\t{decomp_poly(coeff, result[-1])}")
                    result[0].print_params()
                    
                    print('-'*20)
