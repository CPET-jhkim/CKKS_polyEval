# main.py
from algorithm import calculate
from contextlib import redirect_stdout
from print import decomp_poly
from util import make_all_polys
from polynomial import Poly

if __name__ == '__main__':
    made_powers: set[int] = {0, 1}    
    # coeff = [0, 0, 3, 6.9]
    # c = Poly(coeff)
    # result = calculate_step1(c, made_powers)
    # print(result[-1].restore_dcmp())
    
    
    filename = "deg.txt"
    with open(filename, 'w', encoding='utf-8') as f:
        with redirect_stdout(f):
            for deg in range(3, 4):
                print('#' * 20)
                print(f"Degree: {deg}")
                print('#' * 20)

                polys = make_all_polys(deg)
                for coeff in polys:
                    poly = Poly(coeff)
                    poly.print("poly")
                    result = calculate(poly, made_powers)
                    
                    print('-'*50)
