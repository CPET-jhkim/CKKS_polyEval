# all_n_degree_polys.py
try:
    from poly_util import make_type_poly
    from algorithm import cal_polyEval
    from basic_class import Poly, Decomp
except Exception as e:
    from .poly_util import make_type_poly
    from .algorithm import cal_polyEval
    from .basic_class import Poly, Decomp

if __name__ == '__main__': 
    # make arbitrary polynomial by defining each coeff type.
    # "F": float, "0": zero, "I": integer
    coeff = make_type_poly(['F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F'])
    poly = Poly(coeff)
    poly.print()
    
    # Calculate optimal decomposition.
    result: Decomp = cal_polyEval(poly)  
    print(f"Decomp result:\t{result.restore_dcmp()}")
    result.comp.print_params()


        
