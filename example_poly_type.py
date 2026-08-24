# all_n_degree_polys.py
try:
    from poly_util import make_type_poly, make_deg_poly, print_dcmp_detail
    from algorithm import cal_polyEval
    from basic_class import Poly, Decomp
except Exception as e:
    from .poly_util import make_type_poly, make_deg_poly, print_dcmp_detail
    from .algorithm import cal_polyEval
    from .basic_class import Poly, Decomp
        
        
if __name__ == '__main__': 
    # make arbitrary polynomial by defining each coeff type.
    # "F": float, "0": zero, "I": integer
    coeff = make_type_poly(['0', 'F','0', 'F','0', 'F', '0', 'F'])
    # deg = 7
    # coeff = make_deg_poly(deg, "all")
    poly = Poly(coeff)
    poly.print()
    
    # Calculate optimal decomposition.
    result: Decomp = cal_polyEval(poly)  
    print_dcmp_detail(result, 0)


