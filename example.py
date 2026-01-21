# example.py
from algorithm import cal_polyEval
from util import *
from basic_class import Poly, Decomp
from error_bound import EB, cal_bound

if __name__ == '__main__': 
    # Define coefficient, error bound class(EB), input x
    coeff = [1.2, 3.4, 3.6]
    poly = Poly(coeff)
    eb = EB(sigma=3.1, N=pow(2,17), h=192, s=50) 
    x = 1024
    poly.print()
    
    # Calculate optimal decomposition.
    result: Decomp = cal_polyEval(poly)  
    print(f"Decomposition result: {result.restore_dcmp()}")
    
    # Calculate error bound using calculated decomposition above.  
    error = cal_bound(eb, x, result)
    print(f"Error bound: {error:.16f}")