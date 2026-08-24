# print.py

def print_poly(poly: list[float], title: str="poly:") -> None:
    print(f"{title:<8}{poly_to_str(poly)}")
                
def print_poly_type(poly_type: list[str], title: str="type:") -> None:
    print(f"{title:<8}({' '.join(reversed(poly_type))})")

def print_poly_sep(i: int, poly_p: list[float], poly_q: list[float]) -> None:
    print("dcmp: ", end='')
    print(f"(x^{i})", end="")
    print("{", end='')
    print_poly(poly_p, "")
    print("}", end='')
    print_poly(poly_q, "")
        
# def pp(poly: list[float]) -> str:
#     '''
#     Converts polynomial into string.
#     '''
#     res = ""
#     for i in range(len(poly)-1, -1, -1):
#         coeff = poly[i]
#         if coeff == 0:
#             continue

#         if coeff < 0 or i == len(poly)-1:
#             mark = ''
#         elif coeff > 0:
#             mark = '+'
        
#         if coeff.is_integer():
#             coeff = int(coeff)
#         else:
#             coeff = round(coeff, 2)
            
#         if i == 0:
#             res += f"{mark}{coeff}"
#         elif i == 1:
#             res += f"{mark}{coeff}x"
#         else:
#             res += f"{mark}{coeff}(x^{i})"
#     return res

def poly_to_str(poly: list[float]) -> str:
    """
    Converts polynomial into string.
    """
    res = ""
    space = ""
    for i in range(len(poly) - 1, -1, -1):
        coeff = poly[i]
        if coeff == 0:
            continue

        # sign
        if not res:
            mark = "-" if coeff < 0 else ""
        else:
            mark = "-" if coeff < 0 else "+"

        c = abs(coeff)

        # numeric formatting
        if float(c).is_integer():
            c_str = str(int(c))
        else:
            c_str = str(round(c, 2))

        # term formatting (omit 1 for non-constant terms)
        if i == 0:
            term = c_str
        elif i == 1:
            term = ("x" if c == 1 else f"{c_str}x")
        else:
            term = (f"x^{i}" if c == 1 else f"{c_str}x^{i}")

        res += f"{space}{mark}{term}"
        space = " "

    return res

def trim(coeff: list[float]) -> list[float]:
        while coeff and coeff[-1] == 0:
            coeff.pop()
        return coeff
    

def decomp_poly(poly: list[float], dcData: tuple) -> str:
    '''
    Converts decomposition data into string.
    '''
    multA, i, decomp_p, decomp_q = dcData
    if poly == []:
        return ""
    # i=0
    if i == 0:
        return poly_to_str(poly)
    else:
        coeff_p = trim(poly[i:])
        coeff_q = trim(poly[:i])
        if multA and coeff_p:
            leading_coeff = coeff_p[-1]
            if leading_coeff != 0:
                coeff_p = [c / leading_coeff for c in coeff_p]
        coeff = f"{poly[-1]}" if multA else ""
        if i == 1:
            res = f"{coeff}x"
        else:
            res = f"{coeff}x^{i}"
        
        q = f" + ({decomp_poly(coeff_q, decomp_q)})" if coeff_q != [] else ""
        return f"({res})[ {decomp_poly(coeff_p, decomp_p)} ]{q}" 
