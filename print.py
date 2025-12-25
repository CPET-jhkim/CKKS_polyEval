# print.py
from complexity import Complexity

def print_poly(poly: list[float], title: str="다항식:\t") -> None:
    print(f"{title}", end='')
    print(pp(poly))
                
def print_poly_type(poly_type: list[str], title: str="타입:\t") -> None:
    if title:
        print(title, end='')
    print(f"({', '.join(poly_type)})")

def print_poly_sep(i: int, poly_p: list[float], poly_q: list[float]) -> None:
    # 다항식 분해
    print("분해식: ", end='')
    print(f"(x^{i})", end="")
    print("{", end='')
    print_poly(poly_p, "")
    print("}", end='')
    print_poly(poly_q, "")
    

def print_step(poly: list[float], i: int, poly_p: list[float], poly_q: list[float],
               comp_i: Complexity, comp_p: Complexity, comp_q: Complexity, comp_piq: Complexity,
               made_powers: set[int], mp: set[int]):
    print('#'*20)
    print_poly(poly)
    print()
    print_poly_sep(i, poly_p, poly_q)
    print()
    print(f"복잡도: ")
    print(f"Depth:\t{comp_i.depth}\t{comp_p.depth}\t{comp_q.depth}\t=>\t{comp_piq.depth}")
    print(f"CMult:\t{comp_i.cmult}\t{comp_p.cmult}\t{comp_q.cmult}\t=>\t{comp_piq.cmult}")
    print(f"PMult:\t{comp_i.pmult}\t{comp_p.pmult}\t{comp_q.pmult}\t=>\t{comp_piq.pmult}")
    print(f"Add:\t{comp_i.add}\t{comp_p.add}\t{comp_q.add}\t=>\t{comp_piq.add}")
    print(f"생성 차수: {made_powers} => {mp}")
    
# 다항식 출력 문자열
def pp(poly: list[float]) -> str:
    res = ""
    for i in range(len(poly)-1, -1, -1):
        coeff = poly[i]
        if coeff == 0:
            continue

        if coeff < 0 or (i == len(poly)-1 and i != 0):
            mark = ''
        elif coeff > 0:
            mark = '+'
        
        if coeff.is_integer():
            coeff = int(coeff)
            
        if i == 0:
            res += f"{mark}{coeff}"
        elif i == 1:
            res += f"{mark}{coeff}x"
        else:
            res += f"{mark}{coeff}(x^{i})"
    return res

def trim(coeff: list[float]) -> list[float]:
        while coeff and coeff[-1] == 0:
            coeff.pop()
        return coeff
    
# 분해식에 따라 다항식 출력 문자열 반환
def decomp_poly(poly: list[float], dcData: tuple) -> str:
    multA, i, decomp_p, decomp_q = dcData
    ###########        
    # poly가 비어있는 경우
    if poly == []:
        return ""
    # i=0
    if i == 0:
        return pp(poly)
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
            
        return f"{res}[{decomp_poly(coeff_p, decomp_p)}]{decomp_poly(coeff_q, decomp_q)}"
