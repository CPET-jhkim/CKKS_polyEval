# all_n_degree_polys.py
from contextlib import ExitStack, redirect_stdout
import sys
from poly_util import make_deg_poly, make_all_polys, make_type_poly
from algorithm import cal_polyEval, cal_PSMethod
from basic_class import Poly, Decomp


def print_result(result: Decomp):
    print(f"dcmp:\t{result.restore_dcmp()}")
    result.comp.print_params()
    print("\t---")


def run_eval(poly: Poly):
    # PS-Method
    print(f"\tPS-Method.")
    print_result(cal_PSMethod(poly))

    # Our Method
    print(f"\tOur Method.")
    print_result(cal_polyEval(poly, consider_axi=True))


# if __name__ == "__main__":
#     for deg in range(7, 8):
#         option = "even" if deg % 2 == 0 else "odd"

#         print(f"Degree={deg}")

#         for title, kind in [
#             # ("normal Polynomial.", "all"),
#             (f"{option} Polynomial.", option),
#         ]:
#             print(title)

#             coeff = make_deg_poly(deg, kind)
#             poly = Poly(coeff)
#             poly.print()

#             run_eval(poly)

#         print("\n##############\n")

# n차 모든 경우의 수에 대하여 PS-Method와 Our Method와 비교, 오류 확인.
if __name__ == '__main__':
    deg = 9
    coeffs = make_all_polys(deg)

    # 원래 터미널 stdout 저장
    terminal = sys.stdout
    print(f"Degree {deg}, total {len(coeffs)} polynomials.")

    with ExitStack() as stack:
        out_file = stack.enter_context(
            open(f"OPD_deg{deg}.txt", "w", encoding="utf-8")
        )
        stack.enter_context(redirect_stdout(out_file))

        print(f"Degree={deg}")

        for coeff in coeffs:
            poly = Poly(coeff)
            ctype = "".join(poly.coeff_type)
            dcmp_PS = cal_PSMethod(poly)
            dcmp_Ours = cal_polyEval(poly, consider_axi=True)

            flag = "PASS" if dcmp_Ours.comp <= dcmp_PS.comp else "FAIL"

            PS_comp = dcmp_PS.comp.return_params()
            OURS_comp = dcmp_Ours.comp.return_params()

            # 파일 출력
            if PS_comp == OURS_comp:
                print(f"{ctype}\t{flag}\t{PS_comp}")
            else:
                if flag == "FAIL":
                    print(f"{ctype}\t{flag}")
                    print_result(dcmp_PS)
                    print_result(dcmp_Ours)
                else:
                    print(f"{ctype}\t{flag}\t{PS_comp}\t{OURS_comp}")

            # 터미널 출력 (예: FAIL만)
            if flag == "FAIL":
                print(f"{ctype}\t{flag}\t{PS_comp}\t{OURS_comp}", file=terminal)
                    
                    
# 특정 타입 다항식 생성, 비교
# if __name__ == "__main__":
#     ctype = "00FIF0FI"
#     coeff = make_type_poly(list(ctype))
#     poly = Poly(coeff)
#     poly.print()
#     run_eval(poly)
    