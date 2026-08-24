import unittest
import sys
from pathlib import Path

if __package__ and "." in __package__:
    from ..src import NULL_DECOMP, Poly, cal_PSMethod, make_all_poly_types
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src import NULL_DECOMP, Poly, cal_PSMethod, make_all_poly_types


class PSMethodTest(unittest.TestCase):
    def test_every_degree_2_type_is_evaluated(self):
        for coeff_type in make_all_poly_types(2):
            with self.subTest(coeff_type="".join(coeff_type)):
                self.assertIsNot(cal_PSMethod(Poly(coeff_type)), NULL_DECOMP)


if __name__ == "__main__":
    unittest.main()
