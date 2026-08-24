import unittest
import sys
from pathlib import Path

if __package__ and "." in __package__:
    from ..src import make_all_poly_types
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src import make_all_poly_types


class PolynomialTypeGenerationTest(unittest.TestCase):
    def test_degree_2_to_8_contains_every_unique_type(self):
        for degree in range(2, 9):
            poly_types = make_all_poly_types(degree)
            self.assertEqual(len(poly_types), 2 * 3 ** degree)
            self.assertEqual(len({tuple(poly_type) for poly_type in poly_types}), len(poly_types))
            self.assertTrue(all(len(poly_type) == degree + 1 for poly_type in poly_types))
            self.assertTrue(all(poly_type[-1] in {"I", "F"} for poly_type in poly_types))


if __name__ == "__main__":
    unittest.main()
