import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

if __package__ and "." in __package__:
    from ..src.save import get_poly_type_key, load_cache, save_cache
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src.save import get_poly_type_key, load_cache, save_cache


class SaveFormatTest(unittest.TestCase):
    def test_poly_type_key_contains_only_coefficient_types(self):
        '''
        기능:
            Assumption 1 설정과 관계없이 계수타입만 cache key에 포함되는지 검증한다.
        입력:
            계수타입 [I, 0, F]와 두 Assumption 1 설정.
        출력:
            두 key가 모두 I0F이면 테스트 성공.
        '''
        self.assertEqual(get_poly_type_key(["I", "0", "F"], True), "I0F")
        self.assertEqual(get_poly_type_key(["I", "0", "F"], False), "I0F")

    def test_save_cache_uses_selective_line_breaks(self):
        '''
        기능:
            최상위 decomposition과 p/q만 여러 줄로 쓰고 다른 list/dict는
            한 줄로 저장하며 기존 A1| key를 제거하는지 검증한다.
        입력:
            list, 일반 dictionary, p child를 포함한 소형 cache.
        출력:
            JSON 왕복과 key 및 줄바꿈 형식이 모두 일치하면 테스트 성공.
        '''
        cache = {
            "A1|IF": {
                "coeff_type": ["I", "F"],
                "complexity": {"depth": 1, "cmult": 0},
                "term_plans": [
                    {
                        "degree": 1,
                        "coeff_type": "I",
                        "multA": False,
                        "route_ops": [],
                    },
                    {
                        "degree": 2,
                        "coeff_type": "F",
                        "multA": True,
                        "route_ops": [
                            {"kind": "coeff", "lhs": 1, "rhs": 1}
                        ],
                    },
                ],
                "p": {
                    "coeff_type": ["I"],
                    "complexity": {"depth": 0},
                    "p": None,
                    "q": None,
                },
                "q": None,
            }
        }

        with tempfile.TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, "cache.json")
            save_cache(cache, filename)

            with open(filename, "r", encoding="utf-8") as file:
                raw = file.read()

            with open(filename, "r", encoding="utf-8") as file:
                parsed = json.load(file)

            loaded = load_cache(filename)

        self.assertEqual(list(parsed), ["IF"])
        self.assertEqual(list(loaded), ["IF"])
        self.assertIn('"IF": {\n', raw)
        self.assertIn('"p": {\n', raw)
        self.assertIn('"coeff_type": ["I", "F"]', raw)
        self.assertIn('"complexity": {"cmult": 0, "depth": 1}', raw)
        self.assertIn(
            '"term_plans": [\n'
            '      {"coeff_type": "I", "degree": 1, "multA": false, "route_ops": []},\n'
            '      {"coeff_type": "F", "degree": 2, "multA": true, "route_ops": '
            '[{"kind": "coeff", "lhs": 1, "rhs": 1}]}\n'
            '    ]',
            raw,
        )
        self.assertNotIn('"coeff_type": [\n', raw)
        self.assertNotIn('"complexity": {\n', raw)
        self.assertNotIn('"route_ops": [\n', raw)
        self.assertNotIn('A1|', raw)


if __name__ == "__main__":
    unittest.main()
