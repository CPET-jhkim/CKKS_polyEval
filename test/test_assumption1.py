import unittest
import os
import sys
import tempfile
from pathlib import Path

if __package__ and "." in __package__:
    from ..src.algorithm import cal_polyEval
    from ..src.basic_class import Poly, NULL_DECOMP
    from ..src.save import (
        deserialize_search_cache,
        load_opd_caches,
        save_opd_caches,
        select_decomp_for_coeff,
        serialize_search_cache,
    )
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src.algorithm import cal_polyEval
    from src.basic_class import Poly, NULL_DECOMP
    from src.save import (
        deserialize_search_cache,
        load_opd_caches,
        save_opd_caches,
        select_decomp_for_coeff,
        serialize_search_cache,
    )


class Assumption1Test(unittest.TestCase):
    def test_assumption1_forces_nonzero_lower_coefficients_to_float(self):
        '''
        기능:
            Assumption 1 적용 시 정규화된 비영 하위 계수가 모두 F인지 검증한다.
        입력:
            타입 다항식 [I, I, 0, F]와 ax^1 분리.
        출력:
            quotient 타입 경우가 [F, 0, I] 하나이면 테스트 성공.
        '''
        poly = Poly(["I", "I", "0", "F"])
        cases = poly.seperate_cases(1, multA=True, assumption1=True)
        coeff_types = [poly_p.coeff_type for poly_p, _ in cases]

        self.assertEqual(coeff_types, [["F", "0", "I"]])

    def test_without_assumption1_enumerates_integer_and_float_cases(self):
        '''
        기능:
            Assumption 1 제거 시 각 비영 하위 계수의 I/F 경우를 모두 생성하는지 검증한다.
        입력:
            타입 다항식 [I, I, F, F]와 ax^1 분리.
        출력:
            두 비영 하위 항에 대한 네 가지 quotient 타입이 생성되면 테스트 성공.
        '''
        poly = Poly(["I", "I", "F", "F"])
        cases = poly.seperate_cases(1, multA=True, assumption1=False)
        coeff_types = {
            tuple(poly_p.coeff_type)
            for poly_p, _ in cases
        }

        self.assertEqual(
            coeff_types,
            {
                ("I", "I", "I"),
                ("I", "F", "I"),
                ("F", "I", "I"),
                ("F", "F", "I"),
            },
        )

    def test_zero_type_is_not_branched(self):
        '''
        기능:
            ax^i 정규화에서 0 계수가 I/F 분기 대상에 포함되지 않는지 검증한다.
        입력:
            타입 다항식 [I, I, 0, F]와 assumption1=False.
        출력:
            0 위치가 유지된 두 가지 quotient 타입만 생성되면 테스트 성공.
        '''
        poly = Poly(["I", "I", "0", "F"])
        cases = poly.seperate_cases(1, multA=True, assumption1=False)
        coeff_types = {
            tuple(poly_p.coeff_type)
            for poly_p, _ in cases
        }

        self.assertEqual(
            coeff_types,
            {
                ("I", "0", "I"),
                ("F", "0", "I"),
            },
        )

    def test_type_polynomial_runs_opd_with_both_settings(self):
        '''
        기능:
            타입 다항식이 두 Assumption 1 설정에서 OPD 전체 탐색을 완료하는지 검증한다.
        입력:
            타입 다항식 [I, F, I, F].
        출력:
            두 설정 모두 NULL_DECOMP가 아닌 결과를 반환하면 테스트 성공.
        '''
        poly_type = ["I", "F", "I", "F"]
        result_assumed = cal_polyEval(poly_type, assumption1=True)
        result_general = cal_polyEval(poly_type, assumption1=False)

        self.assertIsNot(result_assumed, NULL_DECOMP)
        self.assertIsNot(result_general, NULL_DECOMP)

    def test_numeric_polynomial_remains_supported(self):
        '''
        기능:
            기존 수치 계수 입력 인터페이스가 재설계 후에도 동작하는지 검증한다.
        입력:
            수치 다항식 [1.2, 2.3, 3.4, 4.5].
        출력:
            OPD가 NULL_DECOMP가 아닌 결과를 반환하면 테스트 성공.
        '''
        result = cal_polyEval(
            [1.2, 2.3, 3.4, 4.5],
            assumption1=True,
        )

        self.assertIsNot(result, NULL_DECOMP)

    def test_initial_and_intermediate_caches_store_search_tuples(self):
        '''
        기능:
            초기/중간 power cache가 분리되고 지정된 SearchResult tuple을 저장하는지 검증한다.
        입력:
            타입 다항식 [I, F, I, F, I, F]와 빈 초기/중간 cache dictionary.
        출력:
            초기 cache에는 {0,1}, 중간 cache에는 그 외 power 입력의
            (power cache, decomposition, complexity) tuple이 저장되면 테스트 성공.
        '''
        initial_cache = {}
        intermediate_cache = {}

        cal_polyEval(
            ["I", "F", "I", "F", "I", "F"],
            assumption1=False,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )

        initial_results = [
            item
            for results in initial_cache.values()
            for item in results
        ]
        intermediate_results = [
            item
            for results in intermediate_cache.values()
            for item in results
        ]

        self.assertGreater(len(initial_results), 0)
        self.assertGreater(len(intermediate_results), 0)
        self.assertTrue(
            all(len(item) == 3 for item in initial_results + intermediate_results)
        )
        self.assertTrue(
            all(
                powers == frozenset({(0, 0), (1, 0)})
                for powers, _, _ in initial_results
            )
        )
        self.assertTrue(
            all(
                powers != frozenset({(0, 0), (1, 0)})
                for powers, _, _ in intermediate_results
            )
        )

    def test_cached_search_reuses_existing_candidates(self):
        '''
        기능:
            동일한 타입과 power cache의 반복 탐색이 기존 캐시 후보를 재사용하는지 검증한다.
        입력:
            동일 타입 다항식과 동일한 초기/중간 cache dictionary를 사용한 두 번의 OPD 호출.
        출력:
            두 번째 호출 후 cache tuple 수가 증가하지 않으면 테스트 성공.
        '''
        initial_cache = {}
        intermediate_cache = {}
        poly_type = ["I", "F", "I", "F"]

        cal_polyEval(
            poly_type,
            assumption1=False,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )
        count_before = sum(len(results) for results in initial_cache.values())
        count_before += sum(len(results) for results in intermediate_cache.values())

        cal_polyEval(
            poly_type,
            assumption1=False,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )
        count_after = sum(len(results) for results in initial_cache.values())
        count_after += sum(len(results) for results in intermediate_cache.values())

        self.assertEqual(count_before, count_after)

    def test_search_cache_serialization_round_trip(self):
        '''
        기능:
            SearchResult cache의 메모리 직렬화와 복원이 후보를 보존하는지 검증한다.
        입력:
            assumption1=False로 생성한 초기 SearchResult cache.
        출력:
            복원 전후 상태 수, 후보 수, complexity가 같으면 테스트 성공.
        '''
        initial_cache = {}
        intermediate_cache = {}

        cal_polyEval(
            ["I", "I", "F", "F"],
            assumption1=False,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )
        serialized = serialize_search_cache(initial_cache)
        restored = deserialize_search_cache(serialized)

        self.assertEqual(len(initial_cache), len(restored))
        self.assertEqual(
            sum(len(results) for results in initial_cache.values()),
            sum(len(results) for results in restored.values()),
        )

        original_complexities = sorted(
            comp.return_params()
            for results in initial_cache.values()
            for _, _, comp in results
        )
        restored_complexities = sorted(
            comp.return_params()
            for results in restored.values()
            for _, _, comp in results
        )
        self.assertEqual(original_complexities, restored_complexities)

    def test_opd_cache_file_round_trip(self):
        '''
        기능:
            초기/중간 SearchResult cache가 JSON 파일로 왕복되는지 검증한다.
        입력:
            타입 다항식으로 생성한 초기/중간 cache와 임시 JSON 경로.
        출력:
            저장·복원 후 두 cache의 상태 수와 후보 수가 같으면 테스트 성공.
        '''
        initial_cache = {}
        intermediate_cache = {}

        cal_polyEval(
            ["I", "F", "I", "F"],
            assumption1=False,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, "opd_cache.json")
            save_opd_caches(initial_cache, intermediate_cache, filename)
            restored_initial, restored_intermediate = load_opd_caches(filename)

        self.assertEqual(len(initial_cache), len(restored_initial))
        self.assertEqual(len(intermediate_cache), len(restored_intermediate))
        self.assertEqual(
            sum(len(results) for results in initial_cache.values()),
            sum(len(results) for results in restored_initial.values()),
        )
        self.assertEqual(
            sum(len(results) for results in intermediate_cache.values()),
            sum(len(results) for results in restored_intermediate.values()),
        )

    def test_select_decomp_for_actual_coefficients(self):
        '''
        기능:
            실제 ax^i 정규화 타입과 일치하는 조건부 cache 후보가 선택되는지 검증한다.
        입력:
            타입 I/I/I/F로 생성한 cache와 실제 계수 [1,5,4,2.5].
        출력:
            선택 결과가 존재하고 실제 계수 직접 탐색 complexity와 같으면 테스트 성공.
        '''
        initial_cache = {}
        intermediate_cache = {}
        coeff = [1, 5, 4, 2.5]

        cal_polyEval(
            ["I", "I", "I", "F"],
            assumption1=False,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, "opd_cache.json")
            save_opd_caches(initial_cache, intermediate_cache, filename)
            restored_initial, _ = load_opd_caches(filename)
            selected = select_decomp_for_coeff(
                coeff,
                restored_initial,
                assumption1=False,
            )

        direct = cal_polyEval(coeff, assumption1=False)

        self.assertIsNot(selected, NULL_DECOMP)
        self.assertEqual(selected.comp, direct.comp)


if __name__ == "__main__":
    unittest.main()
