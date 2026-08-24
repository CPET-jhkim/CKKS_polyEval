import os
import sys
import tempfile
import unittest
from pathlib import Path

if __package__ and "." in __package__:
    from ..src.algorithm import cal_polyEval, cal_polyEval_candidates
    from ..src.poly_util import solve_xn_routes, solve_xn_routes_depth_limited
    from ..src.save import load_opd_caches, save_opd_caches
else:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src.algorithm import cal_polyEval, cal_polyEval_candidates
    from src.poly_util import solve_xn_routes, solve_xn_routes_depth_limited
    from src.save import load_opd_caches, save_opd_caches


class DepthLimitedSearchTest(unittest.TestCase):
    def test_power_routes_respect_final_depth_only(self):
        '''
        기능:
            x^i의 자체 최소 depth가 아니라 전달된 최종 depth budget만 적용되는지
            검증한다.
        입력:
            초기 power {0,1}, 목표 x^4, 최대 depth 3.
        출력:
            depth 2와 depth 3의 x^4 경로가 모두 존재하면 성공.
        '''
        routes = solve_xn_routes_depth_limited(
            False,
            4,
            {0, 1},
            3,
        )

        self.assertEqual({xi.depth for xi in routes}, {2, 3})

    def test_x10_includes_nonminimal_x4_intermediate(self):
        '''
        기능:
            최종 x^10의 log-depth 안에서 중간 x^4를 depth 3으로 만드는 경로도
            실제 탐색공간에 포함되는지 검증한다.
        입력:
            초기 power {0,1}, 목표 x^10, 최대 depth ceil(log2(10))=4.
        출력:
            x^2,x^3,x^4,x^6,x^10 경로를 사용하는 depth 4 후보가 있으면 성공.
        '''
        routes = solve_xn_routes_depth_limited(
            False,
            10,
            {0, 1},
            4,
        )
        expected_ops = {
            ("pure", 1, 1),
            ("pure", 1, 2),
            ("pure", 1, 3),
            ("pure", 3, 3),
            ("pure", 4, 6),
        }

        self.assertGreater(len(routes), 0)
        self.assertTrue(all(xi.depth <= 4 for xi in routes))
        self.assertTrue(
            any(
                xi.depth == 4
                and set(xi.route_ops) == expected_ops
                and xi.power_depths[4] == 3
                for xi in routes
            )
        )

    def test_search_cache_distinguishes_power_depth(self):
        '''
        기능:
            동일한 power 집합이라도 실제 생성 depth가 다르면 서로 다른 OPD
            검색상태로 저장되는지 검증한다.
        입력:
            x^4 depth가 각각 2와 3인 cache 및 동일한 타입 다항식.
        출력:
            intermediate cache에 두 power-depth key가 모두 존재하면 성공.
        '''
        initial_cache = {}
        intermediate_cache = {}
        poly_type = ["I", "0", "0", "0", "I"]

        cal_polyEval_candidates(
            poly_type,
            made_powers={0: 0, 1: 0, 4: 2},
            depth_budget=3,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )
        cal_polyEval_candidates(
            poly_type,
            made_powers={0: 0, 1: 0, 4: 3},
            depth_budget=3,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )

        stored_states = {
            powers
            for results in intermediate_cache.values()
            for powers, _, _ in results
        }
        self.assertIn(
            frozenset({(0, 0), (1, 0), (4, 2)}),
            stored_states,
        )
        self.assertIn(
            frozenset({(0, 0), (1, 0), (4, 3)}),
            stored_states,
        )

    def test_all_decompositions_respect_only_depth_budget(self):
        '''
        기능:
            OPD 후보가 CMult 상한 없이 생성되고 depth 상한만 만족하는지 검증한다.
        입력:
            타입 다항식 FFFFF와 depth budget 3.
        출력:
            복수 후보, depth 제한, degree-1보다 큰 CMult 후보가 확인되면 성공.
        '''
        candidates = cal_polyEval_candidates(
            ["F", "F", "F", "F", "F"],
            depth_budget=3,
            assumption1=True,
            consider_axi=True,
        )

        self.assertGreater(len(candidates), 1)
        self.assertTrue(all(dcmp.comp.depth <= 3 for dcmp in candidates))
        self.assertTrue(
            any(dcmp.comp.cmult > 3 for dcmp in candidates),
            [dcmp.comp.return_params() for dcmp in candidates],
        )

    def test_persisted_intermediate_cache_is_reused_by_another_root(self):
        '''
        기능:
            첫 초기 다항식의 탐색 cache를 JSON으로 복원한 뒤 다른 초기
            다항식이 동일한 중간 상태를 중복 저장하지 않고 재사용하는지 검증한다.
        입력:
            첫 입력 FFF, 두 번째 입력 IFF, 공통 Assumption 1 search cache.
        출력:
            공통 FF 상태의 후보 수가 두 번째 탐색 전후 동일하면 성공.
        '''
        initial_cache = {}
        intermediate_cache = {}
        cal_polyEval(
            ["F", "F", "F"],
            assumption1=True,
            consider_axi=True,
            initial_cache=initial_cache,
            intermediate_cache=intermediate_cache,
        )

        shared_key = (
            ("F", "F"),
            True,
            True,
            True,
            2,
        )
        self.assertIn(shared_key, initial_cache)

        with tempfile.TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, "search_cache.json")
            save_opd_caches(
                initial_cache,
                intermediate_cache,
                filename,
            )
            restored_initial, restored_intermediate = load_opd_caches(filename)

        count_before = len(restored_initial[shared_key])
        cal_polyEval(
            ["I", "F", "F"],
            assumption1=True,
            consider_axi=True,
            initial_cache=restored_initial,
            intermediate_cache=restored_intermediate,
        )
        count_after = len(restored_initial[shared_key])

        self.assertEqual(count_before, count_after)


if __name__ == "__main__":
    unittest.main()
