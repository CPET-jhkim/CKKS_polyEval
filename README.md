# polyEval_v3

`polyEval_v2`의 알고리즘 동작을 유지하면서 알고리즘 코드는 `src/`, 실행 및 회귀 테스트는 `test/`로 분리한 버전이다. JSON cache는 별도 metadata wrapper 없이 계수 타입을 최상위 key로 직접 사용한다.

아래 명령은 저장소 루트에서 실행한다.

```powershell
# 기존 단위 테스트
python -m unittest discover -s test -p "test_*.py"

# PS-Method와 2~8차 전체 타입 비교
python -m test.test_compare_ps_method

# 2~8차 전체 타입의 최적 분해식을 JSON으로 저장
python -m test.test_save_all_json

# 실제 계수 또는 계수 타입으로 최적 분해식 계산
python -m test.test_optimal_decomposition --coeff "1,2.5,0,4"
python -m test.test_optimal_decomposition --type "IF0I"

# 특정 차수 전체 타입의 계산복잡도 분포 출력
python -m test.test_degree_complexity 4
python -m test.test_degree_complexity 4 --output data/degree4_summary.json
```

각 차수 `n`의 계수 타입 수는 최고차항 `I/F` 두 경우와 나머지 `n`개 항의 `0/I/F` 경우를 합쳐 `2 * 3^n`개이다. 2~8차 전체는 19,674개이다. 장시간 작업인 JSON 저장 테스트는 각 차수가 끝날 때 파일을 갱신하며, 기본적으로 기존 결과를 이어서 계산한다.
