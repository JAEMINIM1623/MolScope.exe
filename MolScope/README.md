# MolScope — ChemDraw 구조 기반 분자 모델링/분석 도구

ChemDraw에서 그린 분자 구조를 불러와 **분석 → 3D 생성 → Gaussian 입력(.gjf) →
계산결과 파싱 → 사슬길이 외삽**까지 하나의 GUI에서 처리하는 독립 프로그램입니다.
특정 화학계에 종속되지 않으며, 임의의 유기분자/고분자 반복단위에 적용됩니다.

```
ChemDraw (.cdx/.cdxml) ──▶ MolScope ──▶ Gaussian(.gjf) ──▶ 계산 ──▶ MolScope(.log 파싱/외삽)
                                                        └─▶ Multiwfn (파동함수 해석)
```

---

## 1. ChemDraw에서 구조 가져오기 (3가지 경로)

| 경로 | ChemDraw 조작 | 안정성 |
|---|---|---|
| **① .cdxml 저장 (권장)** | File → Save As → **ChemDraw XML (*.cdxml)** | ★★★ |
| ② .cdx 네이티브 | 그냥 저장한 .cdx 그대로 열기 | ★★ (실패 시 ①로) |
| ③ 클립보드 SMILES | 구조 선택 → Edit → **Copy As → SMILES** → 탭1에 붙여넣기 | ★★★ |

- 한 페이지에 **여러 구조**를 그려두면 전부 개별 분자로 들어옵니다.
- `.mol`, `.sdf`, `.smi`도 지원합니다.
- **올리고머용 반복단위**: 양 끝에 attachment point 2개를 찍거나,
  SMILES에 `[*]` 두 개를 넣으세요. 예) `[*]c1ccc([*])s1`

## 2. 탭별 기능

| 탭 | 기능 |
|---|---|
| 1. 구조 불러오기 | .cdx/.cdxml/.mol/.sdf/.smi 열기, SMILES 붙여넣기, 목록 관리 |
| 2. 분석 | 기술자(MW, TPSA, logP, Fsp³ 등) + **공액 분석**(공액계 개수/크기, 최장 공액 경로, 고리-고리 단일결합) + 이면각 스캔 후보 + 2D 미리보기 + 전체 CSV 내보내기 |
| 3. Gaussian 입력 생성 | ETKDGv3+MMFF 3D 생성 → route 프리셋(opt/freq, TDDFT, SMD, 스캔, NMR, 양이온) → 이면각 modredundant 자동/수동 → **전기장 시리즈**(MV/m → `Field=축+N` 자동 변환) |
| 4. 올리고머 시리즈 | 반복단위 → n-mer 일괄 생성 + opt/TDDFT .gjf. 파일명에 `_n숫자` 태그 자동 부여 |
| 5. 계산결과 파싱 | .log 일괄 파싱: 정상종료/허수진동, SCF/G, HOMO/LUMO/gap, 쌍극자, S1. 스캔 PES → CSV + Boltzmann 확률분포 그림 |
| 6. 사슬길이 외삽 | `_n` 태그 잡 자동 수집 → **1/n 선형** + **Kuhn 식** → 극한값·유효공액길이(ECL)·그래프 |

## 3. Windows 실행파일 만들기

.exe는 Windows에서 빌드해야 합니다(크로스컴파일 불가). 두 가지 방법:

### 방법 A — GitHub Actions 자동 빌드 (권장, 좌석배치 프로그램 때와 동일)
1. 이 폴더를 GitHub 저장소로 push
2. Actions 탭 → **Build Windows EXE** 완료 대기 (~10분)
3. Artifacts에서 `MolScope-windows` 다운로드 → `MolScope.exe`

### 방법 B — 내 PC에서 빌드
Python 3.10+ 설치된 Windows에서 `build_exe.bat` 더블클릭 → `dist\MolScope.exe`

> 빌드된 exe는 **단일 파일**이며, 다른 PC에 복사만 하면 Python 없이 실행됩니다.
> RDKit 포함으로 용량이 200–300 MB 정도 나오는 것이 정상입니다.

### 소스로 바로 실행
```
pip install -r requirements.txt
python MolScope.py
```

## 4. 권장 워크플로 예시 (임의의 공액 고분자)

1. ChemDraw에서 반복단위를 그리고 양 끝에 attachment point → `.cdxml` 저장
2. 탭1에서 열기 → 탭2에서 공액 분석으로 구조 확인
3. 탭4에서 n = 1 2 3 4 5 6 8 시리즈 생성 (TDDFT 체크)
4. 서버에서 Gaussian 실행
5. 탭5에서 로그 폴더 파싱 → 실패 잡은 오류 열에서 원인 확인
6. 탭6에서 gap_eV / homo_eV / s1_eV 외삽 → 고분자 극한값 + ECL 보고

이면각 스캔이 필요하면 탭3에서 "이면각 relaxed scan" 프리셋 선택 —
자동 탐색된 고리-고리 이면각으로 `D a b c d S 36 10.0` 블록이 생성되고,
탭5의 스캔 내보내기가 PES + Boltzmann 확률분포(P(θ)) 그림까지 만들어 줍니다.

## 5. 역할 경계

MolScope는 계산의 **준비와 해석**을 담당합니다.
- DFT/TDDFT 계산 자체 → **Gaussian**
- 파동함수 해석(LOL-π, hole-electron Sr/D/Δσ/t, IFCT, ESP 정량) → **Multiwfn** (http://sobereva.com/multiwfn)

TDDFT 프리셋에 들어 있는 `IOp(9/40=4)`는 Multiwfn의 여기상태 해석에 필요한
CI 계수 출력 옵션이므로 지우지 마세요.

## 6. 폴더 구조

```
MolScope/
├─ MolScope.py              # 실행 진입점
├─ molscope/
│  ├─ core.py               # 화학 코어 (GUI 비의존, 스크립트 재사용 가능)
│  └─ gui.py                # Tkinter GUI
├─ requirements.txt
├─ build_exe.bat            # 로컬 Windows 빌드
└─ .github/workflows/build.yml   # GitHub Actions 자동 빌드
```
