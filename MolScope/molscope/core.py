# -*- coding: utf-8 -*-
"""
molscope.core — 화학 계산 코어 (GUI 비의존)

이 모듈은 Tkinter 를 import 하지 않습니다. 따라서 헤드리스 서버,
Jupyter, 배치 스크립트에서도 그대로 재사용 가능합니다.

주요 기능
---------
  * ChemDraw 파일 입력 : .cdx (네이티브), .cdxml, .mol, .sdf, .smi
  * 클립보드 SMILES    : ChemDraw 의 Edit > Copy As > SMILES
  * 분자 기술자        : 조성/크기/극성/공액 관련 지표
  * 공액 분석          : 최대 공액계 크기, 최장 공액 경로, 회전가능 이면각
  * 3D 구조 생성       : ETKDGv3 다중 컨포머 + MMFF/UFF 최적화
  * Gaussian 입력 생성 : 프리셋 기반 route + modredundant + Field 시리즈
  * 로그 파싱          : SCF, HOMO/LUMO, 쌍극자, 진동수, TDDFT, 스캔 PES
  * 사슬길이 외삽      : 1/n 선형 + Kuhn 식 -> 극한값, 유효공액길이
"""

from __future__ import annotations

import csv
import glob
import math
import os
import re
from dataclasses import dataclass, field
from typing import Iterable, Sequence

HARTREE_EV = 27.211386245988
HARTREE_KCAL = 627.5094740631
AU_FIELD_V_PER_M = 5.14220675e11

SUPPORTED_INPUT = [
    ("ChemDraw / 구조 파일", "*.cdx *.cdxml *.mol *.sdf *.smi *.txt"),
    ("ChemDraw 네이티브 (*.cdx)", "*.cdx"),
    ("ChemDraw XML (*.cdxml)", "*.cdxml"),
    ("MDL Molfile (*.mol)", "*.mol"),
    ("SD file (*.sdf)", "*.sdf"),
    ("SMILES (*.smi *.txt)", "*.smi *.txt"),
    ("모든 파일", "*.*"),
]


# ---------------------------------------------------------------------------
# 0. RDKit 지연 로딩
# ---------------------------------------------------------------------------

_RDKIT_ERR = (
    "RDKit 을 불러올 수 없습니다.\n"
    "소스로 실행하는 경우:  pip install rdkit\n"
    "실행파일(.exe)인 경우 빌드가 손상되었을 수 있습니다."
)


def rdkit_available() -> tuple[bool, str]:
    try:
        import rdkit  # noqa: F401
        from rdkit import Chem  # noqa: F401
        return True, rdkit.__version__
    except Exception as exc:  # pragma: no cover
        return False, str(exc)


def _chem():
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        return Chem
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(_RDKIT_ERR) from exc


def chemdraw_cdx_supported() -> bool:
    """RDKit 빌드에 ChemDraw 네이티브(.cdx) 파서가 포함되어 있는지."""
    try:
        from rdkit.Chem import rdmolfiles
        return bool(rdmolfiles.HasChemDrawCDXSupport())
    except Exception:
        return False


# ---------------------------------------------------------------------------
# 1. 구조 레코드 & 입출력
# ---------------------------------------------------------------------------

@dataclass
class MolRecord:
    """분자 1개 + 출처 메타데이터."""
    name: str
    mol: object                     # rdkit.Chem.Mol
    source: str = ""                # 파일 경로 또는 'clipboard'
    origin_format: str = ""         # cdx / cdxml / mol / sdf / smiles
    note: str = ""

    @property
    def smiles(self) -> str:
        Chem = _chem()
        try:
            return Chem.MolToSmiles(self.mol)
        except Exception:
            return ""

    @property
    def n_atoms_heavy(self) -> int:
        return self.mol.GetNumHeavyAtoms()

    @property
    def has_dummy(self) -> bool:
        return any(a.GetAtomicNum() == 0 for a in self.mol.GetAtoms())


def _safe_name(base: str, idx: int, total: int) -> str:
    base = re.sub(r"[^0-9A-Za-z가-힣_.-]", "_", base) or "mol"
    return base if total == 1 else f"{base}_{idx + 1:02d}"


def load_chemdraw(path: str, sanitize: bool = True) -> list[MolRecord]:
    """ChemDraw 파일(.cdx / .cdxml)에서 모든 분자를 읽는다.

    ChemDraw 에서 여러 구조를 한 페이지에 그려두면 전부 개별 분자로 들어옵니다.
    """
    Chem = _chem()
    from rdkit.Chem import rdmolfiles

    ext = os.path.splitext(path)[1].lower()
    base = os.path.splitext(os.path.basename(path))[0]

    params = rdmolfiles.CDXMLParserParams()
    params.sanitize = sanitize
    params.removeHs = True
    # Auto: 내용을 보고 CDX/CDXML 자동 판별 (ChemDraw 확장 지원 빌드에서)
    try:
        params.format = rdmolfiles.CDXMLFormat.Auto
    except Exception:
        pass

    try:
        mols = list(rdmolfiles.MolsFromCDXMLFile(path, params))
    except Exception:
        mols = list(rdmolfiles.MolsFromCDXMLFile(path, sanitize, True))

    mols = [m for m in mols if m is not None and m.GetNumAtoms() > 0]
    if not mols:
        raise ValueError(
            f"{os.path.basename(path)} 에서 구조를 찾지 못했습니다.\n"
            "ChemDraw 에서 'Save As > ChemDraw XML (.cdxml)' 로 다시 저장해 보세요."
        )

    fmt = "cdx" if ext == ".cdx" else "cdxml"
    out = []
    for i, m in enumerate(mols):
        nm = m.GetProp("_Name") if m.HasProp("_Name") else ""
        out.append(MolRecord(
            name=nm.strip() or _safe_name(base, i, len(mols)),
            mol=m, source=path, origin_format=fmt,
        ))
    return out


def load_molfile(path: str) -> list[MolRecord]:
    """MDL Molfile(.mol) / SD file(.sdf).

    ChemDraw: File > Save As > MDL Molfile 또는 MDL SDfile
    """
    Chem = _chem()
    ext = os.path.splitext(path)[1].lower()
    base = os.path.splitext(os.path.basename(path))[0]

    if ext == ".sdf":
        supplier = Chem.SDMolSupplier(path, sanitize=True, removeHs=True)
        mols = [m for m in supplier if m is not None]
        fmt = "sdf"
    else:
        m = Chem.MolFromMolFile(path, sanitize=True, removeHs=True)
        mols = [m] if m is not None else []
        fmt = "mol"

    if not mols:
        raise ValueError(f"{os.path.basename(path)} 를 읽지 못했습니다.")

    out = []
    for i, m in enumerate(mols):
        nm = m.GetProp("_Name") if m.HasProp("_Name") else ""
        out.append(MolRecord(name=nm.strip() or _safe_name(base, i, len(mols)),
                             mol=m, source=path, origin_format=fmt))
    return out


def load_smiles_text(text: str, source: str = "clipboard") -> list[MolRecord]:
    """SMILES 텍스트(여러 줄, 'SMILES 이름' 형식 허용).

    ChemDraw: 구조 선택 후 Edit > Copy As > SMILES -> 붙여넣기
    """
    Chem = _chem()
    out: list[MolRecord] = []
    for lineno, raw in enumerate(text.splitlines()):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        smi = parts[0]
        nm = parts[1].strip() if len(parts) > 1 else ""
        m = Chem.MolFromSmiles(smi)
        if m is None:
            out.append(MolRecord(name=nm or f"line{lineno + 1}", mol=None,
                                 source=source, origin_format="smiles",
                                 note=f"파싱 실패: {smi}"))
            continue
        out.append(MolRecord(name=nm or f"mol{len(out) + 1}", mol=m,
                             source=source, origin_format="smiles"))
    good = [r for r in out if r.mol is not None]
    if not good:
        raise ValueError("유효한 SMILES 를 찾지 못했습니다.")
    return out


def load_any(path: str) -> list[MolRecord]:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".cdx", ".cdxml"):
        return load_chemdraw(path)
    if ext in (".mol", ".sdf"):
        return load_molfile(path)
    if ext in (".smi", ".txt"):
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            return load_smiles_text(fh.read(), source=path)
    raise ValueError(f"지원하지 않는 확장자: {ext}")


# ---------------------------------------------------------------------------
# 2. 2D 렌더링 (미리보기)
# ---------------------------------------------------------------------------

def render_2d_png(mol, path: str, size: tuple[int, int] = (420, 320),
                  legend: str = "") -> str:
    Chem = _chem()
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    m = Chem.Mol(mol)
    if not m.GetNumConformers() or m.GetConformer().Is3D():
        AllChem.Compute2DCoords(m)

    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.dummiesAreAttachments = True
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, m, legend=legend)
    drawer.FinishDrawing()
    with open(path, "wb") as fh:
        fh.write(drawer.GetDrawingText())
    return path


# ---------------------------------------------------------------------------
# 3. 분자 기술자 + 공액 분석
# ---------------------------------------------------------------------------

def describe(mol) -> dict:
    """조성 / 크기 / 극성 / 유연성 기술자."""
    _chem()
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen

    d: dict[str, object] = {}
    d["formula"] = rdMolDescriptors.CalcMolFormula(mol)
    d["MW"] = round(Descriptors.MolWt(mol), 3)
    d["exact_mass"] = round(Descriptors.ExactMolWt(mol), 4)
    d["heavy_atoms"] = mol.GetNumHeavyAtoms()
    d["n_rings"] = rdMolDescriptors.CalcNumRings(mol)
    d["n_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
    d["n_rotatable"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    d["TPSA"] = round(rdMolDescriptors.CalcTPSA(mol), 2)
    d["logP_Crippen"] = round(Crippen.MolLogP(mol), 3)
    d["MR_Crippen"] = round(Crippen.MolMR(mol), 3)
    d["HBD"] = rdMolDescriptors.CalcNumHBD(mol)
    d["HBA"] = rdMolDescriptors.CalcNumHBA(mol)
    d["fraction_Csp3"] = round(rdMolDescriptors.CalcFractionCSP3(mol), 4)
    d["heteroatoms"] = rdMolDescriptors.CalcNumHeteroatoms(mol)
    d["formal_charge"] = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    d["radical_electrons"] = Descriptors.NumRadicalElectrons(mol)
    return d


def conjugation_analysis(mol) -> dict:
    """공액계 분석 — π 전자 비국재화 정도의 위상학적 지표.

    DFT 로 LOL-π 나 밴드갭을 계산하기 전에, 구조만으로 후보를 빠르게
    걸러내는 용도. 공액 결합을 그래프로 보고 연결성분을 찾는다.
    """
    Chem = _chem()

    conj_bonds = [b for b in mol.GetBonds() if b.GetIsConjugated()]
    d: dict[str, object] = {
        "n_conjugated_bonds": len(conj_bonds),
        "frac_conjugated_bonds": (round(len(conj_bonds) / mol.GetNumBonds(), 4)
                                  if mol.GetNumBonds() else 0.0),
    }

    # 공액 결합만으로 부분그래프를 만들고 연결성분(공액계) 추출
    adj: dict[int, set[int]] = {}
    for b in conj_bonds:
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)

    seen: set[int] = set()
    systems: list[set[int]] = []
    for start in adj:
        if start in seen:
            continue
        stack, comp = [start], set()
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u)
            stack.extend(v for v in adj[u] if v not in comp)
        seen |= comp
        systems.append(comp)

    systems.sort(key=len, reverse=True)
    d["n_conjugated_systems"] = len(systems)
    d["largest_system_atoms"] = len(systems[0]) if systems else 0
    d["system_sizes"] = [len(s) for s in systems[:8]]

    # 최대 공액계 안에서 최장 경로 (그래프 지름, BFS 2회)
    def longest_path(comp: set[int]) -> int:
        if not comp:
            return 0

        def bfs(src: int) -> tuple[int, int]:
            dist = {src: 0}
            order = [src]
            while order:
                u = order.pop(0)
                for v in adj[u]:
                    if v in comp and v not in dist:
                        dist[v] = dist[u] + 1
                        order.append(v)
            far = max(dist, key=lambda k: dist[k])
            return far, dist[far]

        a, _ = bfs(next(iter(comp)))
        _, diam = bfs(a)
        return diam + 1  # 원자 개수

    d["longest_conjugated_path_atoms"] = longest_path(systems[0]) if systems else 0

    # 백본 평면성 지표: 공액계 내부의 단일결합(회전 가능 = 공액 단절 위험)
    inter_unit = 0
    for b in conj_bonds:
        if (b.GetBondType() == Chem.BondType.SINGLE and not b.IsInRing()
                and b.GetBeginAtom().GetIsAromatic()
                and b.GetEndAtom().GetIsAromatic()):
            inter_unit += 1
    d["inter_ring_single_bonds"] = inter_unit
    return d


def find_rotatable_dihedrals(mol, max_out: int = 40) -> list[dict]:
    """이면각 스캔 후보 목록 (1-based 원자번호, Gaussian 규약).

    고리에 속하지 않는 단일결합 중 양쪽에 이웃이 있는 것을 모두 나열하고,
    두 방향족 고리를 잇는 결합(= 공액 경로를 지배하는 이면각)을 우선 정렬.
    """
    Chem = _chem()
    out: list[dict] = []
    for b in mol.GetBonds():
        if b.IsInRing() or b.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetDegree() < 2 or a2.GetDegree() < 2:
            continue
        n1 = [x for x in a1.GetNeighbors() if x.GetIdx() != a2.GetIdx()]
        n2 = [x for x in a2.GetNeighbors() if x.GetIdx() != a1.GetIdx()]
        if not n1 or not n2:
            continue
        # 무거운 원자를 이웃으로 선호 (H 를 기준축으로 쓰면 해석이 어렵다)
        n1.sort(key=lambda x: (x.GetAtomicNum() == 1, -x.GetDegree()))
        n2.sort(key=lambda x: (x.GetAtomicNum() == 1, -x.GetDegree()))

        ring_ring = a1.GetIsAromatic() and a2.GetIsAromatic()
        conj = b.GetIsConjugated()
        out.append({
            "atoms": (n1[0].GetIdx() + 1, a1.GetIdx() + 1,
                      a2.GetIdx() + 1, n2[0].GetIdx() + 1),
            "label": (f"{n1[0].GetSymbol()}{n1[0].GetIdx() + 1}-"
                      f"{a1.GetSymbol()}{a1.GetIdx() + 1}-"
                      f"{a2.GetSymbol()}{a2.GetIdx() + 1}-"
                      f"{n2[0].GetSymbol()}{n2[0].GetIdx() + 1}"),
            "aryl_aryl": ring_ring,
            "conjugated": conj,
            "priority": (0 if ring_ring else 1 if conj else 2),
        })
    out.sort(key=lambda r: r["priority"])
    return out[:max_out]


# ---------------------------------------------------------------------------
# 4. 3D 구조 생성
# ---------------------------------------------------------------------------

@dataclass
class Embed3DResult:
    mol: object
    energy: float
    n_conformers_tried: int
    forcefield: str
    converged: bool = True


def embed3d(mol, n_conf: int = 20, seed: int = 0xC0FFEE,
            forcefield: str = "MMFF", max_iters: int = 2000) -> Embed3DResult:
    """ETKDGv3 다중 컨포머 생성 -> 힘장 최적화 -> 최저에너지 1개 반환.

    Gaussian 최적화의 시작 구조를 만드는 단계입니다. ChemDraw 2D 구조는
    3D 좌표가 없으므로 이 과정을 반드시 거쳐야 합니다.
    """
    Chem = _chem()
    from rdkit.Chem import AllChem

    m = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.25
    params.useSmallRingTorsions = True

    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=max(1, n_conf), params=params))
    if not cids:
        if AllChem.EmbedMolecule(m, params) != 0:
            params.useRandomCoords = True
            AllChem.EmbedMolecule(m, params)
        cids = [c.GetId() for c in m.GetConformers()]
    if not cids:
        raise RuntimeError("3D 좌표 생성에 실패했습니다. 구조에 오류가 없는지 확인하세요.")

    ff = forcefield.upper()
    if ff == "MMFF" and AllChem.MMFFHasAllMoleculeParams(m):
        res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=max_iters)
        used = "MMFF94s"
    else:
        res = AllChem.UFFOptimizeMoleculeConfs(m, maxIters=max_iters)
        used = "UFF"

    energies = [e for _, e in res]
    flags = [c for c, _ in res]
    best = int(min(range(len(energies)), key=lambda i: energies[i]))

    keep = Chem.Mol(m)
    keep.RemoveAllConformers()
    conf = m.GetConformer(cids[best])
    keep.AddConformer(conf, assignId=True)

    return Embed3DResult(mol=keep, energy=float(energies[best]),
                         n_conformers_tried=len(cids), forcefield=used,
                         converged=(flags[best] == 0))


def to_xyz_block(mol) -> str:
    conf = mol.GetConformer()
    lines = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f" {atom.GetSymbol():<2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    return "\n".join(lines)


def write_xyz(mol, path: str, comment: str = "") -> str:
    block = to_xyz_block(mol)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{mol.GetNumAtoms()}\n{comment}\n{block}\n")
    return path


def write_molfile(mol, path: str) -> str:
    Chem = _chem()
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(Chem.MolToMolBlock(mol))
    return path


# ---------------------------------------------------------------------------
# 5. Gaussian 입력 생성
# ---------------------------------------------------------------------------

ROUTE_PRESETS: dict[str, str] = {
    "구조 최적화 + 진동수":
        "#p B3LYP/6-31G(d) em=GD3BJ opt freq int=ultrafine",
    "구조 최적화 (큰 분자, 저비용)":
        "#p B3LYP/6-31G em=GD3BJ opt int=ultrafine",
    "구조 최적화 + 진동수 (고정밀)":
        "#p B3LYP/6-311G(d,p) em=GD3BJ opt freq int=superfine",
    "단일점 + 파동함수 출력 (Multiwfn 용)":
        "#p B3LYP/6-311G(d,p) em=GD3BJ out=wfn int=ultrafine",
    "TDDFT (hole-electron 해석용)":
        "#p CAM-B3LYP/6-31G(d) em=GD3BJ td(nstates=10,singlets) "
        "IOp(9/40=4) int=ultrafine",
    "TDDFT (ωB97X-D)":
        "#p wB97XD/6-31G(d) td(nstates=10,singlets) IOp(9/40=4) int=ultrafine",
    "이면각 relaxed scan":
        "#p B3LYP/6-31G(d) em=GD3BJ opt=modredundant int=ultrafine",
    "용매효과 (SMD, 물)":
        "#p B3LYP/6-311G(d,p) em=GD3BJ opt freq scrf=(smd,solvent=water) "
        "int=ultrafine",
    "양이온 (재구성에너지용, charge +1)":
        "#p UB3LYP/6-31G(d) em=GD3BJ opt freq int=ultrafine",
    "NMR 화학이동 (GIAO)":
        "#p B3LYP/6-311+G(2d,p) nmr=giao int=ultrafine",
}


@dataclass
class GaussianJob:
    name: str
    route: str
    coords: str = ""
    charge: int = 0
    mult: int = 1
    nproc: int = 16
    mem: str = "32GB"
    chk: str | None = None          # None -> f"{name}.chk"
    title: str = ""
    tail: str = ""                  # modredundant / genecp / wfn 파일명 등
    read_geom: bool = False         # geom=check guess=read 로 chk 이어받기

    def render(self) -> str:
        chk = self.chk or f"{self.name}.chk"
        route = self.route.strip()
        if self.read_geom and "geom=check" not in route:
            route += " geom=check guess=read"

        head = [f"%chk={chk}",
                f"%nprocshared={self.nproc}",
                f"%mem={self.mem}"]
        body = [route, "",
                self.title.strip() or self.name, "",
                f"{self.charge} {self.mult}"]
        if not self.read_geom:
            body.append(self.coords.rstrip())
        body.append("")
        if self.tail:
            body.append(self.tail.rstrip())
            body.append("")
        return "\n".join(head + body) + "\n"

    def write(self, outdir: str) -> str:
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, f"{self.name}.gjf")
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(self.render())
        return path


def field_units(e_mv_per_m: float) -> int:
    """MV/m -> Gaussian `Field=X+N` 의 정수 N.

    Gaussian 의 Field 는 N x 10^-4 a.u. 단위이며 1 a.u. = 5.1422e11 V/m.
    """
    return int(round(e_mv_per_m * 1e6 / AU_FIELD_V_PER_M * 1e4))


def field_from_units(n: int) -> float:
    """역변환: Gaussian Field 정수 N -> MV/m."""
    return n * 1e-4 * AU_FIELD_V_PER_M / 1e6


def make_dihedral_tail(atoms: Sequence[int], npoints: int, step: float) -> str:
    a, b, c, d = atoms
    return f"D {a} {b} {c} {d} S {npoints} {step:.1f}"


def make_field_series(base_name: str, route: str, chk: str,
                      fields_mv_per_m: Sequence[float], axis: str = "X",
                      nproc: int = 16, mem: str = "32GB") -> list[GaussianJob]:
    """동일 chk 를 이어받아 전기장만 바꾸는 잡 시리즈."""
    jobs = []
    for e in fields_mv_per_m:
        n = field_units(e)
        r = route if n == 0 else f"{route} Field={axis}{'+' if n >= 0 else ''}{n}"
        jobs.append(GaussianJob(
            name=f"{base_name}_F{int(round(e))}", route=r, chk=chk,
            nproc=nproc, mem=mem, read_geom=True,
            title=f"{base_name} @ {e:.0f} MV/m along {axis}",
        ))
    return jobs


# ---------------------------------------------------------------------------
# 6. 올리고머 빌더 (사슬길이 시리즈)
# ---------------------------------------------------------------------------

def _dummy_indices(mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]


def _join(mol_a, mol_b):
    Chem = _chem()
    combo = Chem.RWMol(Chem.CombineMols(mol_a, mol_b))
    offset = mol_a.GetNumAtoms()
    d_a = _dummy_indices(mol_a)
    d_b = [i + offset for i in _dummy_indices(mol_b)]
    if not d_a or not d_b:
        raise ValueError("반복단위 양 끝에 더미원자가 필요합니다.")
    tail, head = d_a[-1], d_b[0]
    nb_a = combo.GetAtomWithIdx(tail).GetNeighbors()[0].GetIdx()
    nb_b = combo.GetAtomWithIdx(head).GetNeighbors()[0].GetIdx()
    combo.AddBond(nb_a, nb_b, Chem.BondType.SINGLE)
    for idx in sorted([tail, head], reverse=True):
        combo.RemoveAtom(idx)
    out = combo.GetMol()
    Chem.SanitizeMol(out)
    return out


CAP_ATOMS = {"H": 1, "F": 9, "Cl": 17, "Br": 35, "CH3": 6}


def build_oligomer(repeat_unit, n: int, cap: str = "H"):
    """더미원자(부착점) 2개를 가진 반복단위로부터 n-mer 생성 후 말단 캡핑.

    ChemDraw 에서 부착점을 그리는 방법
    ---------------------------------
      * R 그룹 라벨(R1, R2) 을 양 끝에 붙이거나
      * Attachment point 도구를 사용하거나
      * SMILES 로 직접 입력할 때는 `[*]` 두 개  예: [*]c1ccc([*])s1
    """
    Chem = _chem()
    ru = Chem.Mol(repeat_unit)
    dummies = _dummy_indices(ru)
    if len(dummies) != 2:
        raise ValueError(
            f"반복단위에 부착점(더미원자)이 정확히 2개 필요합니다. 현재 {len(dummies)}개.\n"
            "예: [*]c1ccc([*])s1"
        )
    if n < 1:
        raise ValueError("n 은 1 이상이어야 합니다.")

    mol = Chem.Mol(ru)
    for _ in range(n - 1):
        mol = _join(mol, Chem.Mol(ru))

    rw = Chem.RWMol(mol)
    for idx in sorted(_dummy_indices(mol), reverse=True):
        atom = rw.GetAtomWithIdx(idx)
        atom.SetAtomicNum(CAP_ATOMS.get(cap, 1))
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


# ---------------------------------------------------------------------------
# 7. Gaussian 로그 파싱
# ---------------------------------------------------------------------------

@dataclass
class LogResult:
    file: str
    name: str = ""
    n: int | None = None
    normal_termination: bool = False
    error_line: str = ""
    n_steps: int = 0
    scf_hartree: float | None = None
    homo_eV: float | None = None
    lumo_eV: float | None = None
    gap_eV: float | None = None
    dipole_debye: float | None = None
    n_imag: int | None = None
    zpe_hartree: float | None = None
    gibbs_hartree: float | None = None
    s1_eV: float | None = None
    s1_nm: float | None = None
    s1_f: float | None = None
    excitations: list[tuple[int, float, float, float]] = field(default_factory=list)
    scan: list[tuple[float, float]] = field(default_factory=list)

    def row(self, cols: Sequence[str]) -> list:
        return [getattr(self, c, None) for c in cols]


TABLE_COLS = ["name", "n", "normal_termination", "n_imag", "scf_hartree",
              "gibbs_hartree", "homo_eV", "lumo_eV", "gap_eV",
              "dipole_debye", "s1_eV", "s1_nm", "s1_f", "error_line"]

_RE_SCF = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
_RE_DIPOLE_TOT = re.compile(r"Tot=\s*(-?\d+\.\d+)")
_RE_EXC = re.compile(
    r"Excited State\s+(\d+):\s+\S+\s+(-?\d+\.\d+)\s*eV\s+(-?\d+\.\d+)\s*nm\s+f=(-?\d+\.\d+)"
)
_RE_ZPE = re.compile(r"Sum of electronic and zero-point Energies=\s*(-?\d+\.\d+)")
_RE_GIBBS = re.compile(r"Sum of electronic and thermal Free Energies=\s*(-?\d+\.\d+)")
_RE_N_TAG = re.compile(r"(?:^|[_-])n(\d+)(?:[_-]|$)", re.IGNORECASE)
_ERR_PAT = ("Error termination", "galloc", "Convergence failure",
            "The combination of multiplicity", "l9999.exe", "l502.exe",
            "l103.exe", "Optimization stopped")


def parse_log(path: str) -> LogResult:
    """Gaussian .log / .out 파싱. 실패한 잡도 사유를 담아 반환."""
    res = LogResult(file=path,
                    name=os.path.splitext(os.path.basename(path))[0])
    m = _RE_N_TAG.search(res.name)
    if m:
        res.n = int(m.group(1))

    occ: list[float] = []
    virt: list[float] = []
    freqs: list[float] = []
    scan_pairs: list[tuple[float, float]] = []
    pending_angle: float | None = None
    in_orbitals = False

    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if "Normal termination" in line:
                res.normal_termination = True
            elif any(p in line for p in _ERR_PAT):
                res.error_line = line.strip()[:180]

            if line.startswith(" SCF Done:") or "SCF Done:" in line:
                mm = _RE_SCF.search(line)
                if mm:
                    res.scf_hartree = float(mm.group(1))
                    res.n_steps += 1
                    if pending_angle is not None:
                        scan_pairs.append((pending_angle, res.scf_hartree))
                        pending_angle = None

            if "Alpha  occ. eigenvalues" in line:
                if not in_orbitals:
                    occ, virt = [], []       # 마지막 SCF 것만 유지
                    in_orbitals = True
                occ.extend(float(x) for x in line.split("--", 1)[1].split())
            elif "Alpha virt. eigenvalues" in line:
                virt.extend(float(x) for x in line.split("--", 1)[1].split())
            else:
                in_orbitals = False

            if "Dipole moment" in line and "field-independent" in line:
                nxt = fh.readline()
                mm = _RE_DIPOLE_TOT.search(nxt)
                if mm:
                    res.dipole_debye = float(mm.group(1))

            if "Frequencies --" in line:
                freqs.extend(float(x) for x in line.split("--", 1)[1].split())

            mm = _RE_ZPE.search(line)
            if mm:
                res.zpe_hartree = float(mm.group(1))
            mm = _RE_GIBBS.search(line)
            if mm:
                res.gibbs_hartree = float(mm.group(1))

            mm = _RE_EXC.search(line)
            if mm:
                res.excitations.append((int(mm.group(1)), float(mm.group(2)),
                                        float(mm.group(3)), float(mm.group(4))))

            # modredundant 스캔 좌표 리포트
            if "Scan" in line and ("!" in line or "D(" in line):
                for tok in line.replace("!", " ").split():
                    try:
                        val = float(tok)
                    except ValueError:
                        continue
                    if -400.0 <= val <= 400.0:
                        pending_angle = val
                        break

    if occ:
        res.homo_eV = occ[-1] * HARTREE_EV
    if virt:
        res.lumo_eV = virt[0] * HARTREE_EV
    if res.homo_eV is not None and res.lumo_eV is not None:
        res.gap_eV = res.lumo_eV - res.homo_eV
    if freqs:
        res.n_imag = sum(1 for f in freqs if f < 0)
    if res.excitations:
        bright = [e for e in res.excitations if e[3] > 0.01] or res.excitations
        _, res.s1_eV, res.s1_nm, res.s1_f = bright[0]
    res.scan = scan_pairs
    return res


def parse_many(paths_or_patterns: Iterable[str]) -> list[LogResult]:
    files: list[str] = []
    for p in paths_or_patterns:
        if any(ch in p for ch in "*?["):
            files.extend(sorted(glob.glob(p)))
        elif os.path.isdir(p):
            for ext in ("log", "out"):
                files.extend(sorted(glob.glob(os.path.join(p, f"*.{ext}"))))
        else:
            files.append(p)
    seen, uniq = set(), []
    for f in files:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    return [parse_log(f) for f in uniq]


def write_results_csv(results: Sequence[LogResult], path: str,
                      cols: Sequence[str] = TABLE_COLS) -> str:
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in results:
            w.writerow(["" if v is None else v for v in r.row(cols)])
    return path


def write_scan_csv(result: LogResult, path: str) -> str:
    if not result.scan:
        raise ValueError(f"{result.name} 에 스캔 데이터가 없습니다.")
    e0 = min(e for _, e in result.scan)
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["dihedral_deg", "E_hartree", "rel_E_kcal_per_mol"])
        for ang, e in result.scan:
            w.writerow([ang, f"{e:.8f}", f"{(e - e0) * HARTREE_KCAL:.4f}"])
    return path


def boltzmann_weights(scan: Sequence[tuple[float, float]],
                      temperature_K: float = 298.15) -> list[tuple[float, float]]:
    """스캔 PES -> Boltzmann 형태 확률분포 P(θ)."""
    if not scan:
        return []
    RT = 0.001987204259 * temperature_K       # kcal/mol
    e0 = min(e for _, e in scan)
    raw = [(a, math.exp(-(e - e0) * HARTREE_KCAL / RT)) for a, e in scan]
    z = sum(w for _, w in raw) or 1.0
    return [(a, w / z) for a, w in raw]


# ---------------------------------------------------------------------------
# 8. 사슬길이 외삽
# ---------------------------------------------------------------------------

def fit_linear_1overn(ns: Sequence[float], ys: Sequence[float]) -> dict:
    """y(n) = y_inf + a/n  (관행적 1차 외삽)."""
    xs = [1.0 / n for n in ns]
    k = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    den = k * sxx - sx * sx
    if abs(den) < 1e-14:
        return {"model": "linear_1/n", "error": "축퇴된 n 값"}
    a = (k * sxy - sx * sy) / den
    b = (sy - a * sx) / k
    ybar = sy / k
    ss_res = sum((y - (a * x + b)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    return {"model": "linear_1/n", "y_inf": b, "slope_a": a,
            "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")}


def fit_kuhn(ns: Sequence[float], ys: Sequence[float]) -> dict:
    """Kuhn 식: y(n) = y_inf + (y_1 - y_inf)*exp(-a*(n-1)).

    포화 거동을 재현하며 유효공액길이(ECL)를 함께 준다.
    SciPy 가 있으면 curve_fit, 없으면 격자탐색 + 선형 최소자승.
    """
    xs = [float(n) for n in ns]
    ys = [float(y) for y in ys]
    ybar = sum(ys) / len(ys)
    ss_tot = sum((y - ybar) ** 2 for y in ys)

    def eval_fit(y_inf, y_1, a):
        pred = [y_inf + (y_1 - y_inf) * math.exp(-a * (x - 1.0)) for x in xs]
        return sum((y - p) ** 2 for y, p in zip(ys, pred))

    best = None
    try:
        import numpy as np
        from scipy.optimize import curve_fit

        def kuhn(n, y_inf, y_1, a):
            return y_inf + (y_1 - y_inf) * np.exp(-a * (n - 1.0))

        p0 = [min(ys), ys[0], 0.35]
        popt, _ = curve_fit(kuhn, np.array(xs), np.array(ys), p0=p0, maxfev=40000)
        best = (float(popt[0]), float(popt[1]), float(popt[2]))
    except Exception:
        # SciPy 없거나 실패 -> a 를 격자로 훑고 나머지는 선형해
        for a in [0.02 * i for i in range(1, 150)]:
            basis = [math.exp(-a * (x - 1.0)) for x in xs]
            k = len(xs)
            sb = sum(basis)
            sbb = sum(b * b for b in basis)
            sy = sum(ys)
            sby = sum(b * y for b, y in zip(basis, ys))
            den = k * sbb - sb * sb
            if abs(den) < 1e-14:
                continue
            slope = (k * sby - sb * sy) / den      # = y_1 - y_inf
            inter = (sy - slope * sb) / k          # = y_inf
            cand = (inter, inter + slope, a)
            if best is None or eval_fit(*cand) < eval_fit(*best):
                best = cand

    if best is None:
        return {"model": "kuhn", "error": "피팅 실패"}

    y_inf, y_1, a = best
    ss_res = eval_fit(y_inf, y_1, a)
    ecl = 1.0 + math.log(20.0) / a if a > 1e-6 else float("nan")
    return {"model": "kuhn", "y_inf": y_inf, "y_1": y_1, "decay_a": a,
            "ECL_units_95pct": ecl,
            "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")}


def extrapolate(ns: Sequence[float], ys: Sequence[float]) -> dict:
    if len(ns) < 3:
        raise ValueError("외삽에는 최소 3개(권장 5개 이상)의 사슬길이가 필요합니다.")
    order = sorted(range(len(ns)), key=lambda i: ns[i])
    ns = [ns[i] for i in order]
    ys = [ys[i] for i in order]
    return {"n": list(ns), "y": list(ys),
            "linear": fit_linear_1overn(ns, ys),
            "kuhn": fit_kuhn(ns, ys)}


def plot_extrapolation(fit: dict, prop_label: str, path: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ns, ys = fit["n"], fit["y"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ax = axes[0]
    ax.plot(ns, ys, "o", ms=7, color="#1f4e79", label="DFT")
    k = fit.get("kuhn", {})
    if "y_inf" in k:
        xs = [min(ns) + (max(ns) * 2.5 - min(ns)) * i / 199 for i in range(200)]
        ys_fit = [k["y_inf"] + (k["y_1"] - k["y_inf"]) * math.exp(-k["decay_a"] * (x - 1))
                  for x in xs]
        ax.plot(xs, ys_fit, "-", lw=1.8, color="#c00000",
                label=f"Kuhn  y($\\infty$)={k['y_inf']:.3f}")
        ax.axhline(k["y_inf"], ls=":", lw=1.1, color="#c00000")
    ax.set_xlabel("repeat units, n")
    ax.set_ylabel(prop_label)
    ax.set_title("chain-length dependence")
    ax.legend(fontsize=8)

    ax = axes[1]
    inv = [1.0 / n for n in ns]
    ax.plot(inv, ys, "o", ms=7, color="#1f4e79", label="DFT")
    lin = fit.get("linear", {})
    if "y_inf" in lin:
        xs = [max(inv) * 1.05 * i / 49 for i in range(50)]
        ax.plot(xs, [lin["y_inf"] + lin["slope_a"] * x for x in xs], "-",
                lw=1.8, color="#2e7d32",
                label=f"1/n  y($\\infty$)={lin['y_inf']:.3f}, "
                      f"R$^2$={lin['r2']:.4f}")
    ax.set_xlabel("1 / n")
    ax.set_ylabel(prop_label)
    ax.set_title("1/n extrapolation")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def plot_scan(scan: Sequence[tuple[float, float]], path: str,
              temperature_K: float = 298.15, title: str = "") -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    e0 = min(e for _, e in scan)
    ang = [a for a, _ in scan]
    rel = [(e - e0) * HARTREE_KCAL for _, e in scan]
    prob = boltzmann_weights(scan, temperature_K)

    fig, ax1 = plt.subplots(figsize=(6.6, 4.2))
    ax1.plot(ang, rel, "o-", ms=4, lw=1.6, color="#1f4e79")
    ax1.set_xlabel("dihedral angle (deg)")
    ax1.set_ylabel("relative energy (kcal/mol)", color="#1f4e79")
    ax1.tick_params(axis="y", labelcolor="#1f4e79")

    ax2 = ax1.twinx()
    ax2.fill_between([a for a, _ in prob], 0, [p for _, p in prob],
                     alpha=0.25, color="#c00000")
    ax2.set_ylabel(f"P($\\theta$) @ {temperature_K:.0f} K", color="#c00000")
    ax2.tick_params(axis="y", labelcolor="#c00000")

    if title:
        ax1.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path
