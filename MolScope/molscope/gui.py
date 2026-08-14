# -*- coding: utf-8 -*-
"""
molscope.gui — Tkinter 기반 데스크톱 UI

워크플로 (탭 순서 그대로)
  1) 구조    : ChemDraw 파일(.cdx/.cdxml) 또는 .mol/.sdf 불러오기, SMILES 붙여넣기
  2) 분석    : 기술자 + 공액 분석 + 이면각 후보 + 2D 미리보기
  3) 계산 생성: 3D 좌표 생성 -> Gaussian 입력(.gjf) — 최적화/TDDFT/스캔/전기장
  4) 올리고머 : 반복단위(부착점 2개) -> n-mer 시리즈 일괄 생성
  5) 결과    : Gaussian 로그 일괄 파싱 -> 표/CSV, 스캔 PES 그림
  6) 외삽    : 사슬길이 시리즈 -> 1/n + Kuhn 외삽, 유효공액길이
"""

from __future__ import annotations

import os
import sys
import tempfile
import threading
import traceback
import webbrowser
from tkinter import (
    Tk, ttk, StringVar, IntVar, DoubleVar, BooleanVar, Text, Menu,
    filedialog, messagebox, END, WORD, VERTICAL, HORIZONTAL,
)

from . import __version__, __app_name__
from . import core

PAD = {"padx": 6, "pady": 4}


def _busy(widget, fn, done=None):
    """긴 작업을 백그라운드 스레드에서 실행하고 완료 시 UI 스레드로 복귀."""
    def run():
        try:
            result = fn()
            err = None
        except Exception as exc:  # noqa: BLE001
            result, err = None, exc
            traceback.print_exc()
        def finish():
            widget.config(cursor="")
            if err is not None:
                messagebox.showerror(__app_name__, f"오류가 발생했습니다:\n{err}")
            elif done is not None:
                done(result)
        widget.after(0, finish)
    widget.config(cursor="watch")
    threading.Thread(target=run, daemon=True).start()


class App(ttk.Frame):
    def __init__(self, master: Tk):
        super().__init__(master)
        master.title(f"{__app_name__} v{__version__} — ChemDraw 구조 분석/계산 준비 도구")
        master.geometry("1180x780")
        master.minsize(980, 640)
        self.pack(fill="both", expand=True)

        self.records: list[core.MolRecord] = []
        self.embedded: dict[int, core.Embed3DResult] = {}   # record index -> 3D
        self.log_results: list[core.LogResult] = []
        self._preview_img = None
        self._tmpdir = tempfile.mkdtemp(prefix="molscope_")

        self._build_menu(master)
        self._build_statusbar()
        self._build_tabs()
        self._check_rdkit()

    # ------------------------------------------------------------------ UI 뼈대
    def _build_menu(self, master: Tk):
        menubar = Menu(master)
        m_file = Menu(menubar, tearoff=0)
        m_file.add_command(label="구조 파일 열기...  (ChemDraw .cdx/.cdxml)",
                           command=self.open_files, accelerator="Ctrl+O")
        m_file.add_separator()
        m_file.add_command(label="종료", command=master.destroy)
        menubar.add_cascade(label="파일", menu=m_file)

        m_help = Menu(menubar, tearoff=0)
        m_help.add_command(label="ChemDraw에서 내보내는 방법", command=self._help_chemdraw)
        m_help.add_command(label="Multiwfn 홈페이지", 
                           command=lambda: webbrowser.open("http://sobereva.com/multiwfn"))
        m_help.add_command(label="정보", command=self._about)
        menubar.add_cascade(label="도움말", menu=m_help)

        master.config(menu=menubar)
        master.bind("<Control-o>", lambda e: self.open_files())

    def _build_statusbar(self):
        self.status = StringVar(value="준비됨")
        bar = ttk.Label(self, textvariable=self.status, relief="sunken", anchor="w")
        bar.pack(side="bottom", fill="x")

    def _set_status(self, msg: str):
        self.status.set(msg)

    def _build_tabs(self):
        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True)
        self.tab_load = ttk.Frame(self.nb)
        self.tab_analyze = ttk.Frame(self.nb)
        self.tab_gjf = ttk.Frame(self.nb)
        self.tab_oligo = ttk.Frame(self.nb)
        self.tab_logs = ttk.Frame(self.nb)
        self.tab_fit = ttk.Frame(self.nb)
        self.nb.add(self.tab_load, text=" 1. 구조 불러오기 ")
        self.nb.add(self.tab_analyze, text=" 2. 분석 ")
        self.nb.add(self.tab_gjf, text=" 3. Gaussian 입력 생성 ")
        self.nb.add(self.tab_oligo, text=" 4. 올리고머 시리즈 ")
        self.nb.add(self.tab_logs, text=" 5. 계산결과 파싱 ")
        self.nb.add(self.tab_fit, text=" 6. 사슬길이 외삽 ")
        self._tab1()
        self._tab2()
        self._tab3()
        self._tab4()
        self._tab5()
        self._tab6()

    def _check_rdkit(self):
        ok, ver = core.rdkit_available()
        if not ok:
            messagebox.showerror(__app_name__, core._RDKIT_ERR)
            self._set_status("RDKit 로드 실패 — 기능이 제한됩니다")
            return
        cdx = "지원" if core.chemdraw_cdx_supported() else "미지원(→ .cdxml 사용)"
        self._set_status(f"RDKit {ver} | ChemDraw 네이티브(.cdx): {cdx}")

    # ------------------------------------------------------------ 탭1: 불러오기
    def _tab1(self):
        f = self.tab_load
        top = ttk.LabelFrame(f, text="ChemDraw / 구조 파일")
        top.pack(fill="x", **PAD)
        ttk.Button(top, text="파일 열기...", command=self.open_files).pack(side="left", **PAD)
        ttk.Label(top, text="지원: .cdx (ChemDraw 네이티브)  .cdxml  .mol  .sdf  .smi").pack(
            side="left", padx=8)

        mid = ttk.LabelFrame(
            f, text="SMILES 붙여넣기  (ChemDraw: 구조 선택 → Edit → Copy As → SMILES)")
        mid.pack(fill="x", **PAD)
        self.txt_smiles = Text(mid, height=4, wrap=WORD)
        self.txt_smiles.pack(fill="x", padx=6, pady=2)
        row = ttk.Frame(mid)
        row.pack(fill="x")
        ttk.Button(row, text="붙여넣은 SMILES 추가", command=self.add_smiles).pack(
            side="left", **PAD)
        ttk.Label(row, text="여러 줄 가능, 형식: SMILES [이름]").pack(side="left")

        lst = ttk.LabelFrame(f, text="불러온 구조 목록")
        lst.pack(fill="both", expand=True, **PAD)
        cols = ("name", "formula", "heavy", "src", "fmt", "note")
        heads = ("이름", "화학식", "중원자", "출처", "형식", "비고")
        widths = (170, 130, 60, 330, 70, 220)
        self.tree_mols = ttk.Treeview(lst, columns=cols, show="headings", height=12)
        for c, h, w in zip(cols, heads, widths):
            self.tree_mols.heading(c, text=h)
            self.tree_mols.column(c, width=w, anchor="w")
        vs = ttk.Scrollbar(lst, orient=VERTICAL, command=self.tree_mols.yview)
        self.tree_mols.configure(yscrollcommand=vs.set)
        self.tree_mols.pack(side="left", fill="both", expand=True)
        vs.pack(side="right", fill="y")

        btns = ttk.Frame(f)
        btns.pack(fill="x", **PAD)
        ttk.Button(btns, text="선택 삭제", command=self.remove_selected).pack(side="left", **PAD)
        ttk.Button(btns, text="전체 비우기", command=self.clear_records).pack(side="left", **PAD)

    def open_files(self):
        paths = filedialog.askopenfilenames(title="구조 파일 선택",
                                            filetypes=core.SUPPORTED_INPUT)
        if not paths:
            return

        def work():
            recs, errors = [], []
            for p in paths:
                try:
                    recs.extend(core.load_any(p))
                except Exception as exc:  # noqa: BLE001
                    errors.append(f"{os.path.basename(p)}: {exc}")
            return recs, errors

        def done(res):
            recs, errors = res
            self.records.extend(recs)
            self._refresh_mol_list()
            self._set_status(f"{len(recs)}개 구조 불러옴 (총 {len(self.records)}개)")
            if errors:
                messagebox.showwarning(__app_name__, "일부 파일 실패:\n" + "\n".join(errors))

        _busy(self, work, done)

    def add_smiles(self):
        text = self.txt_smiles.get("1.0", END)
        if not text.strip():
            return
        try:
            recs = core.load_smiles_text(text)
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror(__app_name__, str(exc))
            return
        bad = [r for r in recs if r.mol is None]
        self.records.extend(r for r in recs if r.mol is not None)
        self._refresh_mol_list()
        self.txt_smiles.delete("1.0", END)
        msg = f"{len(recs) - len(bad)}개 추가"
        if bad:
            msg += f", {len(bad)}개 파싱 실패"
        self._set_status(msg)

    def remove_selected(self):
        sel = [int(i) for i in self.tree_mols.selection()]
        for i in sorted(sel, reverse=True):
            self.records.pop(i)
            self.embedded.pop(i, None)
        # 인덱스 재정렬: 임베딩 캐시는 단순화를 위해 초기화
        self.embedded.clear()
        self._refresh_mol_list()

    def clear_records(self):
        self.records.clear()
        self.embedded.clear()
        self._refresh_mol_list()

    def _refresh_mol_list(self):
        for t in (self.tree_mols,):
            t.delete(*t.get_children())
        from rdkit.Chem import rdMolDescriptors
        for i, r in enumerate(self.records):
            formula = rdMolDescriptors.CalcMolFormula(r.mol) if r.mol is not None else "-"
            self.tree_mols.insert("", END, iid=str(i), values=(
                r.name, formula, r.n_atoms_heavy if r.mol is not None else "-",
                os.path.basename(r.source) if r.source else "-",
                r.origin_format, r.note))
        self._refresh_pickers()

    def _selected_record(self) -> tuple[int, core.MolRecord] | None:
        # 어느 탭에서든 현재 선택된 분자
        cmb = getattr(self, "cmb_target", None)
        if cmb is not None and cmb.current() >= 0 and cmb.current() < len(self.records):
            return cmb.current(), self.records[cmb.current()]
        if self.records:
            return 0, self.records[0]
        return None

    def _refresh_pickers(self):
        names = [f"{i + 1}. {r.name}" for i, r in enumerate(self.records)]
        for attr in ("cmb_target", "cmb_gjf_target", "cmb_ru_target"):
            cmb = getattr(self, attr, None)
            if cmb is not None:
                cur = cmb.current()
                cmb["values"] = names
                if names:
                    cmb.current(min(max(cur, 0), len(names) - 1))

    # ------------------------------------------------------------ 탭2: 분석
    def _tab2(self):
        f = self.tab_analyze
        top = ttk.Frame(f)
        top.pack(fill="x", **PAD)
        ttk.Label(top, text="분자 선택:").pack(side="left")
        self.cmb_target = ttk.Combobox(top, state="readonly", width=42)
        self.cmb_target.pack(side="left", padx=6)
        ttk.Button(top, text="분석 실행", command=self.run_analysis).pack(side="left", **PAD)
        ttk.Button(top, text="결과를 CSV로 저장 (전체 분자)",
                   command=self.export_analysis_csv).pack(side="left", **PAD)

        body = ttk.Panedwindow(f, orient=HORIZONTAL)
        body.pack(fill="both", expand=True, **PAD)

        left = ttk.LabelFrame(body, text="기술자 / 공액 분석")
        self.txt_analysis = Text(left, wrap=WORD, width=58)
        self.txt_analysis.pack(fill="both", expand=True, padx=4, pady=4)
        body.add(left, weight=1)

        right = ttk.LabelFrame(body, text="2D 미리보기")
        self.lbl_preview = ttk.Label(right, text="(분석 실행 시 표시)", anchor="center")
        self.lbl_preview.pack(fill="both", expand=True, padx=4, pady=4)
        body.add(right, weight=1)

    def run_analysis(self):
        sel = self._selected_record()
        if sel is None:
            messagebox.showinfo(__app_name__, "먼저 구조를 불러오세요 (탭 1).")
            return
        idx, rec = sel

        def work():
            d = core.describe(rec.mol)
            c = core.conjugation_analysis(rec.mol)
            dih = core.find_rotatable_dihedrals(rec.mol, max_out=12)
            png = os.path.join(self._tmpdir, f"prev_{idx}.png")
            core.render_2d_png(rec.mol, png, legend=rec.name)
            return d, c, dih, png

        def done(res):
            d, c, dih, png = res
            t = self.txt_analysis
            t.delete("1.0", END)
            t.insert(END, f"◆ {rec.name}\n  SMILES: {rec.smiles}\n\n")
            t.insert(END, "── 기본 기술자 ─────────────────────\n")
            for k, v in d.items():
                t.insert(END, f"  {k:24s} {v}\n")
            t.insert(END, "\n── 공액(π) 분석 ────────────────────\n")
            for k, v in c.items():
                t.insert(END, f"  {k:28s} {v}\n")
            t.insert(END, "\n── 이면각 스캔 후보 (Gaussian 1-based) ─\n")
            if not dih:
                t.insert(END, "  (회전 가능한 이면각 없음)\n")
            for x in dih:
                tag = "★고리-고리" if x["aryl_aryl"] else ("공액" if x["conjugated"] else "")
                t.insert(END, f"  D {x['atoms'][0]:>3} {x['atoms'][1]:>3} "
                              f"{x['atoms'][2]:>3} {x['atoms'][3]:>3}   {x['label']}  {tag}\n")
            t.insert(END, "\n※ '고리-고리' 이면각이 공액 경로/평면성을 지배하는 "
                          "핵심 좌표입니다.\n")
            try:
                from PIL import Image, ImageTk
                img = Image.open(png)
                img.thumbnail((520, 420))
                self._preview_img = ImageTk.PhotoImage(img)
                self.lbl_preview.configure(image=self._preview_img, text="")
            except Exception:
                self.lbl_preview.configure(text="(미리보기 렌더링 실패 — Pillow 필요)")
            self._set_status(f"{rec.name} 분석 완료")

        _busy(self, work, done)

    def export_analysis_csv(self):
        if not self.records:
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")],
                                            initialfile="molscope_analysis.csv")
        if not path:
            return

        def work():
            import csv as _csv
            rows = []
            keys: list[str] = []
            for r in self.records:
                d = {"name": r.name, "smiles": r.smiles, "source": r.source}
                d.update(core.describe(r.mol))
                c = core.conjugation_analysis(r.mol)
                c.pop("system_sizes", None)
                d.update(c)
                rows.append(d)
                for k in d:
                    if k not in keys:
                        keys.append(k)
            with open(path, "w", newline="", encoding="utf-8-sig") as fh:
                w = _csv.DictWriter(fh, fieldnames=keys)
                w.writeheader()
                w.writerows(rows)
            return path

        _busy(self, work, lambda p: self._set_status(f"저장됨: {p}"))

    # ------------------------------------------------------------ 탭3: Gaussian
    def _tab3(self):
        f = self.tab_gjf
        top = ttk.Frame(f)
        top.pack(fill="x", **PAD)
        ttk.Label(top, text="분자 선택:").pack(side="left")
        self.cmb_gjf_target = ttk.Combobox(top, state="readonly", width=42)
        self.cmb_gjf_target.pack(side="left", padx=6)

        opt = ttk.LabelFrame(f, text="계산 설정")
        opt.pack(fill="x", **PAD)

        r0 = ttk.Frame(opt); r0.pack(fill="x", **PAD)
        ttk.Label(r0, text="Route 프리셋:").pack(side="left")
        self.cmb_preset = ttk.Combobox(r0, state="readonly", width=44,
                                       values=list(core.ROUTE_PRESETS))
        self.cmb_preset.current(0)
        self.cmb_preset.pack(side="left", padx=6)
        self.cmb_preset.bind("<<ComboboxSelected>>", self._preset_changed)

        r1 = ttk.Frame(opt); r1.pack(fill="x", **PAD)
        ttk.Label(r1, text="Route:").pack(side="left")
        self.var_route = StringVar(value=core.ROUTE_PRESETS[list(core.ROUTE_PRESETS)[0]])
        ttk.Entry(r1, textvariable=self.var_route, width=100).pack(
            side="left", padx=6, fill="x", expand=True)

        r2 = ttk.Frame(opt); r2.pack(fill="x", **PAD)
        self.var_charge = IntVar(value=0)
        self.var_mult = IntVar(value=1)
        self.var_nproc = IntVar(value=16)
        self.var_mem = StringVar(value="32GB")
        self.var_nconf = IntVar(value=20)
        for lbl, var, w in (("전하", self.var_charge, 4), ("다중도", self.var_mult, 4),
                            ("nproc", self.var_nproc, 5), ("메모리", self.var_mem, 7),
                            ("컨포머 수", self.var_nconf, 5)):
            ttk.Label(r2, text=lbl).pack(side="left", padx=(10, 2))
            ttk.Entry(r2, textvariable=var, width=w).pack(side="left")

        scan = ttk.LabelFrame(f, text="이면각 스캔 (route 에 opt=modredundant 선택 시)")
        scan.pack(fill="x", **PAD)
        rs = ttk.Frame(scan); rs.pack(fill="x", **PAD)
        self.var_scan_atoms = StringVar(value="")
        self.var_scan_np = IntVar(value=36)
        self.var_scan_step = DoubleVar(value=10.0)
        ttk.Label(rs, text="이면각 원자 4개(공백 구분, 비우면 자동):").pack(side="left")
        ttk.Entry(rs, textvariable=self.var_scan_atoms, width=18).pack(side="left", padx=4)
        ttk.Label(rs, text="점 수").pack(side="left", padx=(10, 2))
        ttk.Entry(rs, textvariable=self.var_scan_np, width=5).pack(side="left")
        ttk.Label(rs, text="간격(°)").pack(side="left", padx=(10, 2))
        ttk.Entry(rs, textvariable=self.var_scan_step, width=6).pack(side="left")

        fld = ttk.LabelFrame(f, text="외부 전기장 시리즈 (TDDFT 프리셋과 조합, chk 이어받기)")
        fld.pack(fill="x", **PAD)
        rf_ = ttk.Frame(fld); rf_.pack(fill="x", **PAD)
        self.var_use_field = BooleanVar(value=False)
        self.var_fields = StringVar(value="0 400 800 1200 1600")
        self.var_axis = StringVar(value="X")
        ttk.Checkbutton(rf_, text="전기장 시리즈 생성", variable=self.var_use_field).pack(side="left")
        ttk.Label(rf_, text="전기장(MV/m):").pack(side="left", padx=(10, 2))
        ttk.Entry(rf_, textvariable=self.var_fields, width=28).pack(side="left")
        ttk.Label(rf_, text="축:").pack(side="left", padx=(10, 2))
        ttk.Combobox(rf_, textvariable=self.var_axis, values=["X", "Y", "Z"],
                     width=3, state="readonly").pack(side="left")
        ttk.Label(rf_, text="  ※ Field=축+N, N×10⁻⁴ a.u. 자동 변환").pack(side="left")

        run = ttk.Frame(f); run.pack(fill="x", **PAD)
        ttk.Button(run, text="3D 생성 + .gjf 저장...", command=self.generate_gjf).pack(
            side="left", **PAD)
        ttk.Button(run, text="3D 좌표만 저장 (.xyz/.mol)...", command=self.save_3d).pack(
            side="left", **PAD)

        self.txt_gjf_log = Text(f, height=10, wrap=WORD)
        self.txt_gjf_log.pack(fill="both", expand=True, **PAD)

    def _preset_changed(self, _evt=None):
        name = self.cmb_preset.get()
        self.var_route.set(core.ROUTE_PRESETS.get(name, self.var_route.get()))
        if "charge +1" in name:
            self.var_charge.set(1)
            self.var_mult.set(2)

    def _pick_gjf_record(self):
        i = self.cmb_gjf_target.current()
        if i < 0 or i >= len(self.records):
            if not self.records:
                messagebox.showinfo(__app_name__, "먼저 구조를 불러오세요 (탭 1).")
                return None
            i = 0
        return i, self.records[i]

    def _get_embedded(self, idx: int, rec: core.MolRecord) -> core.Embed3DResult:
        if idx not in self.embedded:
            self.embedded[idx] = core.embed3d(rec.mol, n_conf=self.var_nconf.get())
        return self.embedded[idx]

    def generate_gjf(self):
        sel = self._pick_gjf_record()
        if sel is None:
            return
        idx, rec = sel
        outdir = filedialog.askdirectory(title=".gjf 저장 폴더 선택")
        if not outdir:
            return
        route = self.var_route.get().strip()
        use_scan = "modredundant" in route
        use_field = self.var_use_field.get()

        def work():
            lines = []
            emb = self._get_embedded(idx, rec)
            lines.append(f"3D 생성: {emb.forcefield}, 컨포머 {emb.n_conformers_tried}개, "
                         f"E={emb.energy:.2f} kcal/mol, 수렴={emb.converged}")

            tail = ""
            if use_scan:
                atoms_str = self.var_scan_atoms.get().split()
                if len(atoms_str) == 4:
                    atoms = tuple(int(a) for a in atoms_str)
                else:
                    cands = core.find_rotatable_dihedrals(emb.mol)
                    if not cands:
                        raise RuntimeError("스캔 가능한 이면각이 없습니다. 원자를 직접 지정하세요.")
                    atoms = cands[0]["atoms"]
                    lines.append(f"이면각 자동 선택: {cands[0]['label']}")
                tail = core.make_dihedral_tail(atoms, self.var_scan_np.get(),
                                               self.var_scan_step.get())
                lines.append(f"modredundant: {tail}")

            base = core.GaussianJob(
                name=rec.name, route=route, coords=core.to_xyz_block(emb.mol),
                charge=self.var_charge.get(), mult=self.var_mult.get(),
                nproc=self.var_nproc.get(), mem=self.var_mem.get(),
                title=f"{rec.name} | from ChemDraw via MolScope", tail=tail)
            p = base.write(outdir)
            lines.append(f"저장: {p}")

            if use_field:
                fields = [float(x) for x in self.var_fields.get().split()]
                jobs = core.make_field_series(
                    rec.name, route, f"{rec.name}.chk", fields,
                    axis=self.var_axis.get(), nproc=self.var_nproc.get(),
                    mem=self.var_mem.get())
                for j in jobs:
                    j.charge, j.mult = self.var_charge.get(), self.var_mult.get()
                    lines.append(f"저장: {j.write(outdir)}   "
                                 f"({j.route.split('Field=')[-1] if 'Field=' in j.route else '0 MV/m'})")
                lines.append("※ 전기장 잡들은 기본 잡의 .chk 를 이어받습니다. "
                             "기본 잡을 먼저 돌리세요.")
            return "\n".join(lines)

        def done(msg):
            self.txt_gjf_log.insert(END, msg + "\n" + "─" * 70 + "\n")
            self.txt_gjf_log.see(END)
            self._set_status("Gaussian 입력 생성 완료")

        _busy(self, work, done)

    def save_3d(self):
        sel = self._pick_gjf_record()
        if sel is None:
            return
        idx, rec = sel
        path = filedialog.asksaveasfilename(
            defaultextension=".xyz", initialfile=f"{rec.name}.xyz",
            filetypes=[("XYZ", "*.xyz"), ("MDL Molfile", "*.mol")])
        if not path:
            return

        def work():
            emb = self._get_embedded(idx, rec)
            if path.lower().endswith(".mol"):
                return core.write_molfile(emb.mol, path)
            return core.write_xyz(emb.mol, path, comment=f"{rec.name} ({emb.forcefield})")

        _busy(self, work, lambda p: self._set_status(f"저장됨: {p}"))

    # ------------------------------------------------------------ 탭4: 올리고머
    def _tab4(self):
        f = self.tab_oligo
        info = ttk.Label(f, justify="left", text=(
            "반복단위 → n-mer 시리즈 일괄 생성.  반복단위에는 부착점(더미원자) 2개가 필요합니다.\n"
            "  · ChemDraw: 양 끝에 attachment point 를 찍거나, Copy As → SMILES 로 [*] 포함 SMILES 복사\n"
            "  · 예시 SMILES:  [*]c1ccc([*])s1   (폴리티오펜 반복단위)"))
        info.pack(fill="x", **PAD)

        top = ttk.Frame(f); top.pack(fill="x", **PAD)
        ttk.Label(top, text="반복단위:").pack(side="left")
        self.cmb_ru_target = ttk.Combobox(top, state="readonly", width=36)
        self.cmb_ru_target.pack(side="left", padx=6)
        ttk.Label(top, text="또는 SMILES 직접:").pack(side="left", padx=(12, 2))
        self.var_ru_smiles = StringVar(value="")
        ttk.Entry(top, textvariable=self.var_ru_smiles, width=32).pack(side="left")

        r = ttk.Frame(f); r.pack(fill="x", **PAD)
        self.var_oligo_ns = StringVar(value="1 2 3 4 5 6 8")
        self.var_oligo_cap = StringVar(value="H")
        self.var_oligo_tag = StringVar(value="oligo")
        self.var_oligo_td = BooleanVar(value=True)
        ttk.Label(r, text="n 목록:").pack(side="left")
        ttk.Entry(r, textvariable=self.var_oligo_ns, width=20).pack(side="left", padx=4)
        ttk.Label(r, text="말단 캡:").pack(side="left", padx=(10, 2))
        ttk.Combobox(r, textvariable=self.var_oligo_cap, values=list(core.CAP_ATOMS),
                     width=5, state="readonly").pack(side="left")
        ttk.Label(r, text="이름 태그:").pack(side="left", padx=(10, 2))
        ttk.Entry(r, textvariable=self.var_oligo_tag, width=14).pack(side="left")
        ttk.Checkbutton(r, text="TDDFT 후속 입력도 생성 (chk 이어받기)",
                        variable=self.var_oligo_td).pack(side="left", padx=10)

        ttk.Button(f, text="시리즈 생성 + .gjf 저장...",
                   command=self.generate_oligomers).pack(anchor="w", **PAD)
        self.txt_oligo_log = Text(f, height=16, wrap=WORD)
        self.txt_oligo_log.pack(fill="both", expand=True, **PAD)

    def generate_oligomers(self):
        smi = self.var_ru_smiles.get().strip()
        ru = None
        if smi:
            from rdkit import Chem
            ru = Chem.MolFromSmiles(smi)
            if ru is None:
                messagebox.showerror(__app_name__, f"SMILES 파싱 실패: {smi}")
                return
            ru_name = self.var_oligo_tag.get()
        else:
            i = self.cmb_ru_target.current()
            if i < 0 or i >= len(self.records):
                messagebox.showinfo(__app_name__,
                                    "반복단위를 선택하거나 SMILES 를 입력하세요.")
                return
            ru = self.records[i].mol
            ru_name = self.var_oligo_tag.get() or self.records[i].name

        try:
            ns = sorted({int(x) for x in self.var_oligo_ns.get().split()})
        except ValueError:
            messagebox.showerror(__app_name__, "n 목록은 정수 공백 구분입니다. 예: 1 2 3 4 6 8")
            return
        outdir = filedialog.askdirectory(title=".gjf 저장 폴더 선택")
        if not outdir:
            return

        route = self.var_route.get().strip()
        make_td = self.var_oligo_td.get()

        def work():
            lines = [f"반복단위: {ru_name}, n = {ns}, 캡 = {self.var_oligo_cap.get()}"]
            for n in ns:
                mol = core.build_oligomer(ru, n, cap=self.var_oligo_cap.get())
                emb = core.embed3d(mol, n_conf=max(6, self.var_nconf.get() // 2))
                name = f"{ru_name}_n{n}"
                job = core.GaussianJob(
                    name=name, route=route, coords=core.to_xyz_block(emb.mol),
                    charge=self.var_charge.get(), mult=self.var_mult.get(),
                    nproc=self.var_nproc.get(), mem=self.var_mem.get(),
                    title=f"{name} | oligomer series via MolScope")
                p = job.write(outdir)
                lines.append(f"  n={n:<2d} 원자 {emb.mol.GetNumAtoms():>4d}  -> {p}")
                if make_td:
                    td = core.GaussianJob(
                        name=f"{name}_td",
                        route=core.ROUTE_PRESETS["TDDFT (hole-electron 해석용)"],
                        chk=f"{name}.chk", read_geom=True,
                        charge=self.var_charge.get(), mult=self.var_mult.get(),
                        nproc=self.var_nproc.get(), mem=self.var_mem.get(),
                        title=f"{name} TDDFT")
                    lines.append(f"        TDDFT -> {td.write(outdir)}")
            lines.append("※ 파일명의 _n숫자 태그를 결과 파싱 탭이 자동 인식해 "
                         "외삽에 사용합니다.")
            return "\n".join(lines)

        def done(msg):
            self.txt_oligo_log.insert(END, msg + "\n" + "─" * 70 + "\n")
            self.txt_oligo_log.see(END)
            self._set_status("올리고머 시리즈 생성 완료")

        _busy(self, work, done)

    # ------------------------------------------------------------ 탭5: 로그
    def _tab5(self):
        f = self.tab_logs
        top = ttk.Frame(f); top.pack(fill="x", **PAD)
        ttk.Button(top, text="로그 파일 선택...", command=self.pick_logs).pack(side="left", **PAD)
        ttk.Button(top, text="폴더 전체 파싱...", command=self.pick_log_dir).pack(side="left", **PAD)
        ttk.Button(top, text="표를 CSV 저장...", command=self.save_results_csv).pack(side="left", **PAD)
        ttk.Button(top, text="선택 스캔 → CSV+그림...", command=self.export_scan).pack(side="left", **PAD)

        cols = ("name", "n", "ok", "imag", "scf", "gibbs", "homo", "lumo", "gap",
                "dip", "s1", "err")
        heads = ("이름", "n", "정상종료", "허수", "SCF (Ha)", "G (Ha)", "HOMO (eV)",
                 "LUMO (eV)", "gap (eV)", "μ (D)", "S1 (eV)", "오류")
        widths = (180, 40, 60, 44, 110, 110, 80, 80, 74, 60, 66, 240)
        self.tree_logs = ttk.Treeview(f, columns=cols, show="headings", height=16)
        for c, h, w in zip(cols, heads, widths):
            self.tree_logs.heading(c, text=h)
            self.tree_logs.column(c, width=w, anchor="w")
        vs = ttk.Scrollbar(f, orient=VERTICAL, command=self.tree_logs.yview)
        self.tree_logs.configure(yscrollcommand=vs.set)
        self.tree_logs.pack(side="left", fill="both", expand=True, **PAD)
        vs.pack(side="right", fill="y")

    def _load_logs(self, paths):
        def work():
            return core.parse_many(paths)

        def done(results):
            self.log_results = results
            t = self.tree_logs
            t.delete(*t.get_children())
            for i, r in enumerate(results):
                fmt = lambda v, n=4: ("" if v is None else f"{v:.{n}f}")
                t.insert("", END, iid=str(i), values=(
                    r.name, r.n if r.n is not None else "",
                    "O" if r.normal_termination else "X",
                    r.n_imag if r.n_imag is not None else "",
                    fmt(r.scf_hartree, 6), fmt(r.gibbs_hartree, 6),
                    fmt(r.homo_eV, 3), fmt(r.lumo_eV, 3), fmt(r.gap_eV, 3),
                    fmt(r.dipole_debye, 2), fmt(r.s1_eV, 3), r.error_line))
            n_ok = sum(1 for r in results if r.normal_termination)
            self._set_status(f"로그 {len(results)}개 파싱 (정상 {n_ok}, 실패 {len(results) - n_ok})")

        _busy(self, work, done)

    def pick_logs(self):
        paths = filedialog.askopenfilenames(
            title="Gaussian 로그 선택",
            filetypes=[("Gaussian log", "*.log *.out"), ("모든 파일", "*.*")])
        if paths:
            self._load_logs(list(paths))

    def pick_log_dir(self):
        d = filedialog.askdirectory(title="로그 폴더 선택")
        if d:
            self._load_logs([d])

    def save_results_csv(self):
        if not self.log_results:
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")],
                                            initialfile="molscope_results.csv")
        if path:
            core.write_results_csv(self.log_results, path)
            self._set_status(f"저장됨: {path}")

    def export_scan(self):
        sel = self.tree_logs.selection()
        if not sel:
            messagebox.showinfo(__app_name__, "스캔 로그 행을 선택하세요.")
            return
        r = self.log_results[int(sel[0])]
        if not r.scan:
            messagebox.showinfo(__app_name__, f"{r.name} 에 스캔 데이터가 없습니다.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")],
                                            initialfile=f"{r.name}_scan.csv")
        if not path:
            return

        def work():
            core.write_scan_csv(r, path)
            png = os.path.splitext(path)[0] + ".png"
            core.plot_scan(r.scan, png, title=r.name)
            return path, png

        _busy(self, work,
              lambda res: self._set_status(f"저장됨: {res[0]} + {res[1]}"))

    # ------------------------------------------------------------ 탭6: 외삽
    def _tab6(self):
        f = self.tab_fit
        info = ttk.Label(f, justify="left", text=(
            "탭 5에서 파싱한 결과 중 파일명에 _n숫자 태그가 있는 잡들을 자동으로 모아 "
            "사슬길이 외삽을 수행합니다.\n"
            "모델:  ① 1/n 선형  ② Kuhn 식 y(n)=y∞+(y₁−y∞)·exp(−a(n−1)) → 유효공액길이(ECL)"))
        info.pack(fill="x", **PAD)

        top = ttk.Frame(f); top.pack(fill="x", **PAD)
        ttk.Label(top, text="물성:").pack(side="left")
        self.cmb_prop = ttk.Combobox(top, state="readonly", width=16,
                                     values=["gap_eV", "homo_eV", "lumo_eV",
                                             "s1_eV", "dipole_debye"])
        self.cmb_prop.current(0)
        self.cmb_prop.pack(side="left", padx=6)
        ttk.Button(top, text="외삽 실행", command=self.run_fit).pack(side="left", **PAD)
        ttk.Button(top, text="그래프 저장...", command=self.save_fit_plot).pack(side="left", **PAD)

        self.txt_fit = Text(f, wrap=WORD)
        self.txt_fit.pack(fill="both", expand=True, **PAD)
        self._last_fit = None

    def _series(self, prop: str):
        ns, ys = [], []
        for r in self.log_results:
            v = getattr(r, prop, None)
            if r.n is not None and v is not None and r.normal_termination:
                ns.append(float(r.n))
                ys.append(float(v))
        return ns, ys

    def run_fit(self):
        prop = self.cmb_prop.get()
        ns, ys = self._series(prop)
        if len(ns) < 3:
            messagebox.showinfo(__app_name__,
                                f"'{prop}' 값이 있는 _n 태그 잡이 3개 이상 필요합니다. "
                                f"(현재 {len(ns)}개)")
            return
        try:
            fit = core.extrapolate(ns, ys)
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror(__app_name__, str(exc))
            return
        self._last_fit = (fit, prop)
        t = self.txt_fit
        t.delete("1.0", END)
        t.insert(END, f"◆ {prop} 사슬길이 외삽\n")
        t.insert(END, f"  데이터: n = {fit['n']}\n           y = "
                      f"{[round(y, 4) for y in fit['y']]}\n\n")
        lin = fit["linear"]
        if "y_inf" in lin:
            t.insert(END, f"  [1/n 선형]  y(∞) = {lin['y_inf']:.4f}    "
                          f"기울기 = {lin['slope_a']:.4f}    R² = {lin['r2']:.5f}\n")
        k = fit["kuhn"]
        if "y_inf" in k:
            t.insert(END, f"  [Kuhn]      y(∞) = {k['y_inf']:.4f}    "
                          f"y(1) = {k['y_1']:.4f}    a = {k['decay_a']:.4f}    "
                          f"R² = {k['r2']:.5f}\n")
            t.insert(END, f"              유효공액길이(95% 포화) ≈ "
                          f"{k['ECL_units_95pct']:.1f} 반복단위\n")
        t.insert(END, "\n※ 포화 구간이 보이면 Kuhn 값을, 데이터가 짧으면 두 모델의 "
                      "차이를 불확도로 보고하세요.\n")
        self._set_status(f"{prop} 외삽 완료")

    def save_fit_plot(self):
        if self._last_fit is None:
            messagebox.showinfo(__app_name__, "먼저 외삽을 실행하세요.")
            return
        fit, prop = self._last_fit
        path = filedialog.asksaveasfilename(defaultextension=".png",
                                            filetypes=[("PNG", "*.png")],
                                            initialfile=f"extrapolation_{prop}.png")
        if not path:
            return
        _busy(self, lambda: core.plot_extrapolation(fit, prop, path),
              lambda p: self._set_status(f"저장됨: {p}"))

    # ------------------------------------------------------------ 도움말
    def _help_chemdraw(self):
        messagebox.showinfo("ChemDraw에서 구조 가져오기", (
            "권장 경로 3가지:\n\n"
            "① 파일 저장 (가장 안정적)\n"
            "   ChemDraw: File → Save As → ChemDraw XML (*.cdxml)\n"
            "   → MolScope 탭1에서 파일 열기\n\n"
            "② 네이티브 .cdx\n"
            "   그대로 열 수 있습니다 (이 빌드는 .cdx 지원).\n"
            "   읽기 실패 시 .cdxml 로 다시 저장하세요.\n\n"
            "③ 클립보드 SMILES\n"
            "   구조 선택 → Edit → Copy As → SMILES\n"
            "   → MolScope 탭1의 SMILES 칸에 붙여넣기\n\n"
            "올리고머용 반복단위는 양 끝에 attachment point 2개를 "
            "그리거나 SMILES 의 [*] 를 사용하세요."))

    def _about(self):
        messagebox.showinfo("정보", (
            f"{__app_name__} v{__version__}\n\n"
            "ChemDraw 구조 → 분석 → 3D → Gaussian 입력 → 결과 파싱 → 외삽\n"
            "엔진: RDKit / matplotlib / Pillow\n\n"
            "이 도구는 계산을 '준비/해석'합니다. DFT 자체는 Gaussian 에서, "
            "파동함수 해석(LOL-π, hole-electron 등)은 Multiwfn 에서 수행하세요."))


def main():
    root = Tk()
    try:
        from tkinter import font
        for name in ("TkDefaultFont", "TkTextFont", "TkMenuFont"):
            font.nametofont(name).configure(size=10)
    except Exception:
        pass
    style = ttk.Style()
    if "vista" in style.theme_names():
        style.theme_use("vista")
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
