"""Druggable pocket detection via fpocket, plus a docking-box convention for
downstream ensemble docking against the ligand-free (apo) structures this
project produces.

fpocket (https://github.com/Discngine/fpocket, added to the `mpro` conda env
via `conda install -c conda-forge fpocket`) is invoked as a CLI subprocess --
no Python bindings exist, and none are needed. Given `<name>.pdb`, it writes
`<name>_out/<name>_info.txt` (one text block per pocket, in fpocket's own
detection order) plus, per pocket, `<name>_out/pockets/pocketN_atm.pdb` (the
receptor atoms actually lining that pocket's cavity) and
`pocketN_vert.pqr` (the alpha-sphere centers that define the cavity volume
itself). We re-rank fpocket's own per-pocket order by Druggability Score
(fpocket's own [0, 1] estimate of how likely a pocket is to bind a
drug-like molecule) rather than trusting its default order, since that
default order is fpocket's general cavity "Score", which does not always
agree with Druggability Score (see `rank_pockets`).
"""
from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from .progress import PocketProgress

Coord = Tuple[float, float, float]

_INFO_FLOAT_FIELDS = {
    "Score": "score",
    "Druggability Score": "druggability_score",
    "Total SASA": "total_sasa",
    "Volume": "volume",
}
_INFO_INT_FIELDS = {
    "Number of Alpha Spheres": "n_alpha_spheres",
}
_POCKET_HEADER_RE = re.compile(r"^Pocket (\d+) :\s*$")


@dataclass(frozen=True)
class Residue:
    chain: str
    resnum: int

    def __str__(self) -> str:
        return f"{self.chain}:{self.resnum}"


def parse_residue_list(text: str) -> List[Residue]:
    """Parse `A:42,A:87,B:23` (the same `chain:resnum` notation as
    `dd_docking/pocket.py`'s `Residue`/`format_flexres`) into `Residue`s."""
    residues = []
    for tok in text.split(","):
        tok = tok.strip()
        if not tok:
            continue
        chain, resnum = tok.split(":")
        residues.append(Residue(chain, int(resnum)))
    return residues


def format_residue_list(residues: Sequence[Residue]) -> str:
    return ",".join(str(r) for r in residues)


def run_fpocket(pdb_path: Path, work_dir: Path, *, retries: int = 1) -> Path:
    """Run fpocket on a copy of `pdb_path` inside `work_dir` (fpocket always
    writes its `<stem>_out/` output directory next to its input, so the
    input is copied into `work_dir` first rather than littering wherever
    `pdb_path` lives). Returns the `<stem>_out` directory. Cached: if that
    directory already exists, fpocket is not re-run.

    Both paths are resolved to absolute before use: the subprocess call
    below sets `cwd=work_dir` (so fpocket's own `<stem>_out/` side effect
    lands in the right place) while also passing the `-f` filename
    argument -- if either `pdb_path`/`work_dir` were relative, the child
    process would resolve that filename argument against its *new* (post-
    chdir) working directory, not the caller's, silently looking for the
    file one level too deep and failing with "File ... does not exist".
    `retries` (default 1) retries once more on any other, genuinely
    transient subprocess failure; set to 0 to disable.
    """
    pdb_path = Path(pdb_path).resolve()
    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    local_pdb = work_dir / pdb_path.name
    if not local_pdb.exists():
        shutil.copy(pdb_path, local_pdb)

    out_dir = work_dir / f"{local_pdb.stem}_out"
    if not out_dir.exists():
        attempts = retries + 1
        last_result = None
        for attempt in range(attempts):
            result = subprocess.run(
                ["fpocket", "-f", str(local_pdb)], cwd=str(work_dir),
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                break
            last_result = result
        else:
            raise RuntimeError(
                f"fpocket failed after {attempts} attempt(s) on {local_pdb} "
                f"(exit {last_result.returncode}):\nSTDOUT:\n{last_result.stdout}\nSTDERR:\n{last_result.stderr}"
            )
    return out_dir


def parse_info_txt(out_dir: Path) -> pd.DataFrame:
    """Parse `<stem>_info.txt` into one row per pocket (`fpocket_id` = the
    number fpocket itself assigned, matching its `pocketN_*` file names --
    NOT re-ranked yet)."""
    out_dir = Path(out_dir)
    info_files = list(out_dir.glob("*_info.txt"))
    if not info_files:
        raise FileNotFoundError(f"{out_dir}: no *_info.txt found (fpocket run failed?)")

    text = info_files[0].read_text()
    blocks = re.split(r"\n(?=Pocket \d+ :)", text.strip())

    rows: List[Dict] = []
    for block in blocks:
        header = _POCKET_HEADER_RE.match(block.splitlines()[0])
        if not header:
            continue
        row: Dict = {"fpocket_id": int(header.group(1))}
        for line in block.splitlines()[1:]:
            if ":" not in line:
                continue
            label, value = line.split(":", 1)
            label = label.strip()
            # Exact match on the text before ':' -- a startswith check would
            # let e.g. "Volume score:" shadow "Volume :" since both begin
            # with "Volume ".
            if label in _INFO_INT_FIELDS:
                row[_INFO_INT_FIELDS[label]] = int(value.strip())
            elif label in _INFO_FLOAT_FIELDS:
                row[_INFO_FLOAT_FIELDS[label]] = float(value.strip())
        rows.append(row)

    if not rows:
        raise ValueError(f"{info_files[0]}: found no 'Pocket N :' blocks")
    return pd.DataFrame(rows)


def rank_pockets(df: pd.DataFrame) -> pd.DataFrame:
    """Re-sort by Druggability Score descending and assign a 1-based `rank`
    column -- the numbering `--pocket-rank` refers to (distinct from
    fpocket's own `fpocket_id`/detection order)."""
    out = df.sort_values("druggability_score", ascending=False).reset_index(drop=True)
    out.insert(0, "rank", range(1, len(out) + 1))
    return out


def _parse_atom_lines(pdb_text: str):
    """Yield (chain, resnum, resname, atom_name, coord, element) for every
    ATOM/HETATM line."""
    for ln in pdb_text.splitlines():
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        chain = ln[21]
        try:
            resnum = int(ln[22:26])
            coord = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        resname = ln[17:20].strip()
        atom_name = ln[12:16].strip()
        element = ln[76:78].strip() if len(ln) >= 78 else ""
        yield chain, resnum, resname, atom_name, coord, element


def lining_residues(out_dir: Path, fpocket_id: int) -> List[Residue]:
    """Residues appearing in `pocketN_atm.pdb` -- the receptor atoms fpocket
    itself reports as lining pocket `fpocket_id`'s cavity. No separate
    distance-cutoff judgment is needed: fpocket already tells us the
    contact residues."""
    atm_pdb = Path(out_dir) / "pockets" / f"pocket{fpocket_id}_atm.pdb"
    text = atm_pdb.read_text()
    seen = set()
    residues = []
    for chain, resnum, _resname, _atom_name, _coord, _element in _parse_atom_lines(text):
        key = (chain, resnum)
        if key not in seen:
            seen.add(key)
            residues.append(Residue(chain, resnum))
    return sorted(residues, key=lambda r: (r.chain, r.resnum))


def vert_pqr_path(out_dir: Path, fpocket_id: int) -> Path:
    """Path to `pocketN_vert.pqr` -- the alpha-sphere centers (Voronoi
    vertices) defining pocket `fpocket_id`'s cavity volume, as written by
    fpocket next to `pocketN_atm.pdb`."""
    return Path(out_dir) / "pockets" / f"pocket{fpocket_id}_vert.pqr"


def pocket_center(out_dir: Path, fpocket_id: int) -> Coord:
    """Centroid (Angstrom, receptor coordinate frame) of the alpha-sphere
    centers defining pocket `fpocket_id`'s cavity, from `pocketN_vert.pqr`
    (free-format PQR text: `ATOM serial name resName resSeq x y z q r`)."""
    vert_pqr = vert_pqr_path(out_dir, fpocket_id)
    xs, ys, zs = [], [], []
    for ln in vert_pqr.read_text().splitlines():
        if not ln.startswith("ATOM"):
            continue
        fields = ln.split()
        xs.append(float(fields[5]))
        ys.append(float(fields[6]))
        zs.append(float(fields[7]))
    if not xs:
        raise ValueError(f"{vert_pqr}: no alpha-sphere coordinates found")
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


# Standard three-letter -> one-letter amino-acid codes, for compact residue
# labels (e.g. "K33") in build_pocket_candidates_view; a residue whose
# resname isn't in this table (a non-standard/modified residue) falls back
# to its raw three-letter code instead (see residue_labels).
_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def residue_labels(receptor_pdb: Path, residues: Sequence[Residue]) -> Dict[Residue, Tuple[str, Coord]]:
    """For each of `residues`, read `receptor_pdb`'s own ATOM records once
    and return `{residue: (label, anchor_coord)}` -- `label` is the
    residue's one-letter amino-acid code plus its number (e.g. "K33"),
    falling back to the raw three-letter PDB resname for anything outside
    `_AA3TO1`; `anchor_coord` is that residue's CA atom position, or (for a
    residue with no CA, e.g. a HETATM ligand fragment) the centroid of every
    atom belonging to it. Used to place one 3D label marker per lining
    residue in `visualize.build_pocket_candidates_view`."""
    wanted = {(r.chain, r.resnum) for r in residues}
    resnames: Dict[Tuple[str, int], str] = {}
    ca_coord: Dict[Tuple[str, int], Coord] = {}
    sums: Dict[Tuple[str, int], List[float]] = {}
    counts: Dict[Tuple[str, int], int] = {}

    for chain, resnum, resname, atom_name, coord, _element in _parse_atom_lines(Path(receptor_pdb).read_text()):
        key = (chain, resnum)
        if key not in wanted:
            continue
        resnames.setdefault(key, resname)
        if atom_name == "CA" and key not in ca_coord:
            ca_coord[key] = coord
        s = sums.setdefault(key, [0.0, 0.0, 0.0])
        s[0] += coord[0]
        s[1] += coord[1]
        s[2] += coord[2]
        counts[key] = counts.get(key, 0) + 1

    out: Dict[Residue, Tuple[str, Coord]] = {}
    for r in residues:
        key = (r.chain, r.resnum)
        if key not in resnames:
            continue
        label = f"{_AA3TO1.get(resnames[key], resnames[key])}{r.resnum}"
        if key in ca_coord:
            coord = ca_coord[key]
        else:
            n = counts[key]
            s = sums[key]
            coord = (s[0] / n, s[1] / n, s[2] / n)
        out[r] = (label, coord)
    return out


def compute_box(receptor_pdb: Path, residues: Sequence[Residue], padding: float = 5.0):
    """Docking box center/size (each axis) spanning the heavy-atom
    coordinates of `residues` in `receptor_pdb`, plus `padding` Angstrom on
    every side -- the same "bounding box of a point cloud plus padding"
    convention as `compute_box` in `dd_docking/pocket.py`, applied to
    pocket-lining-residue atoms instead of ligand atoms (there is no ligand:
    dd_afpocket's output structures are apo)."""
    wanted = {(r.chain, r.resnum) for r in residues}
    xs, ys, zs = [], [], []
    for chain, resnum, _resname, _atom_name, coord, element in _parse_atom_lines(Path(receptor_pdb).read_text()):
        if (chain, resnum) in wanted and element != "H":
            xs.append(coord[0])
            ys.append(coord[1])
            zs.append(coord[2])
    if not xs:
        raise ValueError("compute_box: no heavy atoms found for the given residues")
    center = [round((min(v) + max(v)) / 2, 3) for v in (xs, ys, zs)]
    size = [round((max(v) - min(v)) + 2 * padding, 3) for v in (xs, ys, zs)]
    return center, size


@dataclass
class PocketSelection:
    receptor_pdb: str
    fpocket_id: int
    rank: int
    score: float
    druggability_score: float
    n_alpha_spheres: int
    volume: float
    residues: List[Residue]
    center: Coord
    box_center: List[float]
    box_size: List[float]
    # Not part of to_report_dict()'s on-disk JSON schema -- an internal
    # handle back to fpocket's own <stem>_out/ directory (where
    # vert_pqr_path/lining_residues/pocket_center already read from) so
    # visualize.build_pocket_candidates_view can find this pocket's
    # pocketN_vert.pqr without the caller having to re-derive the path.
    fpocket_out_dir: str = ""

    def to_report_dict(self) -> Dict:
        return {
            "receptor_pdb": self.receptor_pdb,
            "fpocket_id": self.fpocket_id,
            "rank": self.rank,
            "score": self.score,
            "druggability_score": self.druggability_score,
            "n_alpha_spheres": self.n_alpha_spheres,
            "volume": self.volume,
            "residues": [str(r) for r in self.residues],
            "center": list(self.center),
        }

    def to_box_dict(self) -> Dict:
        return {"center": self.box_center, "size": self.box_size}


def _selection_from_row(
    receptor_pdb: Path, out_dir: Path, row, *,
    pocket_residues: Optional[Sequence[Residue]] = None, box_padding: float = 5.0,
) -> PocketSelection:
    """Build a `PocketSelection` for one already-ranked fpocket row (a
    `rank_pockets` DataFrame row). Shared by `find_druggable_pocket` (a
    single selected pocket) and `top_pocket_candidates` (the top N) so both
    apply the same residues/center/box-derivation logic."""
    fpocket_id = int(row["fpocket_id"])
    residues = list(pocket_residues) if pocket_residues else lining_residues(out_dir, fpocket_id)
    center = pocket_center(out_dir, fpocket_id)
    box_center, box_size = compute_box(receptor_pdb, residues, padding=box_padding)

    return PocketSelection(
        receptor_pdb=str(receptor_pdb), fpocket_id=fpocket_id, rank=int(row["rank"]),
        score=float(row["score"]), druggability_score=float(row["druggability_score"]),
        n_alpha_spheres=int(row["n_alpha_spheres"]), volume=float(row["volume"]),
        residues=residues, center=center, box_center=box_center, box_size=box_size,
        fpocket_out_dir=str(out_dir),
    )


def find_druggable_pocket(
    receptor_pdb: Path, work_dir: Path, *,
    pocket_rank: int = 1, pocket_residues: Optional[Sequence[Residue]] = None,
    box_padding: float = 5.0, show_progress: bool = True,
) -> PocketSelection:
    """Run fpocket on `receptor_pdb`, print the detected pocket list ranked
    by druggability, and select one pocket.

    If `pocket_residues` is given, fpocket's own residue detection is
    bypassed entirely for the *selected* pocket's residue set (fpocket is
    still run, so its ranking/center/alpha-sphere-count for informational
    pockets remain available) -- for pointing at a known (possibly
    non-top-ranked) site such as an allosteric pocket. Otherwise the
    `pocket_rank`-th pocket (1-indexed, by Druggability Score descending)
    is selected.
    """
    receptor_pdb = Path(receptor_pdb)
    out_dir = run_fpocket(receptor_pdb, work_dir)
    ranked = rank_pockets(parse_info_txt(out_dir))

    progress = PocketProgress(enabled=show_progress)
    for _, row in ranked.iterrows():
        progress.update(int(row["rank"]), row["score"], row["druggability_score"], int(row["n_alpha_spheres"]))

    if pocket_rank < 1 or pocket_rank > len(ranked):
        raise ValueError(f"--pocket-rank {pocket_rank} out of range (1..{len(ranked)} pockets found)")
    chosen = ranked.iloc[pocket_rank - 1]
    return _selection_from_row(receptor_pdb, out_dir, chosen, pocket_residues=pocket_residues, box_padding=box_padding)


def top_pocket_candidates(
    receptor_pdb: Path, work_dir: Path, *, top_n: int = 3, box_padding: float = 5.0, show_progress: bool = True,
) -> List[PocketSelection]:
    """Run fpocket (reusing `work_dir`'s cached output if `find_druggable_pocket`
    already ran there) and return a `PocketSelection` for each of the top
    `top_n` pockets by Druggability Score -- for side-by-side comparison of
    the strongest candidates (see `visualize.build_pocket_candidates_view`),
    distinct from `find_druggable_pocket`'s single-pocket selection. Returns
    fewer than `top_n` entries if fpocket found fewer pockets overall."""
    receptor_pdb = Path(receptor_pdb)
    out_dir = run_fpocket(receptor_pdb, work_dir)
    ranked = rank_pockets(parse_info_txt(out_dir))

    progress = PocketProgress(enabled=show_progress)
    for _, row in ranked.iterrows():
        progress.update(int(row["rank"]), row["score"], row["druggability_score"], int(row["n_alpha_spheres"]))

    n = min(top_n, len(ranked))
    return [_selection_from_row(receptor_pdb, out_dir, ranked.iloc[i], box_padding=box_padding) for i in range(n)]
