"""Self-contained AlphaFold DB fetch + MD-grade structure repair.

Vendored from dd_prep (AFDB-only, MD-repair-only path) so dd_afpocket has
no inter-repo runtime dependency -- dd_prep's general PDB/hetero-
classification/docking-repair machinery (needed by dd_prep's other
consumers) is not needed here: AFDB models are always a single predicted
chain with no ligands and a B-factor column holding pLDDT, so only the
AFDB-specific trimming + PDBFixer MD-protonation path applies.
"""
from __future__ import annotations

import json
import math
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"


# ---- fetch ----

def resolve_afdb_pdb_url(uniprot_id: str) -> str:
    """Look up the current model version's PDB URL for a UniProt accession
    via the AlphaFold DB REST API (model version numbers change over time,
    e.g. v4 -> v6, so this must not be hardcoded)."""
    with urllib.request.urlopen(AFDB_API.format(uniprot_id=uniprot_id.upper())) as fh:
        entries = json.load(fh)
    if not entries:
        raise ValueError(f"AlphaFold DB has no entry for {uniprot_id!r}")
    return entries[0]["pdbUrl"]


def download_afdb(uniprot_id: str, dest: Path) -> str:
    """Fetch the current AlphaFold DB model for a UniProt accession and
    return its text contents. Caches to `dest`: if it already exists, the
    download is skipped and the cached file is used."""
    dest = Path(dest)
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        url = resolve_afdb_pdb_url(uniprot_id)
        urllib.request.urlretrieve(url, dest)
    return dest.read_text()


# ---- minimal fixed-column PDB parsing ----

def _altloc_ok(line: str) -> bool:
    a = line[16] if len(line) > 16 else " "
    return a in (" ", "A")


def chains_present(text: str) -> List[str]:
    """Distinct chain IDs among ATOM records, in first-appearance order."""
    seen: List[str] = []
    for ln in text.splitlines():
        if ln[:6] == "ATOM  " and _altloc_ok(ln):
            c = ln[21]
            if c not in seen:
                seen.append(c)
    return seen


def select_protein(text: str, chains: Optional[Sequence[str]] = None) -> List[str]:
    """ATOM (+ TER) lines for the given chains (None = all chains present),
    keeping only the primary altloc."""
    keep = set(chains) if chains is not None else None
    out: List[str] = []
    for ln in text.splitlines():
        rec = ln[:6]
        if rec == "ATOM  " and _altloc_ok(ln) and (keep is None or ln[21] in keep):
            out.append(ln)
        elif rec == "TER   " and (keep is None or ln[21:22] in keep):
            out.append(ln)
    return out


def residue_b_factors(atom_lines: Sequence[str]) -> Dict[Tuple[str, int], float]:
    """Average B-factor (columns 61-66) per (chain, resseq) -- for
    AlphaFold models this column holds per-residue pLDDT."""
    sums: Dict[Tuple[str, int], List[float]] = {}
    for ln in atom_lines:
        if ln[:6] != "ATOM  ":
            continue
        try:
            key = (ln[21], int(ln[22:26]))
            b = float(ln[60:66])
        except ValueError:
            continue
        sums.setdefault(key, []).append(b)
    return {k: sum(v) / len(v) for k, v in sums.items()}


# ---- pLDDT-based terminal trimming ----

def _ordered_residues(atom_lines: Sequence[str]) -> List[Tuple[str, int]]:
    seen: List[Tuple[str, int]] = []
    seen_set = set()
    for ln in atom_lines:
        if ln[:6] != "ATOM  ":
            continue
        key = (ln[21], int(ln[22:26]))
        if key not in seen_set:
            seen.append(key)
            seen_set.add(key)
    return seen


def trim_low_confidence_termini(atom_lines: Sequence[str], cutoff: float = 50.0):
    """Drop residues at the very start/end of each chain whose pLDDT is
    below `cutoff`, stopping at the first residue (from each end) that
    meets it -- internal low-confidence loops are left alone since removing
    them would fragment an otherwise-ordered domain on a heuristic alone.

    Returns (kept_lines, n_trimmed_n_term, n_trimmed_c_term).
    """
    plddt = residue_b_factors(atom_lines)
    by_chain: dict = {}
    for chain, resseq in _ordered_residues(atom_lines):
        by_chain.setdefault(chain, []).append((chain, resseq))

    drop = set()
    n_trim_n = 0
    n_trim_c = 0
    for chain, residues in by_chain.items():
        lo = 0
        while lo < len(residues) and plddt.get(residues[lo], 100.0) < cutoff:
            drop.add(residues[lo])
            lo += 1
            n_trim_n += 1
        hi = len(residues) - 1
        while hi >= lo and plddt.get(residues[hi], 100.0) < cutoff:
            drop.add(residues[hi])
            hi -= 1
            n_trim_c += 1

    kept = [
        ln for ln in atom_lines
        if ln[:6] != "ATOM  " or (ln[21], int(ln[22:26])) not in drop
    ]
    return kept, n_trim_n, n_trim_c


def confidence_summary(atom_lines: Sequence[str]) -> dict:
    """Fraction of residues in each AFDB confidence bin."""
    plddt = residue_b_factors(atom_lines)
    values = list(plddt.values())
    n = len(values) or 1
    return {
        "very_low": sum(1 for v in values if v < 50) / n,
        "low": sum(1 for v in values if 50 <= v < 70) / n,
        "confident": sum(1 for v in values if 70 <= v < 90) / n,
        "very_high": sum(1 for v in values if v >= 90) / n,
    }


# ---- tidy (renumber, TER insertion, disulfide CYX rename) ----

def tidy_structure(protein_lines: Sequence[str], out_pdb: Path, ss_cutoff: float = 2.5) -> int:
    """Renumber residues 1..N per chain, insert TER at chain boundaries and
    backbone breaks, and rename disulfide-bonded CYS to CYX. Returns the
    number of disulfide bonds found."""
    lines = [ln for ln in protein_lines if ln[:6] in ("ATOM  ", "HETATM")]

    new_key = {}
    cur = 0
    prev_rid = None
    prev_chain = None
    sg = []
    bbN, bbC = {}, {}
    for n, ln in enumerate(lines):
        chain = ln[21]
        rid = (chain, ln[22:26], ln[26])
        if chain != prev_chain:
            cur = 0
            prev_rid = None
        if rid != prev_rid:
            cur += 1
            prev_rid = rid
        prev_chain = chain
        key = (chain, cur)
        new_key[n] = key
        name = ln[12:16].strip()
        crd = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        if name == "N":
            bbN[key] = crd
        elif name == "C":
            bbC[key] = crd
        elif name == "SG" and ln[17:20].strip() == "CYS":
            sg.append((key, crd))

    ss_idx = set()
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            if math.dist(sg[i][1], sg[j][1]) < ss_cutoff:
                ss_idx.add(sg[i][0])
                ss_idx.add(sg[j][0])

    by_chain = {}
    for key in new_key.values():
        by_chain.setdefault(key[0], set()).add(key[1])
    breaks = set()
    for chain, idxs in by_chain.items():
        top = max(idxs)
        for idx in range(1, top):
            k1, k2 = (chain, idx), (chain, idx + 1)
            if k1 in bbC and k2 in bbN:
                if math.dist(bbC[k1], bbN[k2]) > 2.0:
                    breaks.add(k1)
            else:
                breaks.add(k1)

    out = []
    prev_key = None
    for n, ln in enumerate(lines):
        key = new_key[n]
        if prev_key is not None and key != prev_key:
            if key[0] != prev_key[0] or prev_key in breaks:
                out.append("TER")
        prev_key = key
        chain, idx = key
        rn = "CYX" if (ln[17:20] == "CYS" and key in ss_idx) else ln[17:20]
        ln = ln[:17] + rn + ln[20:22] + f"{idx:4d}" + " " + ln[27:]
        out.append(ln)
    if out:
        out.append("TER")
    Path(out_pdb).write_text("\n".join(out) + "\nEND\n")
    return len(ss_idx) // 2


# ---- PDBFixer repair (heavy deps imported lazily) ----

def _resnum(residue) -> str:
    return residue.id


def fix_structure(in_pdb: Path, out_pdb: Path) -> dict:
    """Replace nonstandard residues and complete missing heavy atoms
    (never models missing loops -- dd_afpocket's apo AFDB inputs have no
    experimental gaps to model in the first place).

    Returns a stats dict with both summary counts (`missing_residues_found`,
    `missing_atoms_added`) and per-change detail needed for the report:
      - `residue_renames`: [{chain, resnum, from, to}, ...]
      - `gaps`: [{chain, after_residue, before_residue, n_residues,
        residue_names}, ...]
      - `missing_atom_details`: [{chain, resnum, resname, atoms}, ...]
    """
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(in_pdb))

    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    gaps = []
    for (chain_index, index), res_names in sorted(fixer.missingResidues.items()):
        chain = chains[chain_index]
        residues = list(chain.residues())
        gaps.append({
            "chain": chain.id,
            "after_residue": _resnum(residues[index - 1]) if index > 0 else None,
            "before_residue": _resnum(residues[index]) if index < len(residues) else None,
            "n_residues": len(res_names),
            "residue_names": list(res_names),
        })
    n_missing_residues = sum(len(v) for v in fixer.missingResidues.values())
    fixer.missingResidues = {}

    fixer.findNonstandardResidues()
    renames = [
        {"chain": residue.chain.id, "resnum": _resnum(residue), "from": residue.name, "to": new_name}
        for residue, new_name in fixer.nonstandardResidues
    ]
    fixer.replaceNonstandardResidues()

    fixer.findMissingAtoms()
    atoms_by_residue: dict = {}
    for residue, atoms in fixer.missingAtoms.items():
        key = (residue.chain.id, _resnum(residue), residue.name)
        atoms_by_residue.setdefault(key, []).extend(atom.name for atom in atoms)
    for residue, atom_names in fixer.missingTerminals.items():
        # missingTerminals holds plain atom-name strings (e.g. "OXT"), unlike
        # missingAtoms which holds template Atom objects -- an asymmetry in
        # PDBFixer's own API, not a typo here.
        key = (residue.chain.id, _resnum(residue), residue.name)
        atoms_by_residue.setdefault(key, []).extend(atom_names)
    missing_atom_details = [
        {"chain": chain_id, "resnum": resnum, "resname": resname, "atoms": atoms}
        for (chain_id, resnum, resname), atoms in sorted(atoms_by_residue.items())
    ]
    n_missing_atoms = sum(len(atoms) for atoms in atoms_by_residue.values())
    fixer.addMissingAtoms()

    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    return {
        "missing_residues_found": n_missing_residues,
        "missing_atoms_added": n_missing_atoms,
        "residue_renames": renames,
        "gaps": gaps,
        "missing_atom_details": missing_atom_details,
    }


def add_hydrogens_md(in_pdb: Path, out_pdb: Path, *, ph: float = 7.0) -> None:
    """Protonate at the given pH via PDBFixer's built-in pKa rules (His
    tautomer choice, Asp/Glu/Lys/Cys protonation state). A quick MD-ready
    structure, not a substitute for PROPKA/H++ when precise pKa matters."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(in_pdb))
    fixer.addMissingHydrogens(ph)
    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)


# ---- report ----

@dataclass
class PrepReport:
    input_id: str
    chains_kept: List[str] = field(default_factory=list)
    n_residues: int = 0
    n_chain_breaks: int = 0
    n_disulfides: int = 0
    missing_residues_found: int = 0
    missing_atoms_added: int = 0
    residue_renames: List[Dict] = field(default_factory=list)
    gaps: List[Dict] = field(default_factory=list)
    missing_atom_details: List[Dict] = field(default_factory=list)
    plddt_trimmed_n_term: int = 0
    plddt_trimmed_c_term: int = 0
    plddt_bins: Optional[Dict[str, float]] = None
    md_pdb: Optional[str] = None

    def to_json(self, path: Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2, ensure_ascii=False))


def print_report(report: PrepReport) -> None:
    print(f"[{report.input_id}] chains={','.join(report.chains_kept)} "
          f"residues={report.n_residues} chain_breaks={report.n_chain_breaks} "
          f"disulfides={report.n_disulfides}", flush=True)
    if report.missing_residues_found or report.missing_atoms_added:
        print(f"[{report.input_id}] PDBFixer: {report.missing_residues_found} missing residue(s) found, "
              f"{report.missing_atoms_added} missing atom(s) added", flush=True)
    for r in report.residue_renames:
        print(f"[{report.input_id}] PDBFixer: renamed {r['chain']}:{r['resnum']} {r['from']} -> {r['to']}",
              flush=True)
    for g in report.gaps:
        span = f"{g['after_residue']}..{g['before_residue']}" if g["after_residue"] and g["before_residue"] \
            else (f"before {g['before_residue']}" if g["before_residue"] else f"after {g['after_residue']}")
        print(f"[{report.input_id}] PDBFixer: chain {g['chain']} gap {span} "
              f"({g['n_residues']} residue(s): {','.join(g['residue_names'])}) -- left as gap", flush=True)
    for d in report.missing_atom_details:
        print(f"[{report.input_id}] PDBFixer: added atom(s) {','.join(d['atoms'])} "
              f"to {d['chain']}:{d['resnum']} {d['resname']}", flush=True)
    if report.plddt_bins is not None:
        bins = ", ".join(f"{k}={v:.0%}" for k, v in report.plddt_bins.items())
        print(f"[{report.input_id}] pLDDT: {bins}; trimmed {report.plddt_trimmed_n_term} N-term / "
              f"{report.plddt_trimmed_c_term} C-term residue(s)", flush=True)
    if report.md_pdb:
        print(f"[{report.input_id}] -> {report.md_pdb}", flush=True)


# ---- orchestration ----

def prepare_afdb_structure(
    input_pdb: Path, out_dir: Path, *, input_id: str, ph: float = 7.0,
    plddt_cutoff: float = 50.0, ss_cutoff: float = 2.5, show_progress: bool = True,
) -> PrepReport:
    """Turn one already-downloaded raw AFDB model into an MD-grade,
    protonated structure (`<input_id>_md.pdb`)."""
    input_pdb = Path(input_pdb)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    text = input_pdb.read_text()
    chain_list = chains_present(text)
    protein_lines = select_protein(text, chain_list)
    if not protein_lines:
        raise ValueError(f"{input_pdb}: no ATOM lines found for chain(s) {chain_list}")

    report = PrepReport(input_id=input_id, chains_kept=chain_list)

    protein_lines, n_trim_n, n_trim_c = trim_low_confidence_termini(protein_lines, cutoff=plddt_cutoff)
    report.plddt_trimmed_n_term = n_trim_n
    report.plddt_trimmed_c_term = n_trim_c
    report.plddt_bins = confidence_summary(protein_lines)

    tidy_chains = {ln[21] for ln in protein_lines if ln[:6] in ("ATOM  ", "HETATM")}
    tidied_pdb = out_dir / f"{input_id}_tidy.pdb"
    report.n_disulfides = tidy_structure(protein_lines, tidied_pdb, ss_cutoff=ss_cutoff)
    report.n_residues = len({
        (ln[21], ln[22:26]) for ln in protein_lines if ln[:6] in ("ATOM  ", "HETATM")
    })
    n_ter = sum(1 for ln in tidied_pdb.read_text().splitlines() if ln.rstrip() == "TER")
    report.n_chain_breaks = max(0, n_ter - len(tidy_chains))

    fixed_pdb = out_dir / f"{input_id}_fixed.pdb"
    repair_stats = fix_structure(tidied_pdb, fixed_pdb)
    report.missing_residues_found = repair_stats["missing_residues_found"]
    report.missing_atoms_added = repair_stats["missing_atoms_added"]
    report.residue_renames = repair_stats["residue_renames"]
    report.gaps = repair_stats["gaps"]
    report.missing_atom_details = repair_stats["missing_atom_details"]

    md_pdb = out_dir / f"{input_id}_md.pdb"
    add_hydrogens_md(fixed_pdb, md_pdb, ph=ph)
    report.md_pdb = str(md_pdb)

    report.to_json(out_dir / f"{input_id}_report.json")
    if show_progress:
        print_report(report)
    return report


def fetch_and_prepare_afdb(
    uniprot_id: str, raw_dir: Path, out_dir: Path, *, ph: float = 7.0,
    plddt_cutoff: float = 50.0, show_progress: bool = True,
) -> PrepReport:
    """Download the current AlphaFold DB model for a UniProt accession,
    then prepare it into an MD-grade, protonated structure."""
    raw_path = Path(raw_dir) / f"{uniprot_id.upper()}_afdb_raw.pdb"
    download_afdb(uniprot_id, raw_path)
    return prepare_afdb_structure(
        raw_path, out_dir, input_id=uniprot_id.lower(), ph=ph,
        plddt_cutoff=plddt_cutoff, show_progress=show_progress,
    )
