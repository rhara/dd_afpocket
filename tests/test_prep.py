"""Pure-Python parts of prep.py (vendored, AFDB-only subset of dd_prep):
parsing, pLDDT trimming, tidy, and the report. Network fetch (`download_afdb`)
and PDBFixer-backed repair (`fix_structure`/`add_hydrogens_md`) are heavy/
network-dependent, same as upstream dd_prep, and are not unit tested here.
"""
from _helpers import atom_line, ter_line

from dd_afpocket.prep import (
    PrepReport,
    chains_present,
    confidence_summary,
    print_report,
    residue_b_factors,
    select_protein,
    tidy_structure,
    trim_low_confidence_termini,
)


# ---- parse ----

def test_chains_present_lists_distinct_atom_chains_in_order():
    text = "\n".join([
        atom_line(1, "CA", "ALA", "B", 1, 0, 0, 0),
        atom_line(2, "CA", "GLY", "A", 1, 0, 0, 0),
        atom_line(3, "CA", "GLY", "B", 2, 0, 0, 0),
    ])
    assert chains_present(text) == ["B", "A"]


def test_select_protein_filters_by_chain_and_altloc():
    text = "\n".join([
        atom_line(1, "CA", "ALA", "A", 1, 0, 0, 0, altloc="A"),
        atom_line(2, "CA", "ALA", "A", 1, 0.1, 0, 0, altloc="B"),
        atom_line(3, "CA", "ALA", "B", 1, 0, 0, 0),
    ])
    lines = select_protein(text, ["A"])
    assert len(lines) == 1
    assert lines[0][16] == "A"


def test_select_protein_keeps_ter_for_requested_chain():
    text = "\n".join([
        atom_line(1, "CA", "ALA", "A", 1, 0, 0, 0),
        ter_line(2, "ALA", "A", 1),
        atom_line(3, "CA", "GLY", "B", 1, 0, 0, 0),
        ter_line(4, "GLY", "B", 1),
    ])
    lines = select_protein(text, ["A"])
    assert sum(1 for ln in lines if ln.startswith("TER")) == 1


def test_residue_b_factors_averages_per_residue():
    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0, 0, 0, bfactor=40.0),
        atom_line(2, "CA", "ALA", "A", 1, 1, 0, 0, bfactor=60.0),
    ]
    b = residue_b_factors(lines)
    assert b[("A", 1)] == 50.0


# ---- pLDDT trimming ----

def _residue(chain, resnum, plddt, serial_start=1):
    return [
        atom_line(serial_start, "N", "GLY", chain, resnum, 0, 0, 0, bfactor=plddt),
        atom_line(serial_start + 1, "CA", "GLY", chain, resnum, 1, 0, 0, bfactor=plddt),
    ]


def test_trim_low_confidence_termini_removes_n_and_c_term_tails():
    lines = (
        _residue("A", 1, 20.0, 1) + _residue("A", 2, 30.0, 3)
        + _residue("A", 3, 90.0, 5) + _residue("A", 4, 95.0, 7)
        + _residue("A", 5, 25.0, 9)
    )
    kept, n_term, c_term = trim_low_confidence_termini(lines, cutoff=50.0)
    kept_residues = sorted({(ln[21], int(ln[22:26])) for ln in kept if ln[:6] == "ATOM  "})
    assert kept_residues == [("A", 3), ("A", 4)]
    assert n_term == 2
    assert c_term == 1


def test_trim_low_confidence_termini_leaves_internal_loop_alone():
    lines = (
        _residue("A", 1, 90.0, 1) + _residue("A", 2, 10.0, 3) + _residue("A", 3, 90.0, 5)
    )
    kept, n_term, c_term = trim_low_confidence_termini(lines, cutoff=50.0)
    kept_residues = sorted({(ln[21], int(ln[22:26])) for ln in kept if ln[:6] == "ATOM  "})
    assert kept_residues == [("A", 1), ("A", 2), ("A", 3)]
    assert (n_term, c_term) == (0, 0)


def test_trim_low_confidence_termini_no_trim_when_all_confident():
    lines = _residue("A", 1, 90.0) + _residue("A", 2, 95.0)
    kept, n_term, c_term = trim_low_confidence_termini(lines, cutoff=50.0)
    assert len(kept) == len(lines)
    assert (n_term, c_term) == (0, 0)


def test_confidence_summary_bins_sum_to_one():
    lines = _residue("A", 1, 20.0, 1) + _residue("A", 2, 60.0, 3) + _residue("A", 3, 80.0, 5) + _residue("A", 4, 95.0, 7)
    bins = confidence_summary(lines)
    assert abs(sum(bins.values()) - 1.0) < 1e-9
    assert bins["very_low"] == 0.25
    assert bins["very_high"] == 0.25


# ---- tidy ----

def test_tidy_structure_renumbers_residues_sequentially(tmp_path):
    lines = [
        atom_line(1, "CA", "ALA", "A", 5, 0.0, 0.0, 0.0),
        atom_line(2, "CA", "GLY", "A", 9, 1.0, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    tidy_structure(lines, out_pdb)
    resnums = [int(ln[22:26]) for ln in out_pdb.read_text().splitlines() if ln[:6] in ("ATOM  ", "HETATM")]
    assert resnums == [1, 2]


def test_tidy_structure_inserts_ter_at_backbone_break(tmp_path):
    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
        atom_line(2, "C", "ALA", "A", 1, 1.0, 0.0, 0.0),
        atom_line(3, "N", "GLY", "A", 2, 50.0, 0.0, 0.0),
        atom_line(4, "C", "GLY", "A", 2, 51.0, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    tidy_structure(lines, out_pdb)
    out_lines = out_pdb.read_text().splitlines()
    assert out_lines.count("TER") == 2


def test_tidy_structure_renames_disulfide_cys_to_cyx(tmp_path):
    lines = [
        atom_line(1, "SG", "CYS", "A", 1, 0.0, 0.0, 0.0),
        atom_line(2, "SG", "CYS", "A", 2, 2.0, 0.0, 0.0),
        atom_line(3, "SG", "CYS", "A", 3, 100.0, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    n_ss = tidy_structure(lines, out_pdb)
    assert n_ss == 1
    resnames = [ln[17:20] for ln in out_pdb.read_text().splitlines() if ln[:6] in ("ATOM  ", "HETATM")]
    assert resnames == ["CYX", "CYX", "CYS"]


def test_tidy_structure_separates_chains_and_restarts_numbering(tmp_path):
    lines = [
        atom_line(1, "N", "ALA", "A", 5, 0.0, 0.0, 0.0),
        atom_line(2, "C", "ALA", "A", 5, 1.3, 0.0, 0.0),
        atom_line(3, "N", "GLY", "A", 6, 2.6, 0.0, 0.0),
        atom_line(4, "C", "GLY", "A", 6, 3.9, 0.0, 0.0),
        atom_line(5, "CA", "VAL", "B", 100, 0.0, 0.0, 0.0),
    ]
    out_pdb = tmp_path / "tidy.pdb"
    tidy_structure(lines, out_pdb)
    out_lines = [ln for ln in out_pdb.read_text().splitlines() if ln.strip()]
    kept = [ln for ln in out_lines if ln[:6] in ("ATOM  ", "HETATM")]
    resnums_by_chain: dict = {}
    for ln in kept:
        resnums_by_chain.setdefault(ln[21], set()).add(int(ln[22:26]))
    assert resnums_by_chain == {"A": {1, 2}, "B": {1}}
    assert out_lines.count("TER") == 2


# ---- report ----

def _base_report(**overrides):
    kwargs = dict(input_id="test", chains_kept=["A"], n_residues=10)
    kwargs.update(overrides)
    return PrepReport(**kwargs)


def test_print_report_shows_residue_rename_detail(capsys):
    report = _base_report(residue_renames=[
        {"chain": "A", "resnum": "64", "from": "MSE", "to": "MET"},
    ])
    print_report(report)
    out = capsys.readouterr().out
    assert "PDBFixer: renamed A:64 MSE -> MET" in out


def test_print_report_shows_gap_detail(capsys):
    report = _base_report(gaps=[
        {"chain": "A", "after_residue": "32", "before_residue": "35",
         "n_residues": 2, "residue_names": ["LEU", "GLY"]},
    ])
    print_report(report)
    out = capsys.readouterr().out
    assert "chain A gap 32..35 (2 residue(s): LEU,GLY) -- left as gap" in out


def test_print_report_shows_missing_atom_detail(capsys):
    report = _base_report(missing_atom_details=[
        {"chain": "A", "resnum": "145", "resname": "GLY", "atoms": ["CA", "O"]},
    ])
    print_report(report)
    out = capsys.readouterr().out
    assert "PDBFixer: added atom(s) CA,O to A:145 GLY" in out


def test_print_report_no_pdbfixer_detail_lines_when_empty(capsys):
    report = _base_report()
    print_report(report)
    out = capsys.readouterr().out
    assert "PDBFixer" not in out
