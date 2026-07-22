import shutil

import pytest

from dd_afpocket.pocket import (
    Residue,
    compute_box,
    find_druggable_pocket,
    format_residue_list,
    lining_residues,
    parse_info_txt,
    parse_residue_list,
    pocket_center,
    rank_pockets,
    residue_labels,
    run_fpocket,
    top_pocket_candidates,
    vert_pqr_path,
)
from _helpers import atom_line

_INFO_TXT = """Pocket 1 :
\tScore : \t0.056
\tDruggability Score : \t0.497
\tNumber of Alpha Spheres : \t54
\tTotal SASA : \t106.228
\tPolar SASA : \t10.825
\tApolar SASA : \t95.403
\tVolume : \t340.469
\tVolume score: \t 4.533

Pocket 2 :
\tScore : \t0.158
\tDruggability Score : \t0.829
\tNumber of Alpha Spheres : \t70
\tTotal SASA : \t130.438
\tPolar SASA : \t25.870
\tApolar SASA : \t104.569
\tVolume : \t570.388
\tVolume score: \t 4.471
"""

_ATM_PDB = "\n".join([
    "HEADER",
    atom_line(1, "CA", "GLU", "A", 141, -52.5, -7.0, 27.2),
    atom_line(2, "CA", "LEU", "A", 145, -53.0, -8.0, 28.0),
    atom_line(3, "CB", "LEU", "A", 145, -53.5, -8.5, 28.5),
]) + "\n"

_VERT_PQR = "\n".join([
    "HEADER",
    "ATOM      1    C STP     1     -10.000  -20.000   30.000    0.00     3.64",
    "ATOM      2    C STP     1      -8.000  -18.000   32.000    0.00     3.41",
]) + "\n"


def test_parse_info_txt_does_not_confuse_volume_and_volume_score(tmp_path):
    out_dir = tmp_path / "x_out"
    out_dir.mkdir()
    (out_dir / "x_info.txt").write_text(_INFO_TXT)

    df = parse_info_txt(out_dir)
    assert len(df) == 2
    row1 = df[df["fpocket_id"] == 1].iloc[0]
    assert row1["volume"] == pytest.approx(340.469)
    assert row1["druggability_score"] == pytest.approx(0.497)


def test_rank_pockets_sorts_by_druggability_descending(tmp_path):
    out_dir = tmp_path / "x_out"
    out_dir.mkdir()
    (out_dir / "x_info.txt").write_text(_INFO_TXT)

    ranked = rank_pockets(parse_info_txt(out_dir))
    assert ranked.iloc[0]["fpocket_id"] == 2  # higher druggability (0.829) despite being "Pocket 2"
    assert ranked.iloc[0]["rank"] == 1
    assert ranked.iloc[1]["fpocket_id"] == 1


def test_lining_residues_deduplicates_and_sorts(tmp_path):
    out_dir = tmp_path / "x_out"
    (out_dir / "pockets").mkdir(parents=True)
    (out_dir / "pockets" / "pocket1_atm.pdb").write_text(_ATM_PDB)

    residues = lining_residues(out_dir, 1)
    assert residues == [Residue("A", 141), Residue("A", 145)]


def test_pocket_center_is_alpha_sphere_centroid(tmp_path):
    out_dir = tmp_path / "x_out"
    (out_dir / "pockets").mkdir(parents=True)
    (out_dir / "pockets" / "pocket1_vert.pqr").write_text(_VERT_PQR)

    cx, cy, cz = pocket_center(out_dir, 1)
    assert cx == pytest.approx(-9.0)
    assert cy == pytest.approx(-19.0)
    assert cz == pytest.approx(31.0)


def test_compute_box_spans_residue_heavy_atoms_with_padding(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text("\n".join([
        atom_line(1, "CA", "GLU", "A", 141, 0.0, 0.0, 0.0),
        atom_line(2, "CA", "LEU", "A", 145, 2.0, 0.0, 0.0),
        atom_line(3, "H", "LEU", "A", 145, 99.0, 99.0, 99.0, element="H"),  # excluded (hydrogen)
    ]) + "\n")

    residues = [Residue("A", 141), Residue("A", 145)]
    center, size = compute_box(pdb_path, residues, padding=2.0)
    assert center == [1.0, 0.0, 0.0]
    assert size == [2.0 + 4.0, 4.0, 4.0]


def test_residue_list_roundtrip():
    residues = [Residue("A", 10), Residue("A", 20), Residue("B", 5)]
    text = format_residue_list(residues)
    assert text == "A:10,A:20,B:5"
    assert parse_residue_list(text) == residues


def test_vert_pqr_path_points_at_pockets_subdir(tmp_path):
    out_dir = tmp_path / "x_out"
    assert vert_pqr_path(out_dir, 7) == out_dir / "pockets" / "pocket7_vert.pqr"


def test_residue_labels_uses_ca_position_and_one_letter_code(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text("\n".join([
        atom_line(1, "N", "LYS", "A", 33, 0.0, 0.0, 0.0),
        atom_line(2, "CA", "LYS", "A", 33, 1.0, 2.0, 3.0),
        atom_line(3, "C", "LYS", "A", 33, 2.0, 0.0, 0.0),
    ]) + "\n")

    labels = residue_labels(pdb_path, [Residue("A", 33)])
    label, coord = labels[Residue("A", 33)]
    assert label == "K33"
    assert coord == pytest.approx((1.0, 2.0, 3.0))


def test_residue_labels_falls_back_to_atom_centroid_without_ca(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text("\n".join([
        atom_line(1, "O1", "HOH", "A", 500, 0.0, 0.0, 0.0, element="O"),
        atom_line(2, "O2", "HOH", "A", 500, 2.0, 0.0, 0.0, element="O"),
    ]) + "\n")

    labels = residue_labels(pdb_path, [Residue("A", 500)])
    label, coord = labels[Residue("A", 500)]
    assert label == "HOH500"  # not a standard amino acid -- raw resname kept
    assert coord == pytest.approx((1.0, 0.0, 0.0))


def test_residue_labels_omits_residues_absent_from_receptor(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text(atom_line(1, "CA", "GLY", "A", 1, 0.0, 0.0, 0.0) + "\n")

    labels = residue_labels(pdb_path, [Residue("A", 999)])
    assert labels == {}


def _pocket_info_block(fpocket_id, druggability, n_alpha_spheres=10, volume=100.0, score=0.0):
    return (
        f"Pocket {fpocket_id} :\n"
        f"\tScore : \t{score}\n"
        f"\tDruggability Score : \t{druggability}\n"
        f"\tNumber of Alpha Spheres : \t{n_alpha_spheres}\n"
        f"\tVolume : \t{volume}\n"
    )


def _write_cached_fpocket_output(work_dir, receptor_stem, pockets):
    """Pre-populate `work_dir/<receptor_stem>_out/` the way a real
    `run_fpocket` call would, so `run_fpocket`'s own caching (`if not
    out_dir.exists(): ...`) skips invoking the real fpocket binary --
    `pockets` is `{fpocket_id: (druggability, residues, [alpha_sphere_xyz, ...])}`."""
    out_dir = work_dir / f"{receptor_stem}_out"
    (out_dir / "pockets").mkdir(parents=True)
    info_text = "\n".join(_pocket_info_block(fid, drug) for fid, (drug, _residues, _spheres) in pockets.items())
    (out_dir / f"{receptor_stem}_info.txt").write_text(info_text)

    for fid, (_drug, residues, spheres) in pockets.items():
        atm_lines = [atom_line(i + 1, "CA", "ALA", r.chain, r.resnum, 0.0, 0.0, 0.0) for i, r in enumerate(residues)]
        (out_dir / "pockets" / f"pocket{fid}_atm.pdb").write_text("\n".join(atm_lines) + "\n")
        vert_lines = [
            f"ATOM  {i + 1:>5}    C STP{fid:>6}    {x:8.3f}{y:8.3f}{z:8.3f}    0.00     3.50"
            for i, (x, y, z) in enumerate(spheres)
        ]
        (out_dir / "pockets" / f"pocket{fid}_vert.pqr").write_text("\n".join(vert_lines) + "\n")
    return out_dir


def test_top_pocket_candidates_ranks_by_druggability_and_caps_at_top_n(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text("\n".join([
        atom_line(1, "CA", "GLU", "A", 10, 0.0, 0.0, 0.0),
        atom_line(2, "CA", "LYS", "A", 20, 5.0, 0.0, 0.0),
        atom_line(3, "CA", "ASP", "A", 30, 0.0, 5.0, 0.0),
    ]) + "\n")
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    pockets = {
        1: (0.500, [Residue("A", 10)], [(1.0, 1.0, 1.0)]),
        2: (0.900, [Residue("A", 20)], [(2.0, 2.0, 2.0)]),  # most druggable -> rank 1
        3: (0.100, [Residue("A", 30)], [(3.0, 3.0, 3.0)]),
    }
    out_dir = _write_cached_fpocket_output(work_dir, "receptor", pockets)

    top2 = top_pocket_candidates(pdb_path, work_dir, top_n=2, show_progress=False)
    assert [sel.rank for sel in top2] == [1, 2]
    assert [sel.fpocket_id for sel in top2] == [2, 1]
    assert [sel.druggability_score for sel in top2] == pytest.approx([0.900, 0.500])
    assert all(sel.fpocket_out_dir == str(out_dir) for sel in top2)

    # top_n larger than the number of pockets found -> capped, not padded/errored
    all3 = top_pocket_candidates(pdb_path, work_dir, top_n=10, show_progress=False)
    assert len(all3) == 3


def test_find_druggable_pocket_selects_pocket_rank_by_druggability(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text("\n".join([
        atom_line(1, "CA", "GLU", "A", 10, 0.0, 0.0, 0.0),
        atom_line(2, "CA", "LYS", "A", 20, 5.0, 0.0, 0.0),
    ]) + "\n")
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    pockets = {
        1: (0.200, [Residue("A", 10)], [(1.0, 1.0, 1.0)]),
        2: (0.800, [Residue("A", 20)], [(2.0, 2.0, 2.0)]),
    }
    _write_cached_fpocket_output(work_dir, "receptor", pockets)

    selection = find_druggable_pocket(pdb_path, work_dir, pocket_rank=1, show_progress=False)
    assert selection.fpocket_id == 2
    assert selection.residues == [Residue("A", 20)]
    assert selection.druggability_score == pytest.approx(0.800)


def _synthetic_helix_pdb() -> str:
    """A 20-residue idealized poly-ALA alpha helix (3.6 res/turn, 1.5 A
    rise) -- enough real 3D density for fpocket's internal clustering step
    to run without error (a single collinear/coplanar residue makes Qhull
    fail and produces no output directory at all; this is not a docking-
    realistic structure and is not expected to contain an actual druggable
    pocket, just to exercise the subprocess call and output layout)."""
    import math

    lines = []
    serial = 1
    for i in range(20):
        theta = i * 100.0 * math.pi / 180.0
        z = i * 1.5
        x, y = 2.3 * math.cos(theta), 2.3 * math.sin(theta)
        for name, dx, dy, dz in [("N", 0.0, 0.0, 0.0), ("CA", 0.5, 0.3, 0.2),
                                  ("C", 0.9, 0.6, 0.4), ("O", 1.3, 0.9, 0.6), ("CB", -0.4, 0.6, -0.3)]:
            lines.append(atom_line(serial, name, "ALA", "A", i + 1, x + dx, y + dy, z + dz))
            serial += 1
    return "\n".join(lines) + "\n"


@pytest.mark.skipif(shutil.which("fpocket") is None, reason="fpocket binary not installed")
def test_run_fpocket_produces_parseable_output_layout(tmp_path):
    pdb_path = tmp_path / "helix.pdb"
    pdb_path.write_text(_synthetic_helix_pdb())

    out_dir = run_fpocket(pdb_path, tmp_path / "work")
    info_files = list(out_dir.glob("*_info.txt"))
    assert info_files, "fpocket did not write an *_info.txt"
    # A bare synthetic helix may legitimately contain zero real pockets;
    # parse_info_txt should tolerate that empty-block case rather than
    # crashing on it.
    if info_files[0].read_text().strip():
        df = parse_info_txt(out_dir)
        assert "druggability_score" in df.columns
