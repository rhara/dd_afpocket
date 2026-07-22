from dd_afpocket.pocket import PocketSelection, Residue
from dd_afpocket.visualize import build_pocket_candidates_view, write_pocket_candidates_html
from _helpers import atom_line

_RECEPTOR_PDB = "\n".join([
    atom_line(1, "N", "LYS", "A", 33, 0.0, 0.0, 0.0),
    atom_line(2, "CA", "LYS", "A", 33, 1.0, 2.0, 3.0),
    atom_line(3, "C", "LYS", "A", 33, 2.0, 0.0, 0.0),
    atom_line(4, "CA", "ASP", "A", 127, 10.0, 10.0, 10.0),
]) + "\n"


def _make_selection(tmp_path, *, rank, fpocket_id, druggability, residues):
    """A minimal PocketSelection plus its own fpocket_out_dir/pockets/pocketN_vert.pqr,
    matching what top_pocket_candidates would hand to build_pocket_candidates_view."""
    out_dir = tmp_path / f"rank{rank}_out"
    (out_dir / "pockets").mkdir(parents=True)
    vert_pqr = "\n".join([
        f"ATOM  {i + 1:>5}    C STP{fpocket_id:>6}    {x:8.3f}{y:8.3f}{z:8.3f}    0.00     3.50"
        for i, (x, y, z) in enumerate([(1.0, 1.0, 1.0), (2.0, 2.0, 2.0)])
    ]) + "\n"
    (out_dir / "pockets" / f"pocket{fpocket_id}_vert.pqr").write_text(vert_pqr)

    return PocketSelection(
        receptor_pdb="receptor.pdb", fpocket_id=fpocket_id, rank=rank, score=0.0,
        druggability_score=druggability, n_alpha_spheres=2, volume=100.0,
        residues=list(residues), center=(0.0, 0.0, 0.0), box_center=[0.0, 0.0, 0.0], box_size=[10.0, 10.0, 10.0],
        fpocket_out_dir=str(out_dir),
    )


def test_build_pocket_candidates_view_embeds_receptor_and_pocket_models(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text(_RECEPTOR_PDB)
    selections = [
        _make_selection(tmp_path, rank=1, fpocket_id=22, druggability=0.646, residues=[Residue("A", 33)]),
        _make_selection(tmp_path, rank=2, fpocket_id=8, druggability=0.023, residues=[Residue("A", 127)]),
    ]

    view = build_pocket_candidates_view(pdb_path, selections)
    html = view._make_html()
    assert "LYS" in html  # receptor model text made it into the view
    assert "STP" in html  # alpha-sphere (pqr) model text made it into the view


def test_write_pocket_candidates_html_writes_legend_and_labels(tmp_path):
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text(_RECEPTOR_PDB)
    selections = [
        _make_selection(tmp_path, rank=1, fpocket_id=22, druggability=0.646, residues=[Residue("A", 33)]),
    ]

    out_path = write_pocket_candidates_html(pdb_path, selections, tmp_path / "pocket_candidates.html")
    html = out_path.read_text()
    assert "Rank 1" in html
    assert "0.646" in html
    assert "K33" in html  # residue label text


def test_build_pocket_candidates_view_handles_missing_vert_pqr(tmp_path):
    """A selection whose fpocket_out_dir has no matching pocketN_vert.pqr
    (e.g. hand-built PocketSelection, not one from top_pocket_candidates)
    should not raise -- the sphere layer for that pocket is simply skipped."""
    pdb_path = tmp_path / "receptor.pdb"
    pdb_path.write_text(_RECEPTOR_PDB)
    selection = PocketSelection(
        receptor_pdb="receptor.pdb", fpocket_id=1, rank=1, score=0.0, druggability_score=0.5,
        n_alpha_spheres=0, volume=0.0, residues=[Residue("A", 33)], center=(0.0, 0.0, 0.0),
        box_center=[0.0, 0.0, 0.0], box_size=[10.0, 10.0, 10.0], fpocket_out_dir="",
    )

    view = build_pocket_candidates_view(pdb_path, [selection])
    assert "LYS" in view._make_html()
