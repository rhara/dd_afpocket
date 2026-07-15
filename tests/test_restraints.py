from _helpers import atom_line

from dd_af.pocket import Residue
from dd_af.restraints import freeze_non_pocket_residues, mobile_residue_indices

# A tiny synthetic 5-residue chain, spaced 4 A apart along x so distance-
# based mobile-set membership is easy to reason about. Residue 3 (index 2)
# sits at the pocket center.
_PDB_TEXT = "\n".join(
    atom_line(i + 1, "CA", "ALA", "A", i + 1, float(i) * 4.0, 0.0, 0.0)
    for i in range(5)
) + "\n"


def _load_topology(tmp_path):
    from openmm.app import PDBFile

    pdb_path = tmp_path / "chain.pdb"
    pdb_path.write_text(_PDB_TEXT)
    return PDBFile(str(pdb_path))


def test_mobile_residue_indices_includes_pocket_residues_and_radius_neighbors(tmp_path):
    pdb = _load_topology(tmp_path)
    # Pocket-lining residue: resnum 1 (index 0), far from the geometric
    # center used below -- included regardless of distance.
    pocket_residues = [Residue("A", 1)]
    # Center at residue 3's position (index 2, x=8.0 A = 0.8 nm), radius
    # wide enough to reach residue 3 only (residues 2/4 are 4A=0.4nm away,
    # i.e. within a 0.5nm radius too -- use a tight 0.05nm radius so only
    # the exact-center residue qualifies via distance).
    center_nm = (0.8, 0.0, 0.0)

    mobile = mobile_residue_indices(
        pdb.topology, pdb.positions, pocket_residues, center_nm,
        mobile_radius_nm=0.05, mobile_margin_residues=0,
    )
    # index 0 (explicit pocket residue) + index 2 (within radius of center)
    assert mobile == {0, 2}


def test_mobile_residue_indices_applies_sequence_margin(tmp_path):
    pdb = _load_topology(tmp_path)
    pocket_residues = [Residue("A", 3)]  # index 2
    center_nm = (1000.0, 1000.0, 1000.0)  # far away: no distance-based additions

    mobile = mobile_residue_indices(
        pdb.topology, pdb.positions, pocket_residues, center_nm,
        mobile_radius_nm=0.001, mobile_margin_residues=1,
    )
    # index 2 plus its +-1 sequence neighbors (indices 1 and 3)
    assert mobile == {1, 2, 3}


def _build_toy_system(pdb):
    from openmm import NonbondedForce, System

    system = System()
    nb = NonbondedForce()
    for _ in range(pdb.topology.getNumAtoms()):
        system.addParticle(12.0)
        nb.addParticle(0.0, 0.3, 0.0)
    system.addForce(nb)
    return system


def test_freeze_non_pocket_residues_restrains_only_non_mobile_atoms(tmp_path):
    pdb = _load_topology(tmp_path)
    system = _build_toy_system(pdb)
    mobile = {2}  # only residue index 2 (1 atom: CA) stays mobile

    counts = freeze_non_pocket_residues(system, pdb.topology, pdb.positions, mobile)

    assert counts == {"residues_frozen": 4, "residues_mobile": 1}
    # One CustomExternalForce was added, restraining exactly the 4 frozen
    # atoms (1 CA atom per residue in this toy topology).
    forces = [f for f in system.getForces() if f.__class__.__name__ == "CustomExternalForce"]
    assert len(forces) == 1
    assert forces[0].getNumParticles() == 4


def test_freeze_non_pocket_residues_noop_when_everything_is_mobile(tmp_path):
    pdb = _load_topology(tmp_path)
    system = _build_toy_system(pdb)
    mobile = {0, 1, 2, 3, 4}

    counts = freeze_non_pocket_residues(system, pdb.topology, pdb.positions, mobile)

    assert counts == {"residues_frozen": 0, "residues_mobile": 5}
    forces = [f for f in system.getForces() if f.__class__.__name__ == "CustomExternalForce"]
    assert len(forces) == 0  # no force added when there is nothing to restrain
