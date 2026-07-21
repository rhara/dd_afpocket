"""Position restraints: pin every protein atom except those around a
detected druggable pocket to their starting coordinates, so a short MD run
samples pocket-local flexibility without paying (or risking numerical
instability from) whole-protein dynamics.

This is a self-contained generalization of `dd_md/restraints.py`'s
`freeze_distal_and_gap_residues` (dd_afpocket does not import dd_md): that
function defines "pocket" as "close to the docked ligand", which has no
analogue here, since dd_afpocket's inputs are ligand-free (apo) AlphaFold models.
Instead the pocket is defined by `pocket.py`'s fpocket-detected lining
residues and cavity centroid. The harmonic-restraint mechanics (a
`CustomExternalForce` back to each atom's starting position, not
`system.setParticleMass(index, 0)`) and default force constant
(1000 kJ/mol/nm^2) are reused as-is from dd_md, which verified empirically
that this value runs stably combined with a 4 fs HMR timestep (a stiffer
spring, 1e5, produced `Particle coordinate is NaN` within a few ps) and
gives an RMS thermal fluctuation of ~0.86 A at 300 K -- small next to the
scale of pocket-shape differences this project's clustering step looks for.
"""
from __future__ import annotations

from typing import Any, Dict, Sequence, Set, Tuple

from .pocket import Residue

_SOLVENT_RESNAMES = {"HOH", "WAT", "NA", "CL", "K", "MG", "ZN"}


def mobile_residue_indices(
    topology: Any, positions: Any, pocket_residues: Sequence[Residue], pocket_center_nm: Tuple[float, float, float], *,
    mobile_radius_nm: float = 1.0, mobile_margin_residues: int = 2,
) -> Set[int]:
    """0-based residue indices to leave mobile: the fpocket-detected
    `pocket_residues` (matched by chain id + resSeq, which -- since
    `pocket.py` and this function are always run against the very same
    prepped PDB -- line up exactly with the topology's own residue
    numbering) union every residue with a CA atom within
    `mobile_radius_nm` of `pocket_center_nm`, then expanded by
    `mobile_margin_residues` sequence neighbors on each side of every
    mobile residue (so a lone mobile residue is never sandwiched directly
    between frozen neighbors, which would strain the backbone at that
    junction -- the same idea as dd_md's gap-margin residues, applied to a
    different target).
    """
    import numpy as np
    from openmm import unit

    if hasattr(positions, "value_in_unit"):
        positions = positions.value_in_unit(unit.nanometer)
    pos = np.asarray(positions)
    center = np.asarray(pocket_center_nm)

    pocket_keys = {(r.chain, r.resnum) for r in pocket_residues}
    residues = list(topology.residues())

    mobile: Set[int] = set()
    for residue in residues:
        if residue.name in _SOLVENT_RESNAMES:
            continue
        key = (residue.chain.id, int(residue.id))
        if key in pocket_keys:
            mobile.add(residue.index)
            continue
        for atom in residue.atoms():
            if atom.name == "CA" and float(np.linalg.norm(pos[atom.index] - center)) <= mobile_radius_nm:
                mobile.add(residue.index)
                break

    if mobile_margin_residues > 0:
        expanded = set(mobile)
        for idx in mobile:
            lo, hi = max(0, idx - mobile_margin_residues), min(len(residues), idx + mobile_margin_residues + 1)
            expanded.update(range(lo, hi))
        mobile = expanded
    return mobile


def freeze_non_pocket_residues(
    system: Any, topology: Any, positions: Any, mobile_residues: Set[int], *,
    force_constant_kj_per_mol_nm2: float = 1.0e3,
) -> Dict[str, int]:
    """Add a `CustomExternalForce` to `system` (in place) restraining every
    atom of every non-mobile protein residue to its starting position.
    Residues in `mobile_residues` (see `mobile_residue_indices`) are left
    completely free, side chain and backbone alike. Solvent/ion residues
    are never restrained (there normally are none: dd_afpocket samples apo
    protein-only systems in implicit solvent). Returns residue counts for
    reporting.
    """
    import numpy as np
    from openmm import CustomExternalForce, unit

    if hasattr(positions, "value_in_unit"):
        pos_nm = np.asarray(positions.value_in_unit(unit.nanometer))
    else:
        pos_nm = np.asarray(positions)

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", force_constant_kj_per_mol_nm2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    n_frozen = n_mobile = 0
    for residue in topology.residues():
        if residue.name in _SOLVENT_RESNAMES:
            continue
        if residue.index in mobile_residues:
            n_mobile += 1
            continue
        n_frozen += 1
        for atom in residue.atoms():
            x0, y0, z0 = pos_nm[atom.index]
            force.addParticle(atom.index, [x0, y0, z0])

    if force.getNumParticles() > 0:
        system.addForce(force)
    return {"residues_frozen": n_frozen, "residues_mobile": n_mobile}
