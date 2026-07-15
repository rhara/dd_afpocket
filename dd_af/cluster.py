"""Structural clustering: pool frames from every sampling replica, cluster
them by pocket-atom RMSD, and pick one representative (medoid) structure per
cluster -- the "N representative pocket conformations" this project exists
to produce.

Clustering is restricted to the mobile pocket-neighborhood atoms (the same
residue set `restraints.py` left free during sampling): the rest of the
structure was position-restrained during MD and barely moved, so folding it
into a whole-protein RMSD would both dilute the pocket-shape signal that
actually matters and waste compute on a distance matrix over atoms that
carry almost no information.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Sequence

from .pocket import Residue
from .progress import ClusterProgress


@dataclass
class ClusterReportRow:
    cluster_id: int
    n_frames: int
    source_replica: int
    source_time_ps: float
    mean_intra_cluster_rmsd_a: float
    pocket_volume_nm3: float
    out_path: str


def _pocket_atom_indices(topology: Any, pocket_residues: Sequence[Residue], atom_selection: str = "ca") -> List[int]:
    """0-based atom indices for `pocket_residues` (chain id + resSeq match,
    same convention as `restraints.mobile_residue_indices`), restricted to
    CA atoms (`atom_selection="ca"`, default) or all heavy atoms
    (`"heavy"`).

    `topology` here is an **mdtraj** `Topology` (from a loaded
    trajectory), not an OpenMM one -- its `residues`/`atoms` are plain
    generator properties (no `()`), and the chain letter lives in
    `chain.chain_id`, not `chain.id`.
    """
    wanted = {(r.chain, r.resnum) for r in pocket_residues}
    indices = []
    for residue in topology.residues:
        if (residue.chain.chain_id, residue.resSeq) not in wanted:
            continue
        for atom in residue.atoms:
            if atom_selection == "ca":
                if atom.name == "CA":
                    indices.append(atom.index)
            elif atom.element is not None and atom.element.symbol != "H":
                indices.append(atom.index)
    return indices


def load_pooled_trajectory(top_pdb: Path, dcd_paths: Sequence[Path]):
    """Load and concatenate every replica's DCD (same topology) with
    mdtraj. Returns (traj, frame_replica, frame_time_ps) where the latter
    two arrays record, per pooled frame, which replica it came from and its
    time (ps) within that replica's own trajectory -- for `cluster_report.
    csv`'s provenance columns."""
    import mdtraj as md
    import numpy as np

    trajs, frame_replica, frame_time_ps = [], [], []
    for replica, dcd in enumerate(dcd_paths, start=1):
        t = md.load(str(dcd), top=str(top_pdb))
        trajs.append(t)
        frame_replica.extend([replica] * t.n_frames)
        frame_time_ps.extend((t.time if t.time is not None else np.arange(t.n_frames)).tolist())
    pooled = trajs[0] if len(trajs) == 1 else md.join(trajs)
    return pooled, np.array(frame_replica), np.array(frame_time_ps)


def pocket_volume_proxy(traj, atom_indices: Sequence[int]):
    """Per-frame convex-hull volume (nm^3) of `atom_indices`' coordinates --
    a cheap, no-rerun proxy for pocket "openness". `restraints.py` only ever
    leaves the pocket-lining residues (plus their sequence margin) mobile,
    so a hull that has grown relative to the pre-MD reference structure
    means those residues have spread apart, i.e. the pocket has opened up.

    This is *not* fpocket's own cavity volume -- reproducing that would mean
    re-running fpocket on every pooled frame, far too slow for the hundreds
    of frames a single sampling run produces (see README "Performance").
    It's a purely geometric proxy computed directly from coordinates already
    in memory, over the same pocket-lining atom selection `--cluster-atoms`
    picks for RMSD clustering.

    `qhull_options="QJ"` (joggle input) avoids `scipy.spatial.ConvexHull`
    raising `QhullError` on the exact coplanarity/near-degeneracy a handful
    of residues' CA atoms can produce for a small, flat pocket.
    """
    import numpy as np
    from scipy.spatial import ConvexHull

    if len(atom_indices) < 4:
        raise ValueError(f"pocket_volume_proxy: need >=4 pocket atoms for a convex hull, got {len(atom_indices)}")

    sub = traj.atom_slice(atom_indices)
    volumes = np.empty(sub.n_frames, dtype=float)
    for i in range(sub.n_frames):
        volumes[i] = ConvexHull(sub.xyz[i], qhull_options="QJ").volume
    return volumes


def pairwise_rmsd_matrix(traj, atom_indices: Sequence[int]):
    """All-pairs RMSD matrix (Angstrom) over `atom_indices`, each frame
    superposed onto every other frame. O(n_frames^2); fine for the
    hundreds-to-low-thousands of pooled frames this project expects."""
    import mdtraj as md
    import numpy as np

    sub = traj.atom_slice(atom_indices)
    n = sub.n_frames
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        mat[i, :] = md.rmsd(sub, sub, frame=i) * 10.0  # nm -> Angstrom
    return mat


def cluster_frames(rmsd_matrix, n_clusters: int) -> "Any":
    """`AgglomerativeClustering` (precomputed RMSD distances, average
    linkage) into `n_clusters` groups. `n_clusters` is capped at the
    number of frames available, with a warning, if it exceeds that."""
    from sklearn.cluster import AgglomerativeClustering

    n_frames = rmsd_matrix.shape[0]
    if n_clusters > n_frames:
        print(f"[cluster] --n-clusters {n_clusters} exceeds {n_frames} pooled frame(s); capping at {n_frames}", flush=True)
        n_clusters = n_frames
    model = AgglomerativeClustering(n_clusters=n_clusters, metric="precomputed", linkage="average")
    return model.fit_predict(rmsd_matrix)


def pick_medoids(rmsd_matrix, labels) -> Dict[int, int]:
    """For each cluster label, the frame index with minimum mean RMSD to
    every other member of its own cluster (never the coordinate average,
    which would fabricate an unphysical structure)."""
    import numpy as np

    medoids: Dict[int, int] = {}
    for label in sorted(set(labels.tolist())):
        members = np.where(labels == label)[0]
        sub = rmsd_matrix[np.ix_(members, members)]
        mean_dist = sub.mean(axis=1)
        medoids[int(label)] = int(members[int(np.argmin(mean_dist))])
    return medoids


def write_representative_structures(
    traj, medoids: Dict[int, int], labels, rmsd_matrix, frame_replica, frame_time_ps, pocket_volumes, out_dir: Path, *,
    show_progress: bool = True,
) -> List[ClusterReportRow]:
    """Write one representative structure per cluster (`cluster_00.pdb`
    largest cluster ... descending population), restricted to protein atoms
    (dropping solvent/ions when the sampling run used explicit solvent --
    `dd_docking`'s ensemble-member input is meant to be an apo receptor, not
    a receptor-plus-water-box; this is a no-op selection for implicit-
    solvent runs, which never had solvent atoms in the topology to begin
    with), plus `cluster_report.csv`."""
    import numpy as np
    import pandas as pd

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    protein_atoms = traj.topology.select("protein")

    populations = sorted(set(labels.tolist()), key=lambda label: -int((labels == label).sum()))
    progress = ClusterProgress(enabled=show_progress)
    rows: List[ClusterReportRow] = []
    for cluster_id, label in enumerate(populations):
        members = np.where(labels == label)[0]
        medoid_frame = medoids[label]
        mean_rmsd = float(rmsd_matrix[np.ix_(members, members)].mean())
        out_path = out_dir / f"cluster_{cluster_id:02d}.pdb"
        traj[medoid_frame].atom_slice(protein_atoms).save_pdb(str(out_path))
        rows.append(ClusterReportRow(
            cluster_id=cluster_id, n_frames=int(len(members)), source_replica=int(frame_replica[medoid_frame]),
            source_time_ps=float(frame_time_ps[medoid_frame]), mean_intra_cluster_rmsd_a=round(mean_rmsd, 3),
            pocket_volume_nm3=round(float(pocket_volumes[medoid_frame]), 3), out_path=str(out_path),
        ))
        progress.update(cluster_id, len(members), mean_rmsd, str(out_path))

    pd.DataFrame([r.__dict__ for r in rows]).to_csv(out_dir / "cluster_report.csv", index=False)
    return rows


def diagnostics_plot(rmsd_matrix, out_path: Path, *, n_clusters_range: Sequence[int] = range(2, 21)) -> None:
    """Silhouette score swept over `n_clusters_range`, purely informational
    (does not change the chosen `--n-clusters`)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import silhouette_score

    n_frames = rmsd_matrix.shape[0]
    xs, ys = [], []
    for k in n_clusters_range:
        if k < 2 or k >= n_frames:
            continue
        labels = cluster_frames(rmsd_matrix, k)
        xs.append(k)
        ys.append(silhouette_score(rmsd_matrix, labels, metric="precomputed"))

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(xs, ys, marker="o")
    ax.set_xlabel("n_clusters")
    ax.set_ylabel("silhouette score")
    ax.set_title("Cluster-count diagnostics (informational only)")
    fig.tight_layout()
    fig.savefig(str(out_path))
    plt.close(fig)


def cluster_pocket_trajectory(
    top_pdb: Path, dcd_paths: Sequence[Path], pocket_residues: Sequence[Residue], out_dir: Path, *,
    n_clusters: int = 10, cluster_atoms: str = "ca", diagnostics: bool = False,
    pocket_expand_only: bool = False, pocket_expand_margin: float = 0.0, show_progress: bool = True,
) -> List[ClusterReportRow]:
    """End-to-end: load + pool replica DCDs, compute the pocket-atom RMSD
    matrix, cluster into `n_clusters` groups, and write representative
    structures + `cluster_report.csv` (+ optional diagnostic plot).

    `pocket_expand_only` drops pooled frames whose `pocket_volume_proxy` has
    not grown (by at least a `pocket_expand_margin` fraction) relative to
    `top_pdb` -- the pre-MD reference structure -- *before* clustering, so
    every representative structure comes from a frame where the pocket
    opened up rather than closed down. This is a post-hoc selection over an
    otherwise-unbiased restrained-MD ensemble, not a steering force during
    sampling itself: `dd_af-sample` doesn't know about "pocket expansion" at
    all, so nothing about how the frames were generated is affected, only
    which of them are eligible to become a representative structure.

    `pocket_expand_margin=0.0` (the default) was measured to barely filter
    anything on a real sampling run (kept 10,000/10,000 pooled frames): a
    convex hull is an envelope over its points, so thermal noise alone
    tends to inflate it relative to any single static reference frame,
    regardless of any real "opening" signal. A positive margin (e.g.
    0.05-0.1) is needed for this flag to meaningfully narrow the pool --
    see the `[cluster] pocket-expand-only: kept X/Y frame(s)` line this
    prints to calibrate one for a given system.
    """
    import mdtraj as md
    import numpy as np

    # Atom indices are computed from `top_pdb` loaded on its own, not from
    # the pooled trajectory's topology: `mdtraj.join` (used below when
    # there is more than one replica) does not preserve each chain's
    # `chain_id` on the merged topology (it comes back `None`), which
    # would silently make every `(chain, resnum)` lookup fail. Topology is
    # identical across every replica (they all start from the same
    # `complex_top.pdb`), so computing it once, pre-join, is both correct
    # and cheaper.
    ref_traj = md.load(str(top_pdb))
    atom_indices = _pocket_atom_indices(ref_traj.topology, pocket_residues, atom_selection=cluster_atoms)
    if not atom_indices:
        raise ValueError("cluster_pocket_trajectory: no pocket atoms found for the given residues/selection")

    traj, frame_replica, frame_time_ps = load_pooled_trajectory(top_pdb, dcd_paths)
    volumes = pocket_volume_proxy(traj, atom_indices)

    if pocket_expand_only:
        reference_volume = float(pocket_volume_proxy(ref_traj, atom_indices)[0])
        keep = volumes >= reference_volume * (1.0 + pocket_expand_margin)
        n_kept = int(keep.sum())
        if show_progress:
            print(
                f"[cluster] pocket-expand-only: kept {n_kept}/{traj.n_frames} frame(s) with pocket volume "
                f">= {reference_volume * (1.0 + pocket_expand_margin):.2f} nm^3 "
                f"(reference {reference_volume:.2f} nm^3)", flush=True,
            )
        if n_kept == 0:
            raise ValueError(
                "cluster_pocket_trajectory: --pocket-expand-only left 0 frames -- no pooled frame's pocket "
                "grew past the reference volume (try a smaller --pocket-expand-margin, a longer --sample-ns, "
                "or disable --pocket-expand-only)"
            )
        traj = traj.slice(np.where(keep)[0])
        frame_replica = frame_replica[keep]
        frame_time_ps = frame_time_ps[keep]
        volumes = volumes[keep]

    rmsd_matrix = pairwise_rmsd_matrix(traj, atom_indices)
    labels = cluster_frames(rmsd_matrix, n_clusters)
    medoids = pick_medoids(rmsd_matrix, labels)
    rows = write_representative_structures(
        traj, medoids, labels, rmsd_matrix, frame_replica, frame_time_ps, volumes, out_dir, show_progress=show_progress,
    )

    total = sum(r.n_frames for r in rows)
    if total != traj.n_frames:
        raise AssertionError(f"cluster_report frame counts ({total}) != pooled frames ({traj.n_frames})")

    if diagnostics:
        diagnostics_plot(rmsd_matrix, Path(out_dir) / "cluster_diagnostics.png")

    return rows
