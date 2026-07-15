"""End-to-end orchestration of the four stages (fetch -> pocket -> sample ->
cluster), plus each stage as an independently callable function so the
`dd_af-fetch`/`dd_af-pocket`/`dd_af-sample`/`dd_af-cluster` console commands
can run any one of them on its own.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from . import cluster as cluster_mod
from . import pocket as pocket_mod
from . import sample as sample_mod
from . import visualize as visualize_mod
from .pocket import PocketSelection, Residue


def fetch_and_prep(uniprot_id: str, raw_dir: Path, out_dir: Path, *, ph: float = 7.0,
                   plddt_cutoff: float = 50.0, show_progress: bool = True) -> str:
    """UniProt accession -> AFDB fetch + MD-grade structure repair,
    delegated to `dd_prep.pipeline.fetch_and_prepare_afdb` (see PROMPT for
    why dd_af depends on dd_prep here rather than reimplementing structure
    repair: AlphaFold models have none of the real-PDB-deposit quirks
    dd_md's self-contained `receptor_prep.py` exists to handle). Returns
    the path to the protonated `<uniprot>_afdb_raw_md.pdb` output.
    """
    from dd_prep.pipeline import fetch_and_prepare_afdb

    report = fetch_and_prepare_afdb(
        uniprot_id, Path(raw_dir), Path(out_dir),
        repair_mode="md", ph=ph, plddt_cutoff=plddt_cutoff, show_progress=show_progress,
    )
    if not report.md_pdb:
        raise RuntimeError(f"dd_prep produced no MD-grade output for {uniprot_id}")
    return report.md_pdb


def detect_pocket(
    prepped_pdb: Path, work_dir: Path, *, pocket_rank: int = 1,
    pocket_residues: Optional[Sequence[Residue]] = None, box_padding: float = 5.0,
    show_progress: bool = True,
) -> PocketSelection:
    """Run fpocket, select a pocket, and write `pocket_report.json` /
    `pocket_box.json` under `work_dir`."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    selection = pocket_mod.find_druggable_pocket(
        prepped_pdb, work_dir, pocket_rank=pocket_rank, pocket_residues=pocket_residues,
        box_padding=box_padding, show_progress=show_progress,
    )
    (work_dir / "pocket_report.json").write_text(json.dumps(selection.to_report_dict(), indent=2))
    (work_dir / "pocket_box.json").write_text(json.dumps(selection.to_box_dict(), indent=2))
    return selection


def load_pocket_selection(work_dir: Path) -> Dict[str, Any]:
    """Read back `pocket_report.json`/`pocket_box.json` written by
    `detect_pocket`, returning `{"residues": [Residue, ...], "center":
    (x, y, z), "box_center": [...], "box_size": [...]}`."""
    work_dir = Path(work_dir)
    report = json.loads((work_dir / "pocket_report.json").read_text())
    box = json.loads((work_dir / "pocket_box.json").read_text())
    return {
        "residues": pocket_mod.parse_residue_list(",".join(report["residues"])),
        "center": tuple(report["center"]),
        "box_center": box["center"],
        "box_size": box["size"],
    }


def sample_pocket(
    prepped_pdb: Path, pocket_residues: Sequence[Residue], pocket_center_a: Sequence[float], out_dir: Path,
    **kwargs: Any,
) -> sample_mod.SampleResult:
    """Restrained-MD sampling of the pocket neighborhood -- thin
    pass-through to `sample.sample_pocket` (see there for all keyword
    arguments: `n_replicas`, `n_jobs`, `mobile_radius_nm`,
    `mobile_margin_residues`, `force_constant_kj_per_mol_nm2`, `equil_ps`,
    `sample_ns`, `report_ps`, `progress_ps`, `platform_name`,
    `nonbonded_cutoff_nm`, `show_progress`)."""
    return sample_mod.sample_pocket(prepped_pdb, pocket_residues, pocket_center_a, out_dir, **kwargs)


def cluster_pocket(
    complex_top_pdb: Path, replica_dcds: Sequence[Path], pocket_residues: Sequence[Residue], out_dir: Path, *,
    n_clusters: int = 10, cluster_atoms: str = "ca", diagnostics: bool = False,
    pocket_expand_only: bool = False, pocket_expand_margin: float = 0.0, show_progress: bool = True,
) -> List[cluster_mod.ClusterReportRow]:
    """Structural clustering into representative structures (see
    `cluster.cluster_pocket_trajectory`)."""
    return cluster_mod.cluster_pocket_trajectory(
        complex_top_pdb, replica_dcds, pocket_residues, out_dir,
        n_clusters=n_clusters, cluster_atoms=cluster_atoms, diagnostics=diagnostics,
        pocket_expand_only=pocket_expand_only, pocket_expand_margin=pocket_expand_margin, show_progress=show_progress,
    )


def visualize_clusters(out_dir: Path, pocket_residues: Sequence[Residue]) -> str:
    """Overlay every `cluster_NN.pdb` in `out_dir` into one standalone
    `cluster_overlay.html` (see `visualize.write_cluster_overlay_html`)."""
    out_dir = Path(out_dir)
    cluster_pdbs = sorted(out_dir.glob("cluster_*.pdb"))
    if not cluster_pdbs:
        raise FileNotFoundError(f"visualize_clusters: no cluster_*.pdb found in {out_dir}")
    html_path = visualize_mod.write_cluster_overlay_html(cluster_pdbs, pocket_residues, out_dir / "cluster_overlay.html")
    return str(html_path)


def run_end_to_end(
    uniprot_id: str, out_dir: Path, *, pocket_rank: int = 1, box_padding: float = 5.0,
    n_replicas: int = 4, n_jobs: int = 1, mobile_radius_nm: float = 1.0, mobile_margin_residues: int = 2,
    force_constant_kj_per_mol_nm2: float = 1.0e3, equil_ps: float = 20.0, sample_ns: float = 2.0,
    report_ps: float = 5.0, timestep_fs: float = 4.0, progress_ps: float = 100.0, platform_name: str = "CPU",
    nonbonded_cutoff_nm: float = 1.5, protein_forcefield: str = "amber14-all", solvent: str = "implicit",
    implicit_solvent: str = "gbn2", water_model: str = "tip3p", solvent_padding_nm: float = 1.0,
    ion_concentration_molar: float = 0.15, pressure_atm: float = 1.0, n_clusters: int = 10, cluster_atoms: str = "ca",
    diagnostics: bool = False, pocket_expand_only: bool = False, pocket_expand_margin: float = 0.0,
    visualize: bool = False, ph: float = 7.0, plddt_cutoff: float = 50.0, show_progress: bool = True,
) -> Dict[str, Any]:
    """`dd_af-run`: fetch -> pocket -> sample -> cluster, end-to-end, under
    `<out_dir>/<uniprot_id_lower>/`."""
    out_dir = Path(out_dir) / uniprot_id.lower()
    raw_dir = out_dir / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)

    prepped_pdb = fetch_and_prep(uniprot_id, raw_dir, out_dir, ph=ph, plddt_cutoff=plddt_cutoff, show_progress=show_progress)
    selection = detect_pocket(prepped_pdb, out_dir, pocket_rank=pocket_rank, box_padding=box_padding, show_progress=show_progress)
    sample_result = sample_pocket(
        prepped_pdb, selection.residues, selection.center, out_dir,
        n_replicas=n_replicas, n_jobs=n_jobs, mobile_radius_nm=mobile_radius_nm,
        mobile_margin_residues=mobile_margin_residues, force_constant_kj_per_mol_nm2=force_constant_kj_per_mol_nm2,
        equil_ps=equil_ps, sample_ns=sample_ns, report_ps=report_ps, timestep_fs=timestep_fs, progress_ps=progress_ps,
        platform_name=platform_name, nonbonded_cutoff_nm=nonbonded_cutoff_nm, protein_forcefield=protein_forcefield,
        solvent=solvent, implicit_solvent=implicit_solvent, water_model=water_model,
        solvent_padding_nm=solvent_padding_nm, ion_concentration_molar=ion_concentration_molar,
        pressure_atm=pressure_atm, show_progress=show_progress,
    )
    cluster_rows = cluster_pocket(
        sample_result.complex_top_pdb, sample_result.replica_dcds, selection.residues, out_dir / "clusters",
        n_clusters=n_clusters, cluster_atoms=cluster_atoms, diagnostics=diagnostics,
        pocket_expand_only=pocket_expand_only, pocket_expand_margin=pocket_expand_margin, show_progress=show_progress,
    )
    overlay_html = visualize_clusters(out_dir / "clusters", selection.residues) if visualize else None

    return {
        "prepped_pdb": prepped_pdb, "pocket": selection.to_report_dict(), "box": selection.to_box_dict(),
        "sample": sample_result.to_dict(), "clusters": [row.__dict__ for row in cluster_rows],
        "cluster_overlay_html": overlay_html,
    }
