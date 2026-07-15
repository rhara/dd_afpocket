"""Command-line entry points:
  dd_af-fetch    UNIPROT -o out_dir
  dd_af-pocket   prepped.pdb -o out_dir
  dd_af-sample   prepped.pdb pocket_dir -o out_dir
  dd_af-cluster  sample_dir pocket_dir -o out_dir
  dd_af-run      UNIPROT -o out_dir
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict

from . import pipeline
from .pocket import parse_residue_list
from .sample import PROTEIN_FORCEFIELDS

# Sampling-length/replica-count/GB-model bundles for `--preset`, selected by
# `dd_af-sample`/`dd_af-run`. "default" spells out sample_pocket's own
# defaults so `--preset quick`'s savings are legible by contrast; "quick" is
# a coarse, CPU-only-friendly setting for when only rough pocket-shape
# diversity is needed for downstream ensemble docking, not a converged
# trajectory (see README "Performance"). timestep_fs is left at 4 fs in both
# -- that value, not the sampling length, is the one restraints.py's
# docstring says was empirically checked for stability against the default
# 1000 kJ/mol/nm^2 restraint force constant. quick's implicit_solvent="obc2"
# measured 1.37x faster than the default gbn2 on this project's development
# machine (~2300-atom protein, CPU/4 threads) while avoiding hct's known
# tendency to underestimate buried atoms' Born radii -- see
# sample.IMPLICIT_SOLVENT_FILES's comment for the full measured comparison.
PRESETS: Dict[str, Dict[str, Any]] = {
    "default": {"n_replicas": 4, "equil_ps": 20.0, "sample_ns": 2.0, "report_ps": 5.0, "timestep_fs": 4.0, "implicit_solvent": "gbn2"},
    "quick": {"n_replicas": 2, "equil_ps": 5.0, "sample_ns": 0.3, "report_ps": 2.5, "timestep_fs": 4.0, "implicit_solvent": "obc2"},
}
_PRESET_KEYS = ("n_replicas", "equil_ps", "sample_ns", "report_ps", "timestep_fs", "implicit_solvent")


def _add_preset_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--preset", choices=sorted(PRESETS), default="default",
        help="Bundle of sampling-length/replica-count/GB-model defaults ('quick': coarse, "
             "CPU-only-friendly sampling -- see README 'Performance'). Any of "
             "--n-replicas/--equil-ps/--sample-ns/--report-ps/--timestep-fs/--implicit-solvent given "
             "explicitly overrides the preset's value for that one flag.",
    )
    parser.add_argument("--n-replicas", type=int, default=None, help="Independent sampling replicas (default: 4, or preset's value)")
    parser.add_argument("--equil-ps", type=float, default=None, help="Equilibration length before production, in ps (default: 20.0, or preset's value)")
    parser.add_argument("--sample-ns", type=float, default=None, help="Production length per replica, in nanoseconds (default: 2.0, or preset's value)")
    parser.add_argument("--report-ps", type=float, default=None, help="Trajectory frame-writing interval, in ps (default: 5.0, or preset's value)")
    parser.add_argument("--timestep-fs", type=float, default=None, help="Integration timestep, in fs (default: 4.0, or preset's value; relies on constraints=HBonds + 4 amu hydrogen mass repartitioning, see sample.py)")
    parser.add_argument("--implicit-solvent", choices=["gbn2", "gbn", "obc2", "obc1", "hct"], default=None,
                         help="GB implicit-solvent model (default: gbn2, or preset's value; cheapest to most expensive: hct < obc1 < obc2 < gbn < gbn2 -- see sample.IMPLICIT_SOLVENT_FILES)")


def _resolve_preset_args(args: argparse.Namespace) -> Dict[str, Any]:
    """Individually-specified flags (parsed with `default=None` by
    `_add_preset_args`) win over `--preset`'s bundle; anything left unset
    falls back to the preset."""
    preset = PRESETS[args.preset]
    return {key: (getattr(args, key) if getattr(args, key) is not None else preset[key]) for key in _PRESET_KEYS}


# Water models pooled across every `PROTEIN_FORCEFIELDS` entry, for a single
# `choices=` list; `sample.resolve_water_file` raises a clear error at run
# time if a specific --protein-forcefield/--water-model combination wasn't
# one of the ones actually verified (see PROTEIN_FORCEFIELDS's comment).
_ALL_WATER_MODELS = sorted({w for entry in PROTEIN_FORCEFIELDS.values() for w in entry["water_files"]})


def _add_forcefield_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--protein-forcefield", choices=sorted(PROTEIN_FORCEFIELDS), default="amber14-all",
        help="Protein forcefield (default: amber14-all). Compatible with every --implicit-solvent GB model; "
             "--water-model support varies by forcefield (see sample.PROTEIN_FORCEFIELDS).",
    )
    parser.add_argument(
        "--solvent", choices=["implicit", "explicit"], default="implicit",
        help="implicit (default): GB continuum solvent, CutoffNonPeriodic, no periodic box -- fast, the "
             "setting the rest of this project's CPU-friendly defaults assume. explicit: a real periodic "
             "TIP3P/TIP4P-Ew/... water box + ions + PME + MonteCarloBarostat NPT -- much more expensive per "
             "step (a GPU, --platform CUDA, is strongly recommended over CPU for this).",
    )
    parser.add_argument(
        "--water-model", choices=_ALL_WATER_MODELS, default="tip3p",
        help="Explicit-solvent water model (default: tip3p; ignored when --solvent implicit). Must be one of "
             "the --protein-forcefield-specific choices in sample.PROTEIN_FORCEFIELDS.",
    )
    parser.add_argument(
        "--solvent-padding-nm", type=float, default=1.0,
        help="Minimum solute-to-box-edge padding for the periodic water box, in nm (default: 1.0; ignored when --solvent implicit)",
    )
    parser.add_argument(
        "--ion-concentration-molar", type=float, default=0.15,
        help="Neutralizing/background ion (Na+/Cl-) concentration, in mol/L (default: 0.15, physiological; ignored when --solvent implicit)",
    )
    parser.add_argument(
        "--pressure-atm", type=float, default=1.0,
        help="MonteCarloBarostat target pressure, in atm (default: 1.0; ignored when --solvent implicit)",
    )


def build_fetch_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_af-fetch",
        description="Fetch an AlphaFold DB model and repair it to MD grade (delegates to dd_prep).",
    )
    parser.add_argument("uniprot", help="UniProt accession, e.g. O60674")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory (raw/ and the prepped structure go here)")
    parser.add_argument("--ph", type=float, default=7.0, help="Protonation pH")
    parser.add_argument("--plddt-cutoff", type=float, default=50.0, help="pLDDT threshold below which N-/C-terminal tail residues are trimmed")
    parser.add_argument("--no-progress", action="store_true")
    return parser


def main_fetch(argv=None):
    args = build_fetch_parser().parse_args(argv)
    out_dir = Path(args.out_dir)
    raw_dir = out_dir / "raw"
    md_pdb = pipeline.fetch_and_prep(
        args.uniprot, raw_dir, out_dir, ph=args.ph, plddt_cutoff=args.plddt_cutoff,
        show_progress=not args.no_progress,
    )
    print(f"\n[done] {args.uniprot} -> {md_pdb}")


def build_pocket_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_af-pocket",
        description="Detect a druggable pocket with fpocket and write pocket_report.json/pocket_box.json.",
    )
    parser.add_argument("prepped_pdb", help="MD-grade prepped structure (dd_af-fetch or dd_prep-run --repair md output)")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("--pocket-rank", type=int, default=1, help="1-indexed pocket rank by Druggability Score descending")
    parser.add_argument("--pocket-residues", default=None, metavar="A:42,A:87,...", help="Bypass fpocket's residue detection for the selected pocket with a manual chain:resnum list")
    parser.add_argument("--box-padding", type=float, default=5.0, help="Docking box padding around the pocket-lining residues, in Angstrom")
    parser.add_argument("--no-progress", action="store_true")
    return parser


def main_pocket(argv=None):
    args = build_pocket_parser().parse_args(argv)
    residues = parse_residue_list(args.pocket_residues) if args.pocket_residues else None
    selection = pipeline.detect_pocket(
        args.prepped_pdb, args.out_dir, pocket_rank=args.pocket_rank, pocket_residues=residues,
        box_padding=args.box_padding, show_progress=not args.no_progress,
    )
    print(f"\n[done] pocket rank {selection.rank} (druggability={selection.druggability_score:.3f}, "
          f"{len(selection.residues)} residue(s)) -> {args.out_dir}/pocket_report.json")


def build_sample_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_af-sample",
        description="Restrained MD sampling of the pocket neighborhood (dd_af-pocket output required).",
    )
    parser.add_argument("prepped_pdb", help="MD-grade prepped structure (must match the one dd_af-pocket ran on)")
    parser.add_argument("pocket_dir", help="Directory containing pocket_report.json from dd_af-pocket")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    _add_preset_args(parser)
    _add_forcefield_args(parser)
    parser.add_argument("--n-jobs", type=int, default=1, help="Parallel worker processes across replicas (<=0 = all CPU cores, default 1 = sequential)")
    parser.add_argument("--mobile-radius-nm", type=float, default=1.0, help="Residues with a CA within this distance of the pocket centroid are left mobile, in addition to the fpocket-detected lining residues")
    parser.add_argument("--mobile-margin-residues", type=int, default=2, help="Sequence neighbors added on each side of every mobile residue")
    parser.add_argument("--force-constant", type=float, default=1.0e3, help="Position-restraint force constant for frozen residues, kJ/mol/nm^2")
    parser.add_argument("--progress-ps", type=float, default=100.0, help="In-run progress print interval")
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--nonbonded-cutoff-nm", type=float, default=1.5, help="Nonbonded interaction cutoff (CutoffNonPeriodic for implicit, PME for explicit); NoCutoff was measured impractically slow for large multi-domain proteins")
    parser.add_argument("--no-progress", action="store_true")
    return parser


def main_sample(argv=None):
    args = build_sample_parser().parse_args(argv)
    selection = pipeline.load_pocket_selection(args.pocket_dir)
    resolved = _resolve_preset_args(args)
    if not args.no_progress:
        print(f"[preset] {args.preset}: " + " ".join(f"{k}={v}" for k, v in resolved.items()), flush=True)
    result = pipeline.sample_pocket(
        args.prepped_pdb, selection["residues"], selection["center"], args.out_dir,
        n_jobs=args.n_jobs, mobile_radius_nm=args.mobile_radius_nm,
        mobile_margin_residues=args.mobile_margin_residues, force_constant_kj_per_mol_nm2=args.force_constant,
        progress_ps=args.progress_ps, platform_name=args.platform, nonbonded_cutoff_nm=args.nonbonded_cutoff_nm,
        protein_forcefield=args.protein_forcefield, solvent=args.solvent, water_model=args.water_model,
        solvent_padding_nm=args.solvent_padding_nm, ion_concentration_molar=args.ion_concentration_molar,
        pressure_atm=args.pressure_atm, show_progress=not args.no_progress, **resolved,
    )
    print(f"\n[done] {result.n_replicas} replica(s), {result.n_frames_total} pooled frame(s) -> {args.out_dir}/restraint_report.json")


def build_cluster_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_af-cluster",
        description="Cluster pooled sampling replicas into N representative pocket conformations.",
    )
    parser.add_argument("sample_dir", help="Directory containing complex_top.pdb/sample_r*.dcd from dd_af-sample")
    parser.add_argument("pocket_dir", help="Directory containing pocket_report.json from dd_af-pocket")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory for cluster_NN.pdb / cluster_report.csv")
    parser.add_argument("--n-clusters", type=int, default=10, help="Number of representative structures to produce")
    parser.add_argument("--cluster-atoms", default="ca", choices=["ca", "heavy"], help="Pocket atom selection used for the RMSD distance matrix")
    parser.add_argument("--diagnostics", action="store_true", help="Also write cluster_diagnostics.png (silhouette score vs n_clusters, informational only)")
    parser.add_argument(
        "--pocket-expand-only", action="store_true",
        help="Only cluster pooled frames whose pocket_volume_proxy (convex hull over --cluster-atoms) is "
             ">= the pre-MD reference structure's -- i.e. every representative structure comes from a frame "
             "where the pocket opened up, never one where it closed down. A post-hoc selection over the "
             "sampled ensemble, not a bias applied during dd_af-sample itself.",
    )
    parser.add_argument(
        "--pocket-expand-margin", type=float, default=0.0,
        help="With --pocket-expand-only, require pocket volume >= reference * (1 + this fraction) "
             "(default: 0.0). A convex hull is an envelope over its points, so thermal noise alone tends to "
             "inflate it relative to any single static reference frame -- margin 0.0 was measured to keep "
             "essentially every pooled frame on a real run (see README 'Post-hoc pocket-expansion filtering'). "
             "Start around 0.05-0.1 and look at the printed 'kept X/Y frame(s)' line to calibrate.",
    )
    parser.add_argument(
        "--visualize", action="store_true",
        help="Also write cluster_overlay.html: a standalone (no server) py3Dmol overlay of every "
             "cluster_NN.pdb, colored per cluster, for visually comparing pocket-shape diversity",
    )
    parser.add_argument("--no-progress", action="store_true")
    return parser


def main_cluster(argv=None):
    args = build_cluster_parser().parse_args(argv)
    sample_dir = Path(args.sample_dir)
    selection = pipeline.load_pocket_selection(args.pocket_dir)
    dcds = sorted(sample_dir.glob("sample_r*.dcd"))
    if not dcds:
        raise SystemExit(f"dd_af-cluster: no sample_r*.dcd found in {sample_dir}")

    rows = pipeline.cluster_pocket(
        sample_dir / "complex_top.pdb", dcds, selection["residues"], args.out_dir,
        n_clusters=args.n_clusters, cluster_atoms=args.cluster_atoms, diagnostics=args.diagnostics,
        pocket_expand_only=args.pocket_expand_only, pocket_expand_margin=args.pocket_expand_margin,
        show_progress=not args.no_progress,
    )
    print(f"\n[done] {len(rows)} representative structure(s) -> {args.out_dir}/cluster_report.csv")
    if args.visualize:
        html_path = pipeline.visualize_clusters(args.out_dir, selection["residues"])
        print(f"[done] cluster overlay -> {html_path}")


def build_run_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_af-run",
        description="End-to-end: fetch AFDB model -> detect druggable pocket -> restrained MD sampling -> cluster into N representative structures.",
    )
    parser.add_argument("uniprot", help="UniProt accession, e.g. O60674")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory (a <uniprot>/ subdirectory is created)")
    parser.add_argument("--pocket-rank", type=int, default=1)
    parser.add_argument("--box-padding", type=float, default=5.0)
    _add_preset_args(parser)
    _add_forcefield_args(parser)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--mobile-radius-nm", type=float, default=1.0)
    parser.add_argument("--mobile-margin-residues", type=int, default=2)
    parser.add_argument("--force-constant", type=float, default=1.0e3)
    parser.add_argument("--progress-ps", type=float, default=100.0)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--nonbonded-cutoff-nm", type=float, default=1.5)
    parser.add_argument("--n-clusters", type=int, default=10)
    parser.add_argument("--cluster-atoms", default="ca", choices=["ca", "heavy"])
    parser.add_argument("--diagnostics", action="store_true")
    parser.add_argument("--pocket-expand-only", action="store_true")
    parser.add_argument("--pocket-expand-margin", type=float, default=0.0)
    parser.add_argument("--visualize", action="store_true")
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--plddt-cutoff", type=float, default=50.0)
    parser.add_argument("--no-progress", action="store_true")
    return parser


def main_run(argv=None):
    args = build_run_parser().parse_args(argv)
    resolved = _resolve_preset_args(args)
    if not args.no_progress:
        print(f"[preset] {args.preset}: " + " ".join(f"{k}={v}" for k, v in resolved.items()), flush=True)
    result = pipeline.run_end_to_end(
        args.uniprot, args.out_dir, pocket_rank=args.pocket_rank, box_padding=args.box_padding,
        n_jobs=args.n_jobs, mobile_radius_nm=args.mobile_radius_nm,
        mobile_margin_residues=args.mobile_margin_residues, force_constant_kj_per_mol_nm2=args.force_constant,
        progress_ps=args.progress_ps, platform_name=args.platform, nonbonded_cutoff_nm=args.nonbonded_cutoff_nm,
        protein_forcefield=args.protein_forcefield, solvent=args.solvent, water_model=args.water_model,
        solvent_padding_nm=args.solvent_padding_nm, ion_concentration_molar=args.ion_concentration_molar,
        pressure_atm=args.pressure_atm, n_clusters=args.n_clusters, cluster_atoms=args.cluster_atoms,
        diagnostics=args.diagnostics, pocket_expand_only=args.pocket_expand_only,
        pocket_expand_margin=args.pocket_expand_margin, visualize=args.visualize,
        ph=args.ph, plddt_cutoff=args.plddt_cutoff, show_progress=not args.no_progress, **resolved,
    )
    print(f"\n[done] {len(result['clusters'])} representative structure(s) -> {args.out_dir}/{args.uniprot.lower()}/clusters/cluster_report.csv")
    if result.get("cluster_overlay_html"):
        print(f"[done] cluster overlay -> {result['cluster_overlay_html']}")
