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

from . import pipeline
from .pocket import parse_residue_list


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
    parser.add_argument("--n-replicas", type=int, default=4, help="Independent sampling replicas")
    parser.add_argument("--n-jobs", type=int, default=1, help="Parallel worker processes across replicas (<=0 = all CPU cores, default 1 = sequential)")
    parser.add_argument("--mobile-radius-nm", type=float, default=1.0, help="Residues with a CA within this distance of the pocket centroid are left mobile, in addition to the fpocket-detected lining residues")
    parser.add_argument("--mobile-margin-residues", type=int, default=2, help="Sequence neighbors added on each side of every mobile residue")
    parser.add_argument("--force-constant", type=float, default=1.0e3, help="Position-restraint force constant for frozen residues, kJ/mol/nm^2")
    parser.add_argument("--equil-ps", type=float, default=20.0)
    parser.add_argument("--sample-ns", type=float, default=2.0, help="Production length per replica, in nanoseconds")
    parser.add_argument("--report-ps", type=float, default=5.0, help="Trajectory frame-writing interval")
    parser.add_argument("--progress-ps", type=float, default=100.0, help="In-run progress print interval")
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--nonbonded-cutoff-nm", type=float, default=1.5, help="Nonbonded interaction cutoff (CutoffNonPeriodic); NoCutoff was measured impractically slow for large multi-domain proteins")
    parser.add_argument("--no-progress", action="store_true")
    return parser


def main_sample(argv=None):
    args = build_sample_parser().parse_args(argv)
    selection = pipeline.load_pocket_selection(args.pocket_dir)
    result = pipeline.sample_pocket(
        args.prepped_pdb, selection["residues"], selection["center"], args.out_dir,
        n_replicas=args.n_replicas, n_jobs=args.n_jobs, mobile_radius_nm=args.mobile_radius_nm,
        mobile_margin_residues=args.mobile_margin_residues, force_constant_kj_per_mol_nm2=args.force_constant,
        equil_ps=args.equil_ps, sample_ns=args.sample_ns, report_ps=args.report_ps, progress_ps=args.progress_ps,
        platform_name=args.platform, nonbonded_cutoff_nm=args.nonbonded_cutoff_nm, show_progress=not args.no_progress,
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
        show_progress=not args.no_progress,
    )
    print(f"\n[done] {len(rows)} representative structure(s) -> {args.out_dir}/cluster_report.csv")


def build_run_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_af-run",
        description="End-to-end: fetch AFDB model -> detect druggable pocket -> restrained MD sampling -> cluster into N representative structures.",
    )
    parser.add_argument("uniprot", help="UniProt accession, e.g. O60674")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory (a <uniprot>/ subdirectory is created)")
    parser.add_argument("--pocket-rank", type=int, default=1)
    parser.add_argument("--box-padding", type=float, default=5.0)
    parser.add_argument("--n-replicas", type=int, default=4)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--mobile-radius-nm", type=float, default=1.0)
    parser.add_argument("--mobile-margin-residues", type=int, default=2)
    parser.add_argument("--force-constant", type=float, default=1.0e3)
    parser.add_argument("--equil-ps", type=float, default=20.0)
    parser.add_argument("--sample-ns", type=float, default=2.0)
    parser.add_argument("--report-ps", type=float, default=5.0)
    parser.add_argument("--progress-ps", type=float, default=100.0)
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--nonbonded-cutoff-nm", type=float, default=1.5)
    parser.add_argument("--n-clusters", type=int, default=10)
    parser.add_argument("--cluster-atoms", default="ca", choices=["ca", "heavy"])
    parser.add_argument("--diagnostics", action="store_true")
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--plddt-cutoff", type=float, default=50.0)
    parser.add_argument("--no-progress", action="store_true")
    return parser


def main_run(argv=None):
    args = build_run_parser().parse_args(argv)
    result = pipeline.run_end_to_end(
        args.uniprot, args.out_dir, pocket_rank=args.pocket_rank, box_padding=args.box_padding,
        n_replicas=args.n_replicas, n_jobs=args.n_jobs, mobile_radius_nm=args.mobile_radius_nm,
        mobile_margin_residues=args.mobile_margin_residues, force_constant_kj_per_mol_nm2=args.force_constant,
        equil_ps=args.equil_ps, sample_ns=args.sample_ns, report_ps=args.report_ps, progress_ps=args.progress_ps,
        platform_name=args.platform, nonbonded_cutoff_nm=args.nonbonded_cutoff_nm,
        n_clusters=args.n_clusters, cluster_atoms=args.cluster_atoms, diagnostics=args.diagnostics,
        ph=args.ph, plddt_cutoff=args.plddt_cutoff, show_progress=not args.no_progress,
    )
    print(f"\n[done] {len(result['clusters'])} representative structure(s) -> {args.out_dir}/{args.uniprot.lower()}/clusters/cluster_report.csv")
