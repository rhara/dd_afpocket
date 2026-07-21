"""Local restrained-MD sampling: build a protein-only system (implicit
solvent by default, optionally solvated in an explicit periodic water box)
from a prepped apo AlphaFold structure, freeze everything except the
detected pocket's neighborhood (`restraints.py`), and run one or more
independent replicas whose pooled frames feed `cluster.py`.

Plain `openmm.app.ForceField` is used to build the system rather than
`openmmforcefields.generators.SystemGenerator` (the pattern
`dd_docking/refine_md.py` and `dd_md/system_build.py` use): both of those
projects need `SystemGenerator` because their systems include a small-
molecule ligand that needs GAFF/SMIRNOFF parameterization. dd_afpocket's systems
never have a ligand -- these are apo AlphaFold models -- so a plain
`ForceField(...)` is all that is needed, and it happens to sidestep a real
environment issue: constructing a `SystemGenerator` without an explicit
`small_molecule_forcefield="gaff-*"` eagerly triggers `openff.toolkit`'s
SMIRNOFF force-field discovery, which in the current `mpro` env fails with
`ModuleNotFoundError: No module named 'pkg_resources'` (an
`openff-amber-ff-ports` dependency on the `pkg_resources` API that
setuptools >= 81 no longer ships). Confirmed this does not affect
dd_docking/dd_md, since both always pass `small_molecule_forcefield=
"gaff-2.11"` explicitly, which takes a different (GAFF) template-generator
code path that never reaches the broken import.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Sequence

from .pocket import Residue
from .progress import ReplicaProgress, StepProgress

# Every combo below was empirically confirmed to build (`ForceField(...)
# .createSystem(...)`, implicit; `Modeller.addSolvent()` +
# `.createSystem(...)`, explicit) on this project's development machine
# before being added here -- OpenMM ships many more protein/water XML
# files than are listed, but an untested pairing (e.g. an AMBER GB file
# against CHARMM atom types) can fail outright or build a silently
# mismatched system, so the choices exposed here are deliberately a
# curated, verified subset rather than "everything OpenMM has".
PROTEIN_FORCEFIELDS: Dict[str, Dict[str, Any]] = {
    "amber14-all": {
        "files": ["amber14-all.xml"],
        "water_files": {
            "tip3p": "amber14/tip3p.xml", "tip3pfb": "amber14/tip3pfb.xml",
            "tip4pew": "amber14/tip4pew.xml", "tip4pfb": "amber14/tip4pfb.xml",
            "spce": "amber14/spce.xml",
        },
    },
    "amber99sbildn": {
        "files": ["amber99sbildn.xml"],
        "water_files": {"tip3p": "tip3p.xml", "tip4pew": "tip4pew.xml", "tip5p": "tip5p.xml"},
    },
    "amber19-all": {
        "files": ["amber19-all.xml"],
        "water_files": {"tip3p": "amber19/tip3p.xml", "tip3pfb": "amber19/tip3pfb.xml", "tip4pew": "amber19/tip4pew.xml"},
    },
    "charmm36": {
        "files": ["charmm36.xml"],
        "water_files": {"tip3p": "charmm36/tip3p-pme-b.xml", "tip4pew": "charmm36/tip4pew.xml", "tip5p": "charmm36/tip5p.xml"},
    },
}
# `Modeller.addSolvent`'s `model=` argument selects water *geometry*
# (site count/placement) and is independent of which parameter xml is
# loaded: "fb"-suffixed (force-balance-reparametrized) files reuse their
# unsuffixed sibling's geometry, and charmm36's own "-pme-b" TIP3P variant
# reuses the plain tip3p geometry too.
WATER_MODEL_GEOMETRY = {
    "tip3p": "tip3p", "tip3pfb": "tip3p",
    "tip4pew": "tip4pew", "tip4pfb": "tip4pew",
    "tip5p": "tip5p", "spce": "spce",
}
# GB implicit-solvent variants, cheapest to most expensive (all measured on
# this project's development machine, ~2300-atom protein, CPU/4 threads,
# ms/step relative to gbn2): hct 1.75x faster, obc1 1.55x, obc2 1.37x, gbn
# ~1.02x (gbn2 adds a "neck" correction term to gbn's Born-radius sum, which
# is nearly all of gbn2's extra cost over gbn). gbn2 remains the default for
# backward compatibility and its explicit stability check against the
# default restraint force constant (see restraints.py); `--preset quick`
# picks obc2 as a faster-but-still-robust alternative (its Born-radius
# correction avoids hct's known tendency to underestimate buried atoms'
# radii, which matters here since some pocket-lining side chains are
# partially buried). All five were confirmed to build against every
# `PROTEIN_FORCEFIELDS` entry above, including charmm36 (OpenMM's generic
# `implicit/*.xml` GB parameters are keyed by element/class, not by a
# specific protein forcefield's atom-type naming).
IMPLICIT_SOLVENT_FILES = {
    "gbn2": "implicit/gbn2.xml",
    "gbn": "implicit/gbn.xml",
    "obc2": "implicit/obc2.xml",
    "obc1": "implicit/obc1.xml",
    "hct": "implicit/hct.xml",
}


def resolve_protein_forcefield_files(protein_forcefield: str) -> List[str]:
    entry = PROTEIN_FORCEFIELDS.get(protein_forcefield)
    if entry is None:
        raise ValueError(f"Unknown protein forcefield {protein_forcefield!r}; choose one of {sorted(PROTEIN_FORCEFIELDS)}")
    return list(entry["files"])


def resolve_water_file(protein_forcefield: str, water_model: str) -> str:
    entry = PROTEIN_FORCEFIELDS.get(protein_forcefield)
    water_files = entry["water_files"] if entry else {}
    if water_model not in water_files:
        raise ValueError(
            f"--water-model {water_model!r} is not a verified combination with --protein-forcefield "
            f"{protein_forcefield!r}; available for this forcefield: {sorted(water_files)}"
        )
    return water_files[water_model]


@dataclass
class SampleResult:
    complex_top_pdb: str
    replica_dcds: List[str]
    n_replicas: int
    n_frames_total: int
    residues_frozen: int
    residues_mobile: int
    force_constant_kj_per_mol_nm2: float
    mobile_radius_nm: float
    sample_ns: float
    equil_ps: float
    report_ps: float
    timestep_fs: float
    protein_forcefield: str
    solvent: str
    implicit_solvent: str
    water_model: str
    platform: str
    wall_seconds: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return dict(self.__dict__)


def _solvate_and_write(
    protein_pdb: Path, out_path: Path, *, protein_forcefield: str = "amber14-all", water_model: str = "tip3p",
    padding_nm: float = 1.0, ion_concentration_molar: float = 0.15,
) -> Path:
    """Solvate `protein_pdb` once, writing the periodic (CRYST1 record +
    every atom, protein/water/ion alike) structure to `out_path`.

    Must run exactly once and be shared by every replica: `Modeller.
    addSolvent` has no seed parameter, so calling it independently per
    replica (the way each replica's *implicit*-solvent system build calls
    `ForceField.createSystem` independently -- deterministic, no problem
    there) would give each replica's own water/ion placement, and hence a
    different atom count/ordering than whatever `complex_top.pdb`
    `cluster.py` ends up loading every replica's DCD against. Solvating
    once here and having every replica load the resulting file via
    `PDBFile` (as plain, deterministic input, exactly like the apo case)
    sidesteps that entirely.

    `PDBFile` round-trips the periodic box vectors losslessly through the
    written CRYST1 record (confirmed empirically), so no separate
    box-vector bookkeeping is needed once this file exists.
    """
    from openmm import unit
    from openmm.app import ForceField, Modeller, PDBFile

    protein_files = resolve_protein_forcefield_files(protein_forcefield)
    water_file = resolve_water_file(protein_forcefield, water_model)
    ff = ForceField(*protein_files, water_file)

    pdb = PDBFile(str(protein_pdb))
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(
        ff, model=WATER_MODEL_GEOMETRY[water_model], padding=padding_nm * unit.nanometer,
        ionicStrength=ion_concentration_molar * unit.molar,
    )
    out_path = Path(out_path)
    with open(out_path, "w") as fh:
        # keepIds=True: restraints.py/pocket.py match residues by
        # (chain id, resSeq) parsed from the PDB text -- addSolvent doesn't
        # touch the protein's existing chains/residues, but writing without
        # keepIds would still renumber everything from scratch, breaking
        # that matching for no reason.
        PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
    return out_path


def _build_system(protein_pdb: Path, *, platform_name: str = "CPU", nonbonded_cutoff_nm: float = 1.5,
                   cpu_threads: int = 0, protein_forcefield: str = "amber14-all", solvent: str = "implicit",
                   implicit_solvent: str = "gbn2", water_model: str = "tip3p"):
    """Load `protein_pdb` (already protonated -- dd_prep's `*_md.pdb`
    output, or, for `solvent="explicit"`, `_solvate_and_write`'s already-
    periodic output) and build an OpenMM system for it. Returns (pdbfile,
    system, platform, platform_name_used, properties).

    `solvent="implicit"` (default) uses `CutoffNonPeriodic` (default
    1.5 nm), not `NoCutoff`: every GB model's own Born-radius/GB-energy
    terms already approximate the far-field electrostatic effect of atoms
    beyond the cutoff, but the *pairwise* nonbonded loop itself is a naive
    O(n_atoms^2) evaluation regardless -- for a large multi-domain protein
    (e.g. full-length JAK2 at 1132 residues / ~17800 atoms with hydrogens)
    `NoCutoff` was measured to make even a few picoseconds impractically
    slow on CPU. Since restraints.py freezes everything outside a ~1 nm
    pocket neighborhood anyway, a 1.5 nm nonbonded cutoff costs essentially
    no accuracy for the region this project actually samples.
    `implicit_solvent` selects the GB variant (`IMPLICIT_SOLVENT_FILES`,
    default `"gbn2"`).

    `solvent="explicit"` uses `PME` (the periodic box vectors come from
    `protein_pdb`'s own CRYST1 record) and `rigidWater=True` (explicit
    water is normally kept rigid rather than hydrogen-mass-repartitioned --
    HMR only matters for the protein's own 4 fs timestep stability here).

    `protein_forcefield` selects the protein XML(s) (`PROTEIN_FORCEFIELDS`).

    `cpu_threads` (0 = let OpenMM choose, its own default) matters when
    several replicas run concurrently via `--n-jobs`: the CPU platform's
    own default thread count is normally "all logical cores", so N
    processes each independently trying to claim every core oversubscribes
    the machine and every replica gets slower, not faster (the same
    "workers don't oversubscribe" concern `dd_docking/screening.py`
    documents for Vina; `sample_pocket` computes a per-replica thread count
    from `os.cpu_count() // n_workers` the same way).
    """
    from openmm import Platform, app, unit
    from openmm.app import ForceField, PDBFile

    pdb = PDBFile(str(protein_pdb))
    protein_files = resolve_protein_forcefield_files(protein_forcefield)
    if solvent == "explicit":
        ff = ForceField(*protein_files, resolve_water_file(protein_forcefield, water_model))
        nonbonded_method, rigid_water = app.PME, True
    else:
        ff = ForceField(*protein_files, IMPLICIT_SOLVENT_FILES[implicit_solvent])
        nonbonded_method, rigid_water = app.CutoffNonPeriodic, False

    system = ff.createSystem(
        pdb.topology, constraints=app.HBonds, rigidWater=rigid_water, hydrogenMass=4 * unit.amu,
        nonbondedMethod=nonbonded_method, nonbondedCutoff=nonbonded_cutoff_nm * unit.nanometer,
    )
    try:
        platform = Platform.getPlatformByName(platform_name)
        used = platform_name
    except Exception:  # noqa: BLE001
        platform = Platform.getPlatformByName("CPU")
        used = "CPU"
    properties = {"Threads": str(cpu_threads)} if (used == "CPU" and cpu_threads > 0) else {}
    return pdb, system, platform, used, properties


def run_one_replica(
    protein_pdb: Path, pocket_residues: Sequence[Residue], pocket_center_a: Sequence[float], workdir: Path, *,
    replica: int = 1, temperature_k: float = 300.0, timestep_fs: float = 4.0,
    mobile_radius_nm: float = 1.0, mobile_margin_residues: int = 2,
    force_constant_kj_per_mol_nm2: float = 1.0e3, equil_ps: float = 20.0, sample_ns: float = 2.0,
    report_ps: float = 5.0, progress_ps: float = 100.0, platform_name: str = "CPU",
    nonbonded_cutoff_nm: float = 1.5, protein_forcefield: str = "amber14-all", solvent: str = "implicit",
    implicit_solvent: str = "gbn2", water_model: str = "tip3p", pressure_atm: float = 1.0,
    cpu_threads: int = 0, show_progress: bool = True,
) -> Dict[str, Any]:
    """Minimize -> equilibrate -> one production segment for a single
    replica, writing `workdir/complex_top.pdb` (shared across replicas,
    rewritten identically each time) and `workdir/sample_r{replica}.dcd`.

    A different random seed is used per replica (via
    `setVelocitiesToTemperature`'s own RNG, reseeded per replica through
    the integrator) so independent replicas actually diverge instead of
    retracing the same trajectory.

    `solvent="explicit"` adds a `MonteCarloBarostat` (target `pressure_atm`,
    default 1 atm) for NPT-style pressure coupling during both equilibration
    and production; `protein_pdb` must already be `_solvate_and_write`'s
    periodic output in that case (`sample_pocket` handles this).
    """
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, app, unit
    from openmm.app import Simulation

    from .restraints import freeze_non_pocket_residues, mobile_residue_indices

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    pdb, system, platform, platform_used, properties = _build_system(
        protein_pdb, platform_name=platform_name, nonbonded_cutoff_nm=nonbonded_cutoff_nm,
        cpu_threads=cpu_threads, protein_forcefield=protein_forcefield, solvent=solvent,
        implicit_solvent=implicit_solvent, water_model=water_model,
    )
    if solvent == "explicit":
        system.addForce(MonteCarloBarostat(pressure_atm * unit.atmosphere, temperature_k * unit.kelvin))

    center_nm = tuple(c / 10.0 for c in pocket_center_a)
    mobile = mobile_residue_indices(
        pdb.topology, pdb.positions, pocket_residues, center_nm,
        mobile_radius_nm=mobile_radius_nm, mobile_margin_residues=mobile_margin_residues,
    )
    counts = freeze_non_pocket_residues(
        system, pdb.topology, pdb.positions, mobile,
        force_constant_kj_per_mol_nm2=force_constant_kj_per_mol_nm2,
    )

    integrator = LangevinMiddleIntegrator(
        temperature_k * unit.kelvin, 1.0 / unit.picosecond, timestep_fs * unit.femtoseconds,
    )
    integrator.setRandomNumberSeed(1000 + replica)
    sim = Simulation(pdb.topology, system, integrator, platform, properties)
    sim.context.setPositions(pdb.positions)

    top_pdb = workdir / "complex_top.pdb"
    if not top_pdb.exists():
        with open(top_pdb, "w") as fh:
            app.PDBFile.writeFile(pdb.topology, pdb.positions, fh)

    sim.minimizeEnergy(maxIterations=500)

    dt_ps = timestep_fs / 1000.0
    sim.context.setVelocitiesToTemperature(temperature_k * unit.kelvin, 1000 + replica)

    n_progress = max(1, round(progress_ps / dt_ps))
    n_equil = max(1, round(equil_ps / dt_ps))
    sim.reporters.append(StepProgress(f"r{replica} equil", n_progress, dt_ps, n_equil, enabled=show_progress))
    sim.step(n_equil)
    sim.reporters.clear()

    dcd = workdir / f"sample_r{replica}.dcd"
    n_report = max(1, round(report_ps / dt_ps))
    n_prod = max(1, round(sample_ns * 1000.0 / dt_ps))
    sim.reporters.append(app.DCDReporter(str(dcd), n_report))
    sim.reporters.append(StepProgress(f"r{replica} prod", n_progress, dt_ps, n_prod, enabled=show_progress, start_step=sim.currentStep))
    sim.step(n_prod)
    sim.reporters.clear()

    n_frames = n_prod // n_report
    wall = time.time() - t0
    ReplicaProgress(enabled=show_progress).update(replica, replica, n_frames, wall)

    return {
        "replica": replica, "dcd": str(dcd), "top_pdb": str(top_pdb), "n_frames": n_frames,
        "residues_frozen": counts["residues_frozen"], "residues_mobile": counts["residues_mobile"],
        "platform": platform_used, "wall_seconds": round(wall, 1),
    }


def _replica_task(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Picklable wrapper for `run_one_replica`, for `parallel.parallel_map`."""
    return run_one_replica(**kwargs)


def sample_pocket(
    protein_pdb: Path, pocket_residues: Sequence[Residue], pocket_center_a: Sequence[float], out_dir: Path, *,
    n_replicas: int = 4, n_jobs: int = 1, temperature_k: float = 300.0, timestep_fs: float = 4.0,
    mobile_radius_nm: float = 1.0, mobile_margin_residues: int = 2,
    force_constant_kj_per_mol_nm2: float = 1.0e3, equil_ps: float = 20.0, sample_ns: float = 2.0,
    report_ps: float = 5.0, progress_ps: float = 100.0, platform_name: str = "CPU",
    nonbonded_cutoff_nm: float = 1.5, protein_forcefield: str = "amber14-all", solvent: str = "implicit",
    implicit_solvent: str = "gbn2", water_model: str = "tip3p", solvent_padding_nm: float = 1.0,
    ion_concentration_molar: float = 0.15, pressure_atm: float = 1.0, show_progress: bool = True,
) -> SampleResult:
    """Run `n_replicas` independent restrained-MD replicas (each an
    independent task, runnable via `--n-jobs`' `ProcessPoolExecutor` map)
    and write `restraint_report.json` summarizing the frozen/mobile-residue
    counts and run parameters.

    When `n_jobs != 1`, each replica's CPU platform thread count is pinned
    to `os.cpu_count() // n_workers` (the same "don't oversubscribe the
    machine" pattern `dd_docking/screening.py` uses for parallel Vina
    workers): the CPU platform's own default is normally "every logical
    core", so N concurrently-running replica processes each independently
    grabbing every core would fight each other for the same cores instead
    of actually running in parallel. `n_jobs == 1` (sequential, the
    default) leaves the thread count unpinned, so a single replica can use
    the whole machine.

    `solvent="explicit"` solvates `protein_pdb` exactly once (via
    `_solvate_and_write`, writing `out_dir/solvated_input.pdb`) before any
    replica starts, and every replica loads that shared file rather than
    solvating independently -- see `_solvate_and_write`'s docstring for why.
    """
    import json
    import os

    from .parallel import parallel_map

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    if n_jobs == 1:
        cpu_threads = 0
    else:
        n_workers = n_jobs if n_jobs and n_jobs > 0 else (os.cpu_count() or 1)
        cpu_threads = max(1, (os.cpu_count() or 1) // n_workers)

    if solvent == "explicit":
        if show_progress:
            print(
                f"[sample] solvating (protein_forcefield={protein_forcefield}, water_model={water_model}, "
                f"padding={solvent_padding_nm}nm, ionic_strength={ion_concentration_molar}M)...", flush=True,
            )
        replica_protein_pdb = _solvate_and_write(
            protein_pdb, out_dir / "solvated_input.pdb", protein_forcefield=protein_forcefield,
            water_model=water_model, padding_nm=solvent_padding_nm, ion_concentration_molar=ion_concentration_molar,
        )
        if show_progress:
            from openmm.app import PDBFile
            n_atoms = PDBFile(str(replica_protein_pdb)).topology.getNumAtoms()
            print(f"[sample] solvated: {n_atoms} atom(s) total -> {replica_protein_pdb}", flush=True)
    else:
        replica_protein_pdb = Path(protein_pdb)

    tasks = [
        {
            "protein_pdb": replica_protein_pdb, "pocket_residues": list(pocket_residues),
            "pocket_center_a": tuple(pocket_center_a), "workdir": out_dir, "replica": r,
            "temperature_k": temperature_k, "timestep_fs": timestep_fs,
            "mobile_radius_nm": mobile_radius_nm, "mobile_margin_residues": mobile_margin_residues,
            "force_constant_kj_per_mol_nm2": force_constant_kj_per_mol_nm2,
            "equil_ps": equil_ps, "sample_ns": sample_ns, "report_ps": report_ps,
            "progress_ps": progress_ps, "platform_name": platform_name,
            "nonbonded_cutoff_nm": nonbonded_cutoff_nm, "protein_forcefield": protein_forcefield,
            "solvent": solvent, "implicit_solvent": implicit_solvent, "water_model": water_model,
            "pressure_atm": pressure_atm, "cpu_threads": cpu_threads, "show_progress": show_progress,
        }
        for r in range(1, n_replicas + 1)
    ]
    results = sorted(parallel_map(_replica_task, tasks, n_jobs=n_jobs), key=lambda r: r["replica"])

    result = SampleResult(
        complex_top_pdb=results[0]["top_pdb"], replica_dcds=[r["dcd"] for r in results],
        n_replicas=n_replicas, n_frames_total=sum(r["n_frames"] for r in results),
        residues_frozen=results[0]["residues_frozen"], residues_mobile=results[0]["residues_mobile"],
        force_constant_kj_per_mol_nm2=force_constant_kj_per_mol_nm2, mobile_radius_nm=mobile_radius_nm,
        sample_ns=sample_ns, equil_ps=equil_ps, report_ps=report_ps, timestep_fs=timestep_fs,
        protein_forcefield=protein_forcefield, solvent=solvent, implicit_solvent=implicit_solvent,
        water_model=water_model, platform=results[0]["platform"], wall_seconds=round(time.time() - t0, 1),
    )
    (out_dir / "restraint_report.json").write_text(json.dumps(result.to_dict(), indent=2))
    return result
