"""Progress reporting: one concise line per completed item (mirrors
`dd_docking/progress.py`'s convention), plus a periodic in-loop reporter for
the sampling MD itself (`StepProgress`, mirroring `dd_mdstability/progress.py` --
even a CPU-friendly implicit-solvent run is a slow loop that should never go
silent between start and finish).
"""
import time


class StepProgress:
    """OpenMM reporter protocol: prints one line every `interval_steps`
    steps during equilibration/production, so long runs stay visible
    (simulated time, wall-clock elapsed, potential energy, temperature)."""

    def __init__(self, stage: str, interval_steps: int, dt_ps: float, total_steps: int,
                 enabled: bool = True, start_step: int = 0):
        self.stage = stage
        self.interval_steps = max(1, interval_steps)
        self.dt_ps = dt_ps
        self.total_steps = total_steps
        self.enabled = enabled
        self.start_step = start_step
        self._t0 = time.time()

    def describeNextReport(self, simulation):
        steps = self.interval_steps - simulation.currentStep % self.interval_steps
        return {"steps": steps, "periodic": None, "include": ["energy"]}

    def report(self, simulation, state):
        if not self.enabled:
            return
        from openmm import unit

        step = simulation.currentStep - self.start_step
        e_pot = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        # Degrees of freedom, not 3*N: HBonds constraints remove one DOF
        # each, and OpenMM's default CMMotionRemover removes 3 more.
        system = simulation.system
        dof = 3 * system.getNumParticles() - system.getNumConstraints() - 3
        temp = (2 * state.getKineticEnergy() / (dof * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)
        elapsed = time.time() - self._t0
        print(
            f"[sample {self.stage}] {step * self.dt_ps:.1f}/{self.total_steps * self.dt_ps:.1f} ps  "
            f"Epot={e_pot:.1f} kcal/mol  T={temp:.1f} K  elapsed={elapsed:.0f}s",
            flush=True,
        )


class ReplicaProgress:
    """Prints one line per finished sampling replica:
    `[sample] replica <i>/<n>  n_frames=.. wall=..s`"""

    def __init__(self, enabled: bool = True):
        self.enabled = enabled

    def update(self, replica: int, n_replicas: int, n_frames: int, wall_seconds: float) -> None:
        if not self.enabled:
            return
        print(
            f"[sample] replica {replica}/{n_replicas}  n_frames={n_frames}  wall={wall_seconds:.0f}s",
            flush=True,
        )


class PocketProgress:
    """Prints one line per detected fpocket pocket:
    `<rank>  score=.. druggability=.. n_alpha_spheres=..`"""

    def __init__(self, enabled: bool = True):
        self.enabled = enabled

    def update(self, rank: int, score: float, druggability: float, n_alpha_spheres: int) -> None:
        if not self.enabled:
            return
        print(
            f"[pocket {rank}] score={score:.3f}  druggability={druggability:.3f}  "
            f"n_alpha_spheres={n_alpha_spheres}",
            flush=True,
        )


class ClusterProgress:
    """Prints one line per cluster assigned to a representative structure:
    `cluster <id>  n_frames=.. mean_rmsd_A=.. -> <path>`"""

    def __init__(self, enabled: bool = True):
        self.enabled = enabled

    def update(self, cluster_id: int, n_frames: int, mean_rmsd_a: float, out_path: str) -> None:
        if not self.enabled:
            return
        print(
            f"[cluster {cluster_id}] n_frames={n_frames}  mean_rmsd_A={mean_rmsd_a:.2f}  -> {out_path}",
            flush=True,
        )
