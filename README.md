[Japanese version](README.jp.md)

# dd_afpocket — Pocket detection meets restrained-MD sampling: one AlphaFold model in, a docking-ready ensemble out

Turns a static AlphaFold DB predicted structure into a small ensemble of
receptor conformations suitable for ensemble docking: detects a druggable
pocket, samples that pocket's local flexibility with restrained MD (only
the pocket neighborhood moves -- the rest of the protein is position-
restrained), and structurally clusters the resulting trajectory into a
handful of representative conformations. Designed as a reusable package,
not tied to any specific target. Fully self-contained, with no other
runtime dependency: `prep.py` carries its own AFDB-only fetch + MD-grade
repair logic (AlphaFold models have none of the real-PDB-deposit quirks
that would otherwise call for a general-purpose structure-prep toolchain,
so only the narrow AFDB/MD-repair path was implemented), and pocket
detection, restrained-MD sampling, and clustering are all implemented
directly in this package (see "Design notes" below).

- **Fetch (`dd_afpocket-fetch`)**: UniProt accession -> AlphaFold DB model ->
  MD-grade repair, via `prep.py`'s `fetch_and_prepare_afdb`.
- **Pocket (`dd_afpocket-pocket`)**: runs `fpocket` and selects a pocket by
  Druggability Score rank (default: top-ranked), writing the pocket's
  lining residues/center (`pocket_report.json`) and a docking box derived
  from those residues' coordinates (`pocket_box.json`, since there is no
  co-crystal ligand to derive one from -- these are apo structures).
- **Sample (`dd_afpocket-sample`)**: restrained MD, implicit solvent (GBn2 by
  default, `--implicit-solvent`) or an explicit periodic water box
  (`--solvent explicit`, `--water-model`), over a choice of protein
  forcefield (`--protein-forcefield`). Residues outside the pocket
  neighborhood are harmonically restrained to their starting positions;
  only pocket-lining residues (plus anything within `--mobile-radius-nm` of
  the pocket centroid) move freely. Runs several independent replicas
  (`--n-replicas`) in parallel (`--n-jobs`).
- **Cluster (`dd_afpocket-cluster`)**: pools every replica's frames, computes a
  pocket-atom RMSD distance matrix, and clusters it (`AgglomerativeClustering`,
  average linkage) into `--n-clusters` representative structures (default
  10), each the cluster's medoid frame (never a coordinate average).
  `--pocket-expand-only` restricts clustering to frames where the pocket
  has geometrically opened up relative to the pre-MD structure;
  `--visualize` writes a standalone HTML overlay of every representative
  structure.
- **Run (`dd_afpocket-run`)**: the four stages above, end-to-end.

## Installation

Requires rdkit, numpy, pandas, pdbfixer, openmm, mdtraj, matplotlib, scipy,
scikit-learn, py3Dmol (all declared in `pyproject.toml`'s `dependencies`),
plus the `fpocket` CLI binary (for `dd_afpocket-pocket`; not available via
pip, conda-forge only). Uses its own dedicated conda env `dd_afpocket`
(python 3.12):

```bash
mamba create -n dd_afpocket -c conda-forge python=3.12 rdkit numpy pandas pdbfixer openmm mdtraj matplotlib scipy scikit-learn py3dmol pytest fpocket
/opt/miniforge3/envs/dd_afpocket/bin/pip install --no-deps -e .
```

This installs five console commands: `dd_afpocket-fetch`, `dd_afpocket-pocket`,
`dd_afpocket-sample`, `dd_afpocket-cluster`, `dd_afpocket-run`.

## Usage

### 1. Fetch (`dd_afpocket-fetch`)

```bash
dd_afpocket-fetch O60674 -o data/prepped
```

Downloads the current AlphaFold DB model for human JAK2 (UniProt O60674)
and MD-grade-repairs it (pLDDT-based terminal trimming, PDBFixer, pH 7.0
protonation) via `prep.py`, writing `data/prepped/o60674_md.pdb` (`dd_afpocket-fetch`
does not nest output under a `<uniprot>/` subdirectory -- only `dd_afpocket-run`
does, for its own multi-stage output layout below).

### 2. Pocket detection (`dd_afpocket-pocket`)

```bash
dd_afpocket-pocket data/prepped/o60674_md.pdb -o data/prepped/o60674_pocket
```

Runs fpocket and prints every detected pocket ranked by Druggability Score.
Measured output on the JAK2 model above (87 candidate pockets found across
the full multi-domain structure; abbreviated here):

```
[pocket 1] score=0.158  druggability=0.829  n_alpha_spheres=70
[pocket 2] score=0.056  druggability=0.497  n_alpha_spheres=54
[pocket 3] score=0.027  druggability=0.455  n_alpha_spheres=111
...
[done] pocket rank 1 (druggability=0.829, 17 residue(s)) -> data/prepped/o60674_pocket/pocket_report.json
```

`--pocket-rank N` selects a different rank (e.g. a known non-top-ranked
allosteric site); `--pocket-residues A:42,A:87,...` bypasses fpocket's
residue detection entirely for a manually-specified site.

`--visualize` additionally writes `pocket_candidates.html`, a standalone
(no server) py3Dmol view comparing the top `--visualize-top-n` candidate
pockets (default 3) on the receptor structure: translucent cartoon, opaque
sticks over each pocket's lining residues, that pocket's fpocket
alpha-sphere cloud rendered as a solid-colored sphere cluster (the cavity's
actual empty-space volume, not just the residues around it), and one label
per lining residue -- one color per rank, so the (often large) gap between
rank 1 and the rest reads visually instead of only off the printed
druggability scores:

```bash
dd_afpocket-pocket data/prepped/q8izl9_md.pdb -o data/prepped/q8izl9_pocket --visualize
```

Measured on the human CDK20 model (UniProt Q8IZL9): 28 candidate pockets
found, rank 1 (druggability 0.646, the known ATP-binding kinase pocket --
Gly-rich loop residues 10-18, catalytic Lys33, catalytic-loop Asp127) is
~28x more druggable than rank 2 (0.023) and every pocket beyond that is
below ~0.02, i.e. a shallow surface indentation rather than a real
druggable site -- exactly the kind of gap `pocket_candidates.html` is meant
to make visible at a glance.

### 3. Restrained-MD sampling (`dd_afpocket-sample`)

```bash
dd_afpocket-sample data/prepped/o60674_md.pdb data/prepped/o60674_pocket \
    -o data/sample/o60674 --n-replicas 4 --n-jobs 2
```

Position-restrains every residue outside the pocket neighborhood (harmonic,
`k=1000 kJ/mol/nm^2` by default -- verified to give an RMS thermal
fluctuation of ~0.86 A at 300 K, small next to the scale of pocket-shape
differences clustering looks for) and runs
`--n-replicas` independent implicit-solvent MD replicas. On the JAK2 pocket
above, 17 lining residues plus everything within 1 nm of its centroid left
49 of 1097 total residues mobile (`residues_mobile` in
`restraint_report.json`) -- see "Performance" below for realistic timing on
a target this large.

`--preset {default,quick}` bundles `--n-replicas`/`--equil-ps`/`--sample-ns`/
`--report-ps`/`--timestep-fs`/`--implicit-solvent` into one flag; `--preset
quick` is a coarse, CPU-only-friendly setting (2 replicas, 5 ps
equilibration, 300 ps production, 2.5 ps report interval, `obc2` implicit
solvent) for when only rough pocket-shape diversity is needed for downstream
ensemble docking, not a converged trajectory -- reasonable given `dd_afpocket`
produces an ensemble of conformations rather than a single stability verdict
in the first place (see "Limitations"). Any of those six flags given
explicitly still overrides the preset's value for that one flag, e.g.
`--preset quick --n-replicas 4`. `--timestep-fs` (default 4 fs, kept as-is by
both presets) is exposed separately since it is the one value
`restraints.py` empirically checked for stability against the default
`k=1000 kJ/mol/nm^2` restraint -- raising it further is not a preset-driven
trade-off.

`--implicit-solvent {gbn2,gbn,obc2,obc1,hct}` selects the GB model (default
`gbn2`). Measured on this project's development machine (human lysozyme,
P61626, ~2300 atoms, CPU/4 threads, ms/step relative to `gbn2`): `hct` 1.75x
faster, `obc1` 1.55x, `obc2` 1.37x, `gbn` ~1.02x (`gbn2` = `gbn` plus a
"neck" Born-radius correction term, which is nearly all of `gbn2`'s extra
cost over `gbn`). `--preset quick` picks `obc2`: `hct`, the fastest, is the
oldest GB model and is known to underestimate buried atoms' Born radii,
which matters here since some pocket-lining side chains are partially
buried; `obc2` fixes that at a smaller speed cost and is the de facto
standard GB model in AMBER-family tooling (`igb=5`).

### Forcefield and solvent model

`--protein-forcefield {amber14-all,amber99sbildn,amber19-all,charmm36}`
(default `amber14-all`) selects the protein XML(s); every combination with
every `--implicit-solvent` model is verified to build. `--solvent
{implicit,explicit}` (default `implicit`) switches between GB continuum
solvent (the setting this project's other CPU-friendly defaults assume) and
a real periodic water box. `--solvent explicit` needs `--water-model`,
verified per `--protein-forcefield` (`sample.PROTEIN_FORCEFIELDS`):

| `--protein-forcefield` | verified `--water-model` choices |
|---|---|
| `amber14-all` | `tip3p`, `tip3pfb`, `tip4pew`, `tip4pfb`, `spce` |
| `amber99sbildn` | `tip3p`, `tip4pew`, `tip5p` |
| `amber19-all` | `tip3p`, `tip3pfb`, `tip4pew` |
| `charmm36` | `tip3p`, `tip4pew`, `tip5p` |

Explicit solvent also adds `--solvent-padding-nm` (solute-to-box-edge
padding, default 1.0), `--ion-concentration-molar` (Na+/Cl- neutralizing/
background concentration, default 0.15, physiological), and
`--pressure-atm` (`MonteCarloBarostat` target, default 1.0 -- explicit runs
are NPT, not NVT). Solvation (`Modeller.addSolvent`, which has no seed
parameter) happens exactly once per `dd_afpocket-sample` invocation, before any
replica starts, and every replica loads that one shared, already-periodic
structure -- re-solvating independently per replica would give each
replica's DCD a different atom count/ordering than the `complex_top.pdb`
`dd_afpocket-cluster` loads them all against.

**Explicit solvent is dramatically slower than implicit on CPU** -- measured
end-to-end on this project's development machine: human lysozyme solvated
with ~0.6 nm padding (amber14-all + tip3p, 105,716 atoms total including
water/ions, CPU, thread count unpinned) ran a full 1-replica correctness
check (minimize -> 0.08 ps equilibration -> 20 ps production, `--n-jobs 1`)
in 2,304 s wall, with production itself holding steady at ~87 s/ps once
past the first couple of ps -- versus the same protein's own implicit-GBn2
number in "Performance" below (~2,300 atoms, CPU/4 threads pinned,
~0.09 s/step at a 4 fs timestep, i.e. ~22.7 s/ps). Already several-fold
slower per ps even though this isn't a clean apples-to-apples comparison
(the explicit system's 105,716 atoms is ~46x the implicit one's atom
count, mostly water, and thread pinning differed between the two runs).
The real-world takeaway is unambiguous either way: `--platform CUDA` is
strongly recommended for any `--solvent explicit` run beyond a quick
correctness check; `--solvent implicit` (the default) remains the
practical choice for CPU-only environments (see "Performance" below).
`dd_afpocket-cluster` loading a `>99,999`-atom explicit-solvent `complex_top.pdb`
also prints a `mdtraj... Need to guess atom number ...` warning -- benign,
just the fixed-width PDB atom-serial field overflowing past 99,999
(mdtraj tracks atoms by index regardless, not by that text field).

### 4. Clustering (`dd_afpocket-cluster`)

```bash
dd_afpocket-cluster data/sample/o60674 data/prepped/o60674_pocket -o data/clusters/o60674 --n-clusters 10
```

Writes `cluster_00.pdb` (largest cluster) through `cluster_09.pdb`, plus
`cluster_report.csv` (population, source replica/time, mean intra-cluster
RMSD per representative, and each representative's `pocket_volume_nm3` --
see below).

#### Post-hoc pocket-expansion filtering (`--pocket-expand-only`)

`cluster.pocket_volume_proxy` computes a cheap, no-rerun geometric proxy for
"how open is the pocket": the convex-hull volume of the pocket-lining atoms
(`--cluster-atoms`) in each pooled frame -- not fpocket's own cavity volume,
which would mean re-running fpocket per frame, far too slow (see
"Performance"). `--pocket-expand-only` restricts clustering to frames whose
volume is `>= reference * (1 + --pocket-expand-margin)`, where `reference`
is that same proxy computed on the pre-MD structure -- a post-hoc selection
over an already-sampled, otherwise-unbiased ensemble (`dd_afpocket-sample` itself
knows nothing about "expansion"; nothing about how the frames were
generated changes, only which of them are eligible to become a
representative structure).

**Calibrate `--pocket-expand-margin` -- the default (0.0) barely filters
anything.** A convex hull is an envelope over its points, so thermal noise
alone tends to inflate it relative to any single static reference frame,
regardless of any real "opening" signal -- measured on a real (if short)
sampling run, margin 0.0 kept 10,000 of 10,000 pooled frames. `--pocket-
expand-only` always prints `[cluster] pocket-expand-only: kept X/Y
frame(s) ...`; start `--pocket-expand-margin` around 0.05-0.1 and adjust
from that printed ratio rather than trusting the default to do anything.

#### Visual comparison (`--visualize`)

Writes `cluster_overlay.html` in `--out-dir`: a standalone (no server)
py3Dmol scene with every `cluster_NN.pdb` loaded as its own model, each
given a distinct color (translucent cartoon for the whole receptor, opaque
sticks over the pocket-lining residues, since that's the part that actually
varies between clusters -- everything else was frozen by `restraints.py`
during sampling). No structural alignment step is needed first: every
cluster PDB already shares dd_afpocket's restrained-MD reference frame.  Opening
the file needs internet access (3Dmol.js loads from a CDN, same as
py3Dmol's normal Jupyter usage).

### End-to-end (`dd_afpocket-run`)

```bash
dd_afpocket-run O60674 -o data --n-replicas 4 --n-jobs 2 --n-clusters 10
```

Runs all four stages, writing everything under `data/o60674/`.

### Feeding the ensemble into `dd_docking`

The N representative structures are protein-only (apo) PDBs, directly
usable as `dd_docking`'s ensemble-member input. Since there is no
co-crystal ligand to derive a docking box from, use `pocket_box.json`'s
box (center/size) for every member -- they share the same pocket
definition and residue frame, so the box is consistent across all of them.

## Design notes

- **Plain `openmm.app.ForceField`, not `SystemGenerator`.** `SystemGenerator`
  exists for systems that include a small-molecule ligand (GAFF/SMIRNOFF
  parameterization). `dd_afpocket`'s systems are apo -- never a ligand -- so a
  plain `ForceField(...)` is all that's needed. This also sidesteps a real
  environment issue: constructing a `SystemGenerator` without an explicit
  `small_molecule_forcefield="gaff-*"` eagerly triggers `openff.toolkit`'s
  SMIRNOFF force-field discovery, which can fail with
  `ModuleNotFoundError: No module named 'pkg_resources'`
  (`openff-amber-ff-ports` depends on the `pkg_resources` API that
  setuptools >= 81 no longer ships) -- avoided entirely by not needing
  `SystemGenerator` in the first place.
- **A curated, empirically-verified forcefield/water registry
  (`sample.PROTEIN_FORCEFIELDS`), not "every XML OpenMM ships".** OpenMM
  bundles far more protein/water/GB parameter files than dd_afpocket exposes; an
  untested pairing (e.g. an AMBER-parametrized GB file against CHARMM atom
  types) can fail outright or silently build a mismatched system. Every
  `--protein-forcefield` x `--implicit-solvent` combination and every
  `--protein-forcefield` x `--water-model` combination actually listed was
  confirmed (via `ForceField(...).createSystem(...)` / `Modeller.
  addSolvent(...)` + `createSystem(...)`) to build successfully on this
  project's development machine before being added.
- **Explicit solvent is solvated exactly once per `dd_afpocket-sample`
  invocation, shared by every replica.** `Modeller.addSolvent` has no seed
  parameter -- its ion placement is randomized -- so solvating
  independently per replica (safe for the *implicit* system build, which is
  deterministic) would give each replica's DCD a different atom count/
  ordering than whatever `complex_top.pdb` `dd_afpocket-cluster` ends up loading
  every replica's DCD against. `sample.sample_pocket` solvates once,
  writes the periodic structure to `solvated_input.pdb`, and every replica
  loads that file as plain, deterministic input -- structurally identical
  to how the apo (implicit) case already worked.
- **`rigidWater=True` and a `MonteCarloBarostat` for explicit solvent
  only.** Implicit systems have no water to make rigid (the flag has no
  effect); explicit water is conventionally kept rigid (bond/angle
  constraints, not hydrogen-mass-repartitioned) rather than integrated at
  the protein's own 4 fs HMR timestep. The barostat runs NPT throughout
  equilibration and production at `--pressure-atm` (default 1 atm) --
  `restraints.py`'s harmonic (not literally frozen) restraints don't
  prevent the barostat's periodic volume-rescaling trials, they just pull
  restrained atoms back afterward, same as any other perturbation to a
  restrained atom's position.
- **`CutoffNonPeriodic`, not `NoCutoff`.** GBn2's Born-radius/GB-energy
  terms already approximate the far-field electrostatic effect of distant
  atoms, but the pairwise nonbonded loop itself is a naive O(n_atoms^2)
  evaluation regardless of that approximation. For a large multi-domain
  protein (full-length JAK2: 1132 residues, ~17800 atoms with hydrogens),
  `NoCutoff` measured impractically slow even for a handful of picoseconds
  on CPU. Since `restraints.py` freezes everything outside a ~1 nm pocket
  neighborhood anyway, a 1.5 nm nonbonded cutoff (default,
  `--nonbonded-cutoff-nm`) costs essentially no accuracy for the region
  this project actually samples.
- **CPU thread pinning across parallel replicas.** The CPU platform's
  default thread count is normally "every logical core". Running several
  replicas concurrently via `--n-jobs` without limiting each one's thread
  count means N processes fight each other for the same cores instead of
  actually running in parallel -- `sample_pocket` pins each replica's
  thread count to `os.cpu_count() // n_workers` when `n_jobs != 1`
  (`n_jobs == 1`, the default, leaves the thread count unpinned so a
  single replica can use the whole machine).

## Performance

Implicit-solvent (GBn2) CPU MD is not fast, and force-evaluation cost
scales with total atom count even though only the pocket neighborhood is
actually integrated (the restrained region still needs its forces
computed every step). Measured end-to-end (`dd_afpocket-sample`, 2 replicas via
`--n-jobs 2`, ~8 CPU threads pinned per replica) on this project's
development machine (a shared 16-core desktop, not a dedicated compute
node):

- Streptavidin (1STP, 121 residues, 1744 atoms with hydrogens): 5250
  total steps (equilibration + production) per replica in 384 s wall --
  ~0.073 s/step.
- Human lysozyme (P61626, 146 residues, ~2300 atoms): 2100 steps per
  replica in 232 s wall -- ~0.110 s/step.
- Full-length JAK2 (O60674, 1132 residues, ~17800 atoms) is substantially
  slower per step at the same thread count, purely from the larger
  force-evaluation cost, even though `restraints.py` only lets ~50
  residues actually move -- a short correctness-check run (525 steps) did
  not finish within several CPU-minutes per thread when contending with
  other jobs for cores on this shared machine. Budget accordingly (or use
  a GPU) rather than assuming the small-protein numbers above scale
  linearly with atom count; pocket detection (`dd_afpocket-pocket`, fpocket) is
  unaffected by protein size in practice -- it completed on JAK2 (87
  candidate pockets, top-ranked druggability 0.829) in a few seconds.

**Recommendations for practical turnaround:**
- Prefer a CUDA GPU (`--platform CUDA`) for anything beyond a quick
  correctness check -- implicit-solvent GB kernels are dramatically faster
  on GPU than CPU.
- For large multi-domain targets like full-length JAK2, consider providing
  a smaller construct (e.g. an isolated kinase domain) if your workflow
  allows it -- `dd_afpocket` doesn't currently offer residue-range slicing, only
  whole-chain selection via `prep.py`.
- Scale `--sample-ns`/`--equil-ps`/`--n-replicas` down for exploratory runs
  and back up once you've confirmed the pipeline behaves as expected on
  your hardware; there's no dedicated "screen then confirm" gate here (this
  project makes an ensemble of conformations, not a single stability
  verdict, so there's no natural place to cut a bad run short). `--preset quick` (see "Restrained-
  MD sampling" above) is this scaled-down bundle -- including a cheaper GB
  model (`--implicit-solvent obc2`) -- pre-picked for a CPU-only machine;
  pair it with `--n-jobs -1` to also use every core.

## Limitations

- Being apo, `dd_afpocket` cannot reproduce ligand-induced conformational change
  (induced fit) -- only the unbound protein's own local flexibility around
  a computationally-detected pocket.
- Pocket detection quality depends entirely on fpocket; a pocket it misses
  or mis-scores cannot be recovered except via `--pocket-residues`'
  manual override.
- No residue-range slicing: `dd_afpocket` samples whatever chain(s) `prep.py`
  fetched, which for a multi-domain protein may be far larger (and slower)
  than the domain actually relevant to the pocket of interest.

## Module layout

| Module | Purpose |
|---|---|
| `pocket.py` | fpocket subprocess wrapper, Druggability Score ranking, pocket-residue/center extraction, docking-box computation |
| `restraints.py` | Pocket-center-based harmonic position restraints (mobile-residue-set computation + `CustomExternalForce`) |
| `sample.py` | Forcefield/water-model registry, implicit- and explicit-solvent system building, multi-replica restrained-MD sampling |
| `cluster.py` | Pooled-trajectory RMSD clustering, medoid selection, pocket-volume proxy/filtering, representative-structure/report output |
| `visualize.py` | Standalone py3Dmol HTML: cluster-representative-structure overlay, and top-candidate-pocket comparison (lining residues, alpha-sphere cavity volume, residue labels) |
| `pipeline.py` | Per-stage orchestration functions plus `run_end_to_end` |
| `cli.py` | `dd_afpocket-fetch`/`dd_afpocket-pocket`/`dd_afpocket-sample`/`dd_afpocket-cluster`/`dd_afpocket-run` argparse entry points |
| `progress.py` | Progress-line printing (per-replica, per-pocket, per-cluster, in-run OpenMM step reporter) |
| `parallel.py` | `ProcessPoolExecutor`-based parallel map |
