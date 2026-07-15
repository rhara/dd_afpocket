# dd_af

Turns a static AlphaFold DB predicted structure into a small ensemble of
receptor conformations suitable for ensemble docking: detects a druggable
pocket, samples that pocket's local flexibility with restrained MD (only
the pocket neighborhood moves -- the rest of the protein is position-
restrained), and structurally clusters the resulting trajectory into a
handful of representative conformations. Designed as a reusable package,
not tied to any specific target (same philosophy as `dd_prep` / `dd_docking`
/ `dd_overlay` / `dd_viewer` / `dd_confgen`) -- conceptually the middle
stage of a `dd_prep` (AFDB fetch + MD-grade repair) -> **`dd_af`** (pocket
detection + local restrained-MD sampling + clustering) -> `dd_docking`
(ensemble docking against the generated conformations) pipeline. `dd_af`
imports `dd_prep` directly for the fetch/repair step (AlphaFold models have
none of the real-PDB-deposit quirks `dd_md`'s self-contained receptor prep
exists to handle, so there was no reason to reimplement it here); it does
not import `dd_docking`/`dd_md`, only mirrors their docking-box convention
and harmonic-restraint mechanics respectively (see "Design notes" below).

- **Fetch (`dd_af-fetch`)**: UniProt accession -> AlphaFold DB model ->
  MD-grade repair, delegating to `dd_prep.pipeline.fetch_and_prepare_afdb`.
- **Pocket (`dd_af-pocket`)**: runs `fpocket` and selects a pocket by
  Druggability Score rank (default: top-ranked), writing the pocket's
  lining residues/center (`pocket_report.json`) and a docking box derived
  from those residues' coordinates (`pocket_box.json`, since there is no
  co-crystal ligand to derive one from -- these are apo structures).
- **Sample (`dd_af-sample`)**: implicit-solvent (GBn2) restrained MD.
  Residues outside the pocket neighborhood are harmonically restrained to
  their starting positions; only pocket-lining residues (plus anything
  within `--mobile-radius-nm` of the pocket centroid) move freely. Runs
  several independent replicas (`--n-replicas`) in parallel (`--n-jobs`).
- **Cluster (`dd_af-cluster`)**: pools every replica's frames, computes a
  pocket-atom RMSD distance matrix, and clusters it (`AgglomerativeClustering`,
  average linkage) into `--n-clusters` representative structures (default
  10), each the cluster's medoid frame (never a coordinate average).
- **Run (`dd_af-run`)**: the four stages above, end-to-end.

## Installation

Requires rdkit, numpy, pandas, pdbfixer, openmm, openmmforcefields, mdtraj,
scikit-learn, and the `dd_prep` package (for `dd_af-fetch`), plus the
`fpocket` CLI binary (for `dd_af-pocket`). Best installed via conda-forge
(the `mpro` env already has everything except `fpocket`):

```bash
conda install -n mpro -c conda-forge fpocket
cd dd_prep && pip install -e . && cd ..   # if not already installed
cd dd_af && pip install -e .
```

This installs five console commands: `dd_af-fetch`, `dd_af-pocket`,
`dd_af-sample`, `dd_af-cluster`, `dd_af-run`.

## Usage

### 1. Fetch (`dd_af-fetch`)

```bash
dd_af-fetch O60674 -o data/prepped
```

Downloads the current AlphaFold DB model for human JAK2 (UniProt O60674)
and MD-grade-repairs it (pLDDT-based terminal trimming, PDBFixer, pH 7.0
protonation) via `dd_prep`, writing `data/prepped/o60674_md.pdb` (`dd_af-fetch`
does not nest output under a `<uniprot>/` subdirectory -- only `dd_af-run`
does, for its own multi-stage output layout below).

### 2. Pocket detection (`dd_af-pocket`)

```bash
dd_af-pocket data/prepped/o60674_md.pdb -o data/prepped/o60674_pocket
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

### 3. Restrained-MD sampling (`dd_af-sample`)

```bash
dd_af-sample data/prepped/o60674_md.pdb data/prepped/o60674_pocket \
    -o data/sample/o60674 --n-replicas 4 --n-jobs 2
```

Position-restrains every residue outside the pocket neighborhood (harmonic,
`k=1000 kJ/mol/nm^2` by default -- the same value `dd_md/restraints.py`
verified gives an RMS thermal fluctuation of ~0.86 A at 300 K, small next
to the scale of pocket-shape differences clustering looks for) and runs
`--n-replicas` independent implicit-solvent MD replicas. On the JAK2 pocket
above, 17 lining residues plus everything within 1 nm of its centroid left
49 of 1097 total residues mobile (`residues_mobile` in
`restraint_report.json`) -- see "Performance" below for realistic timing on
a target this large.

### 4. Clustering (`dd_af-cluster`)

```bash
dd_af-cluster data/sample/o60674 data/prepped/o60674_pocket -o data/clusters/o60674 --n-clusters 10
```

Writes `cluster_00.pdb` (largest cluster) through `cluster_09.pdb`, plus
`cluster_report.csv` (population, source replica/time, mean intra-cluster
RMSD per representative).

### End-to-end (`dd_af-run`)

```bash
dd_af-run O60674 -o data --n-replicas 4 --n-jobs 2 --n-clusters 10
```

Runs all four stages, writing everything under `data/o60674/`.

### Feeding the ensemble into `dd_docking`

The N representative structures are protein-only (apo) PDBs, directly
usable as `dd_docking`'s ensemble-member input. Since there is no
co-crystal ligand to derive a docking box from, use `pocket_box.json`'s
box (center/size) for every member -- they share the same pocket
definition and residue frame, so the box is consistent across all of them.

## Design notes

- **Plain `openmm.app.ForceField`, not `SystemGenerator`.** `dd_docking/
  refine_md.py` and `dd_md/system_build.py` need `SystemGenerator` because
  their systems include a small-molecule ligand (GAFF/SMIRNOFF
  parameterization). `dd_af`'s systems are apo -- never a ligand -- so a
  plain `ForceField("amber14-all.xml", "implicit/gbn2.xml")` is all that's
  needed. This also sidesteps a real environment issue: constructing a
  `SystemGenerator` without an explicit `small_molecule_forcefield="gaff-*"`
  eagerly triggers `openff.toolkit`'s SMIRNOFF force-field discovery,
  which fails in the `mpro` env with `ModuleNotFoundError: No module named
  'pkg_resources'` (`openff-amber-ff-ports` depends on the `pkg_resources`
  API that setuptools >= 81 no longer ships). Confirmed this does not
  affect `dd_docking`/`dd_md`, since both always pass
  `small_molecule_forcefield="gaff-2.11"` explicitly, a different
  (GAFF) code path that never reaches the broken import.
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
  thread count to `os.cpu_count() // n_workers` when `n_jobs != 1`, the
  same pattern `dd_docking/screening.py` uses for parallel Vina workers
  (`n_jobs == 1`, the default, leaves the thread count unpinned so a
  single replica can use the whole machine).

## Performance

Implicit-solvent (GBn2) CPU MD is not fast, and force-evaluation cost
scales with total atom count even though only the pocket neighborhood is
actually integrated (the restrained region still needs its forces
computed every step). Measured end-to-end (`dd_af-sample`, 2 replicas via
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
  linearly with atom count; pocket detection (`dd_af-pocket`, fpocket) is
  unaffected by protein size in practice -- it completed on JAK2 (87
  candidate pockets, top-ranked druggability 0.829) in a few seconds.

**Recommendations for practical turnaround:**
- Prefer a CUDA GPU (`--platform CUDA`) for anything beyond a quick
  correctness check -- implicit-solvent GB kernels are dramatically faster
  on GPU than CPU.
- For large multi-domain targets like full-length JAK2, consider providing
  a smaller construct (e.g. an isolated kinase domain) if your workflow
  allows it -- `dd_af` doesn't currently offer residue-range slicing, only
  whole-chain selection via `dd_prep`.
- Scale `--sample-ns`/`--equil-ps`/`--n-replicas` down for exploratory runs
  and back up once you've confirmed the pipeline behaves as expected on
  your hardware; there's no dedicated "screen then confirm" gate here (this
  project makes an ensemble of conformations, not a single stability
  verdict, so there's no natural place to cut a bad run short the way
  `dd_md`'s screen-then-confirm flow does).

## Limitations

- Being apo, `dd_af` cannot reproduce ligand-induced conformational change
  (induced fit) -- only the unbound protein's own local flexibility around
  a computationally-detected pocket. Combine with `dd_docking`'s own
  induced-fit refinement (`dd_docking-refine`) downstream if that matters
  for your target.
- Pocket detection quality depends entirely on fpocket; a pocket it misses
  or mis-scores cannot be recovered except via `--pocket-residues`'
  manual override.
- No residue-range slicing: `dd_af` samples whatever chain(s) `dd_prep`
  fetched, which for a multi-domain protein may be far larger (and slower)
  than the domain actually relevant to the pocket of interest.

## Module layout

| Module | Purpose |
|---|---|
| `pocket.py` | fpocket subprocess wrapper, Druggability Score ranking, pocket-residue/center extraction, docking-box computation |
| `restraints.py` | Pocket-center-based harmonic position restraints (mobile-residue-set computation + `CustomExternalForce`) |
| `sample.py` | Implicit-solvent system building and multi-replica restrained-MD sampling |
| `cluster.py` | Pooled-trajectory RMSD clustering, medoid selection, representative-structure/report output |
| `pipeline.py` | Per-stage orchestration functions plus `run_end_to_end` |
| `cli.py` | `dd_af-fetch`/`dd_af-pocket`/`dd_af-sample`/`dd_af-cluster`/`dd_af-run` argparse entry points |
| `progress.py` | Progress-line printing (per-replica, per-pocket, per-cluster, in-run OpenMM step reporter) |
| `parallel.py` | `ProcessPoolExecutor`-based parallel map, copied from `dd_docking/parallel.py` |
