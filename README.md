# rxflow

A small transition-metal reaction workflow using `ligands.csv`, `rxembed` and `rxts`.

Clone rxflow and install it into `.venv` (Python 3.13):

```sh
git clone https://github.com/HenryTeahan/rxflow.git
cd rxflow
uv sync --locked
source .venv/bin/activate
rxflow --help
```

`uv sync --locked` installs rxflow and all Python dependencies, including GitHub
`rxembed`, `rxts` and `xyzrender`; no sibling checkouts are needed. `rxembed` is pinned to
[PR #28](https://github.com/aligfellow/rxembed/pull/28), commit
`2c60b21f0eb0aebb5361bae2a968050b616c4dd9`. `rxts` uses
[HenryTeahan/rxts](https://github.com/HenryTeahan/rxts), pinned to
`8eac9f1b2a2d13e52b47a1438961e0fe68a45f47`, containing molecule-input and GFN-FF
graph-topology fixes, without rotation search. Git dependencies use HTTPS so
installation needs no SSH credentials. Use the project install above so uv
applies the same rxembed pin to rxts's transitive dependency. GFN-FF's source
build requires a Fortran compiler and BLAS/LAPACK on the cluster.

```sh
source .venv/bin/activate
jupyter lab manual_workflow.ipynb
```

Select **Python (rxflow)**. The notebook checks its interpreter, prints the installed
coembedding source and revision, and uses one ligand. Screening results are reused
when present. Minimum optimization uses three conformers;
coembedding requests **10 pairs and selects the lowest TS RMSD pair** for one TS guess.

All outputs use CSV ligand IDs, padded to at least two digits:

```text
results/
  ligand_01/
    reaction_pairs.csv
    coordination_modes.json
    input.json
    failures.json
    coordination/
      mode_01_isomer_01_input.xyz
      mode_01_isomer_01_stdout.log
    reactant/
      pair_001_seating_03_reactive.xyz
      pair_001_seating_03_conformer_001_input.xyz
      pair_001_seating_03_conformer_001_xtbopt.xyz
      pair_001_seating_03_minimum_001.xyz
      pair_001_seating_03_minima.json
    product/
      ...
    ts_guess/
      pair_001_seating_03.xyz
      pair_001_seating_03.inp
      pair_001_seating_03_path.xyz
      pair_001_seating_03_selection.json
  ligand_02/
    ...
```

Pair, seating and conformer identifiers are in filenames. Molecular files and
calculation logs are directly inside their role folder. In the notebook,
`RUN_NAME` adds a filename prefix so another manual run can coexist.

Run the full CSV in parallel (one worker per ligand, one CPU thread each):

```sh
XTB_EXE=/path/to/xtb python run_workflow.py
```

This uses 3 conformers per endpoint optimization and 10 coembedded pairs per
reaction. Progress is written to `results/run.json` and each ligand's `status.json`.

The main functions are:

- `generate_reaction_pairs`: enumerate donor subsets and coordination isomers,
  score one conformer per isomer with gas-phase g-xTB, retain within 20 kcal/mol,
  and generate mapped H-transfer pairs over substrate sites and resonance forms.
- `create_endpoint_isomers`: retain the screened ligand seating while adding the substrate.
- `optimize_minima`: optimize endpoint conformers and retain valid minima within
  3 kcal/mol. Reject a graph change only if the associated atom-pair distance
  changes by **more than 2 Å from its starting value**.
- `generate_ts_guess`: request 10 coembedded endpoint pairs, select the minimum
  `ts_rmsds` conformer ID, then run one GFN-FF `rxts` scan and write a TS input.
- `sample_reaction_endpoints`: run minimum optimization and reactive sampling for one pair.

```sh
# Screen ligands into results/ligand_01/, results/ligand_02/, ...
XTB_EXE=/path/to/xtb python rxflow.py --max-ligands 1

# Sample minima and generate TS guesses in the same ligand folder.
XTB_EXE=/path/to/xtb python rxflow_endpoints.py results/ligand_01/reaction_pairs.csv

# Generate only the TS guess for the first pair; 10 coembedded conformers by default.
python rxflow_endpoints.py results/ligand_01/reaction_pairs.csv --stage ts --max-pairs 1
```

`--output-dir` defaults to `results`. `--conformers` controls minimum sampling
(default 3); `--coembed-conformers` controls reactive sampling (default 10).
Embedding may return fewer valid pairs. Selection JSON records all TS RMSDs and
the selected ID. Existing calculation names are protected against overwriting.
`--enumerate-only` skips embedding and scoring; those graph candidates cannot
supply endpoint sampling until their coordination isomers have been screened.

The CSV requires a `substrate` SMILES column and accepts `central_metal` or `metal`.
Each ligand row can specify its own substrate. The supplied CSV uses the original
trimethoxybenzene substrate. All symmetry-distinct aromatic **C–H** positions are
enumerated, with one representative per equivalent group; alkyl C–H positions
are excluded. `site_001` etc. identify the bound carbon's one-based atom index in
the parsed substrate (before adding hydrogens). Input atom-map labels do not
affect symmetry detection. Use a fresh results directory after changing substrate.

Charges come from the graph (+1 for
these Pd complexes). Atom maps are one-based; atom indices and TS indices are
zero-based. The notebook's broad heteroatom donor screen is retained.

Local calculations omit solvent. `run_gxtb(..., solvent="h2o")` remains available.
TS inputs default to gas-phase `r2SCAN-3c`, configurable through `--orca-method`.
They include `Calc_Hess true`, `NumHess true`, a reaction-specific `TS_Mode` and
`NumFreq` for the frequency calculation. See the [ORCA TS manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations_TS.html)
and [frequency manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/frequencies.html).
Inputs are written without submitting ORCA jobs. FF guesses require saddle-point verification.

The original `reaction_workflow.py` is a reference requiring the absent
`workflow_store` module. The local runner uses CSV/JSON files.

## HPC / ORCA

The installed `rxflow` command provides separate HPC sampling and ORCA stages.
Set `paths.ligands_csv`, `paths.source` (default `hpc_seeds`) and `paths.results`
(default `hpc_results`) in `hpc.toml`. Use fresh directories to start a new run.
Sampling begins with the CSV, screens coordination modes with gas-phase g-xTB
and the 20 kcal/mol cutoff, embeds three initial conformers per endpoint, and
creates one rxts TS guess from the best of ten coembedded pairs. It does not run
gas-phase endpoint optimization. Those endpoints are optimized by ORCA later.

Supply all external executable paths yourself. `ORCA_EXTERNAL_TOOLS` is the full
path to the **oet_gxtb_alpb executable**, wherever it is installed (including opt).
No installation directory is assumed and no separate wrapper installation in
rxflow's environment is required when your executable already works.

```sh
export ORCA_EXE=/your/path/to/orca
export XTB_EXE=/your/path/to/gxtb-enabled/xtb
export ORCA_EXTERNAL_TOOLS=/your/path/to/oet_gxtb_alpb
```

Load your cluster's ORCA/MPI modules before submitting, then start from the CSV:

```sh
rxflow submit --account YOUR_ACCOUNT --partition YOUR_PARTITION
```

Omit account/partition options when your cluster supplies defaults. `submit`
determines the sampling array size from the CSV and runs one CPU per ligand.
Database preparation waits for sampling to finish (`afterany`), including failed
tasks; the ORCA worker array waits for successful preparation (`afterok`).
It prints all three Slurm job IDs and exports the current environment to every job.
Existing runs are protected: choose fresh source/results paths to start again.

To resume pending jobs in an existing database, including jobs reset with `retry`,
submit only the ORCA worker array:

```sh
rxflow resume --workers 100 --account YOUR_ACCOUNT --partition YOUR_PARTITION
```

Omit account/partition when supplied by cluster defaults. Use `--config other.toml`
for another database. This skips sampling and preparation and leaves succeeded,
failed and running rows unchanged; workers claim only pending rows. Failed rows
must first be reset with `rxflow retry --id ID`. Rows left running by killed workers
also require an explicit retry after the old workers have stopped.

Sampling coembeds 10 pairs and creates TS guesses from the five lowest TS RMSDs
(`sampling.ts_candidates = 5`), or fewer if fewer pairs survive. Each candidate
gets a separate TS row: `conformer=1..5` is its RMSD rank, `ts_rmsd` records the
score, and `ts_strategy='direct'` identifies a direct OptTS job. A separate
`conformer=0`, `ts_strategy='neb'` row is the conditional NEB fallback.
Existing single-guess seeds retain their previous behavior; retrying old rows
does not generate additional candidates. Multiple candidates require new sampling.
Geometry-optimization cycle limits use ORCA's defaults; inputs do not set `MaxIter`.

The default is 100 independent ORCA workers (`--workers 100`), one CPU and 2 GB
RAM per worker as configured in `hpc.toml`, with a 24-hour Slurm time limit per
stage. Endpoint embedding and TS-guess failures are recorded in each ligand's
`sampling.json` and as database rows with `status='failed'`, `stage='sampling'`,
an error message and empty energy fields. Other conformers and reaction groups
continue. Inspect `sampling-*.out` for failures before a reaction group was
created; those tasks also do not hold up the other ligands.
ORCA logs are `orca-*.out` beside the config file.
Scratch uses Python's temporary directory (`TMPDIR`); set it to node-local storage
as required by your cluster. Use the cluster's shared filesystem for the database.

`hpc_sample.sbatch` and `hpc_worker.sbatch` remain available for manual submission;
the sample template's `0-12` range matches the supplied 13-row CSV. Those templates
can load cluster-specific modules when that is needed inside the batch job.

The stages are also available individually, for a compute-node allocation:

```sh
rxflow sample --ligand-index 0   # Zero-based CSV row; omit to sample all sequentially.
rxflow prepare
rxflow worker
rxflow status
```

To reuse previously sampled local results instead, set `paths.source = "results"`
and run `rxflow prepare`, then submit only the ORCA worker array. Preparation
imports initial endpoint XYZs and selected TS guesses, never local energies or
old ORCA inputs. It is repeatable without resetting existing jobs. After
preparation, only the full `hpc_results/` directory is needed, not the seed source.
Calculation settings are frozen per row during preparation; use a new results
directory when changing the protocol. The local runner remains separate.

Each worker atomically claims a row in `hpc_results/workflow.sqlite`. Its single
`jobs` table holds one row per endpoint conformer and per TS candidate/fallback per ligand/pair/seating,
including status, stage, attempts, geometries, electronic/Gibbs energies in Eh,
frequencies and errors. Paths are relative to `hpc_results/`; XYZ, INP, OUT and Hessian
files stay directly in `ligand_XX/{reactant,product,ts_guess}/`. Use shared storage
that supports SQLite file locking; the database uses rollback journaling rather
than WAL across nodes.

Rows also record `substrate`, `binding_permutation`, `binding_mode`,
`coordination_isomer`, `donor_maps`, `substrate_site` and `resonance`. Binding
permutations are numbered from 1 per ligand: each donor subset / coordination
isomer / substrate seating has its own number, shared across its substrate sites,
endpoint conformers and TS rows. Numbers remain stable when preparation is repeated.
The 20 kcal/mol coordination-screening cutoff still applies; arrangements outside
it are not submitted. Every successfully completed ORCA endpoint/TS has a Gibbs
energy; failed or pending rows keep their status rather than implying a result.

```sql
SELECT ligand_id, binding_permutation, substrate, substrate_site,
       pair_index, endpoint, conformer, status, gibbs_free_energy_eh
FROM jobs
ORDER BY ligand_id, binding_permutation, substrate_site, pair_index, endpoint, conformer;
```

Opening an older database adds these columns; rerun `prepare` to backfill them
from the original sampling records. Historical alpha/beta runs retain their actual
original substrate and geometries.

Endpoint jobs run ORCA optimization and a **separate NumFreq**, reject only graph
changes accompanied by a distance shift greater than 2 Å, and retain valid minima
within 3 kcal/mol of the lowest electronic energy. Minima require zero imaginary
frequencies. If the first endpoint NumFreq finds any imaginary frequencies, run
VeryTightOpt from that optimized geometry, check connectivity again, and repeat
NumFreq once. TS jobs do not use this recovery.
A TS candidate waits for its endpoint jobs and runs direct OptTS followed by
separate NumFreq, which must give exactly one imaginary frequency. Direct candidates
use up to `calculation.ts_attempts` OptTS attempts (two by default), stopping on
convergence. Each repeat uses the last saved geometry, or the previous input if
no geometry was saved, and calculates a fresh numerical Hessian. If these attempts
fail, or the subsequent NumFreq fails or does not give exactly one imaginary mode,
restart from the **original TS guess** and run constrained optimization, fixing the
Cartesian positions of the donor, transferring hydrogen and acceptor. Repeat that
optimization once if needed, continuing from its last geometry. Only after it
converges, run up to two further OptTS attempts, then NumFreq on convergence.
The retry budgets for all three phases use `calculation.ts_attempts`. Numbered
OptTS and NumFreq artifacts preserve both the initial and recovery calculations.
Among successful
candidates for each ligand/pair/seating, only the lowest Gibbs-energy TS has
`retained=1`; all candidate results remain in the database. The fallback waits
until every candidate finishes. If any succeeds, the fallback is marked `skipped`;
otherwise it runs NEB-TS followed by NumFreq. Historical single-guess jobs still try direct OptTS,
then constrained relaxation/OptTS retries, then NEB.
NEB uses the original coembedded reactant and product geometries with the rxts
TS guess, preserving their shared atom mapping and alignment.

All ORCA stages use `oet_gxtb_alpb -e XTB_EXE -s h2o`: gas-phase g-xTB plus the
GFN2(ALPB water) minus GFN2(gas) correction. Gibbs energies come exclusively from
ORCA NumFreq thermochemistry; no solution
standard-state correction is added. The parser takes the reference electronic
energy from thermochemistry, avoiding displaced NumFreq single-point energies.
See the [ORCA external-program tutorial](https://www.faccts.de/docs/orca/6.1/tutorials/workflows/extopt.html)
and [thermochemistry manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html).

Failed rows are left for inspection. For sampling failures, correct the missing
or invalid seeds and rerun `rxflow prepare`; it reactivates repaired sampling rows
without resetting completed calculations. For ORCA failures, stop any previous
worker for the row and retry explicitly; attempt-specific filenames preserve
earlier artifacts:

```sh
rxflow retry --id 42
rxflow resume --workers 100
```

The separate precursor workflow prepares ligand–metal complexes and free,
fully protonated neutral substrates for the eventual formation-energy steps:

```sh
rxflow --config hpc_ED_dF.toml precursors submit
rxflow --config hpc_ED_dF.toml precursors status
```

Use `hpc_ED_pyr.toml`, `hpc_ni.toml` or `hpc.toml` for the other runs. The reaction
workflow must already be prepared: its database supplies the binding-permutation
IDs, donor atoms and coordination-isomer identities. Each binding permutation gets
fresh embedding with **10 conformers per state** by default (`--conformers 10`).
The deprotonated ligand–metal complex has charge **+1**; protonating each explicitly
anionic ligand site in the input CSV gives a separate **+2** complex. The TS
H-acceptor does not determine this protonation site. Complexes are built from
ligand and metal graphs, without extracting or reusing reaction geometries.
Each distinct neutral substrate is freshly embedded once per run with 10 conformers
and shared across the relevant ligands. Missing embeddings are recorded as sampling
failures rather than silently reducing the requested number of jobs.

Submission creates a preparation job followed by a dependent worker array (100
workers by default; `--workers N` changes this). Workers use the configured ORCA
minimum-optimization and separate NumFreq protocol, including optimization
restarts and the VeryTightOpt/NumFreq recovery for imaginary modes. The solvent,
external method and thermochemistry settings are the same as the reaction run.

Outputs stay under the existing results directory:

```text
ED_dF_hpc_results/
  precursors.sqlite
  precursors.json
  ligand_01/complex/
    binding_001_deprotonated_conformer_001_input.xyz
    binding_001_protonated_site_017_conformer_001_input.xyz
  substrate/
    substrate_001_neutral_conformer_001_input.xyz
```

The reaction `workflow.sqlite` and its results are untouched. In the separate
`precursors.sqlite`, `precursor_links` records every existing ligand/pair/seating,
binding permutation, substrate site and resonance route. The `precursor_minima`
view selects the **lowest NumFreq Gibbs energy among all succeeded conformers**
for each precursor species, independent of the electronic-energy retention window.
The `precursor_diagram_references` view joins those energies back to every reaction
pair and exposes both complex-plus-substrate sums in Eh. Shared precursors reuse
one energy while distinct reaction routes remain separate; missing results are NULL.
The reaction-diagram notebook reads these links and adds the precursor and acid/base
levels described below. Isolated free ligands and bare metals are not calculated.

For individual stages in a compute-node allocation, or to repair/restart jobs:

```sh
rxflow --config hpc_ED_dF.toml precursors prepare
rxflow --config hpc_ED_dF.toml precursors worker
rxflow --config hpc_ED_dF.toml precursors retry --id 42
rxflow --config hpc_ED_dF.toml precursors resume --workers 100
```

`prepare` is repeatable and repairs sampling failures without resetting completed
jobs. Keep the same conformer count and protocol when rerunning it. Stop a row's
old worker before retrying a failed or interrupted ORCA calculation; `resume`
submits only workers for pending rows. Precursor IDs belong to `precursors.sqlite`,
so use the `precursors` commands when checking or retrying these jobs.

Acetic acid and acetate provide the shared acid/base reference pair. Prepare them
once, then submit their pending optimization/NumFreq jobs:

```sh
rxflow --config hpc_ED_dF.toml acid-base prepare
rxflow --config hpc_ED_dF.toml acid-base resume --workers 20
rxflow --config hpc_ED_dF.toml acid-base status
```

Alternatively, `acid-base submit --workers 100 --partition kemi1` submits preparation
and workers together. It can be repeated after earlier acid/base runs have finished.
Each species gets 10 fresh
conformers: acetic acid `CC(=O)O` has charge 0; acetate `CC(=O)[O-]` has charge -1;
both use singlet multiplicity. Outputs are shared in `solution_references/`, with
`references.sqlite` and separate `acetic_acid/` and `acetate/` directories. All
configs in the same project use this set; an optional `paths.references` can point
to a different reference directory. Calculation settings are frozen and must
match the precursor settings before the energies can be combined. Preparing these
references leaves running precursor and reaction jobs untouched.

The same `acid-base` command also prepares **protonated bound reactants** for
the selected run: ligand + metal + the fully protonated, coordinated substrate,
with charge **+2**. It freshly embeds 10 conformers from the native protonated
CSV ligand graph, retaining the recorded donor slots and substrate-to-metal
bond. It runs the existing gxTB/ALPB(water) minimum optimization and NumFreq
workflow. These calculations are stored in the run's `acid_base.sqlite`, with
files under `ligand_XX/complex/protonated_bound/`. Different deprotonated resonance
pathways for the same binding arrangement, substrate site and native protonation
site share one protonated bound state. `bound_reactant_links` records every
original pair/resonance mapping; none of those pathways is discarded.

To add this state to an already completed run, use:

```sh
rxflow --config hpc_ED_dF.toml acid-base submit --workers 100 --partition kemi1
rxflow --config hpc_ED_dF.toml acid-base status
# After the new gxTB optimization/NumFreq jobs finish:
rxflow --config hpc_ED_dF.toml singlepoint submit --workers 100 --partition kemi1
```

Completed acid/base references and existing single points are reused. Repeating
`acid-base submit` discovers missing bound-state jobs; `resume` only drains
already prepared pending jobs. Use `acid-base retry --bound --id ID` to retry a
failed bound complex; without `--bound`, retry IDs refer to the shared references.

`precursors status` exposes the balanced totals for every pair, selecting the
lowest successful NumFreq Gibbs energy for each reference:

| State | Total Gibbs energy |
| --- | --- |
| Protonated ligand–metal complex + free substrate + acetate | G(protonated complex) + G(substrate) + G(acetate) |
| Deprotonated ligand–metal complex + free substrate + acetic acid | G(deprotonated complex) + G(substrate) + G(acetic acid) |
| Protonated bound reactant + acetate | G(protonated bound reactant) + G(acetate) |
| Bound reactant complex + acetic acid | G(reactant) + G(acetic acid) |
| Transition state + acetic acid | G(TS) + G(acetic acid) |
| Bound proton-transfer product + acetic acid | G(product) + G(acetic acid) |

Each state has total charge +1 and the same atom inventory. The substrate is
already included in the bound reactant, TS and product; no free-substrate term is
added to those states. Missing energies remain NULL. These totals are thermodynamic
states; the acid/base reference calculations do not locate a proton-transfer
transition state involving acetate. `hpc_precursors.balanced_diagram_references`
provides the energy links; `load_diagram_inputs` reads all three databases for the
notebook without changing their results. All totals use the existing ORCA thermochemistry convention,
without a separate solution standard-state or concentration correction.

Open `reaction_diagrams.ipynb`, set `DB_PATH` to the run's `workflow.sqlite`, and
rerun all cells to refresh the results. `PRECURSOR_DB_PATH` defaults to the adjacent
`precursors.sqlite`; `REFERENCE_DB_PATH` points to the shared acid/base database
and can be overridden for a different reference protocol. The notebook draws the
balanced states, preserving every complete reaction pair/resonance and every
recorded native ligand protonation site. Each ligand uses a single zero: the lowest
available protonated LM + substrate + acetate total across its binding permutations.
Once protonated bound states are prepared, the plot adds the alternative
association-before-deprotonation route. Both formation routes connect the same
unbound protonated start and bound deprotonated reactant; the two alternative
intermediates are not connected to each other. Missing levels remain blank.
Tables outside the axes list selected job IDs, conformers, absolute Gibbs energies,
calculation statuses and missing precursor inputs. The numerical plot data are also
available in the `diagram_rows` variable.

For the bound proton-transfer step, forward barriers are G(TS) − G(reactant) and
reverse barriers are G(TS) − G(product). A pathway qualifies when both are at most
25 kcal/mol (`BARRIER_LIMIT_KCAL` in the notebook). The report lists qualifying
pairs/permutations with both values; contrasting plot colors distinguish each
site and barrier class. Plot legends and ligand-grid captions give both barriers
without repeating a threshold tag. Tables have no colored shading, grid captions
use black text, and the grid section headings identify the threshold. The
blue/red/purple atom-role highlights are preserved. All other pathways remain visible. Classification uses
unrounded energies; the numeric reports are also available as `barrier_report`
and `qualifying_pathways`.

Missing precursor levels stay blank, without connectors across gaps. Until a
protonated starting total is available, the common zero explicitly falls back to
the lowest bound reactant + acetic acid. If acetic acid itself is unavailable, the
notebook shows the original three-state diagram relative to the lowest reactant.
Permutation coverage continues to count complete reactant–TS–product pathways
independently of precursor readiness. Pending candidates do not hide successful
results; rerunning may update the minima while calculations are still active.

The single-point workflow evaluates **PBE0/def2-TZVP-D4 with SMD(water)** on the
currently selected lowest-G(gxTB) geometry for every state. It includes reactants,
transition states, products, both ligand–metal protonation states, protonated
bound reactants from `acid_base.sqlite`, free substrates,
acetic acid and acetate. It uses all succeeded NumFreq results when selecting a
geometry, ignoring `retained`; it does not evaluate every conformer or re-optimize
the selected structures.

```sh
rxflow --config hpc_ED_dF.toml singlepoint plan
rxflow --config hpc_ED_dF.toml singlepoint submit --workers 10
rxflow --config hpc_ED_dF.toml singlepoint status
```

Use the corresponding config for each other run. Shared acid/base single points
are prepared once in `solution_references/singlepoints.sqlite` and reused across
runs. Each run has its own `singlepoints.sqlite` for reaction and precursor jobs.
The source gxTB databases and calculated energies are unchanged. ORCA inputs,
outputs and frozen XYZ snapshots are kept in a `singlepoint/` subfolder beside
each source state's optimized geometry. Scratch wavefunctions and integral files
are not copied back. The input uses ORCA's [SMD water model](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html):

```text
! PBE0 def2-TZVP D4 SP
%cpcm
 smd true
 SMDsolvent "water"
end
```

SCF convergence settings and iteration limits use ORCA's defaults, with no
explicit SCF override block. Charge and multiplicity come from the selected gxTB job. An SP is accepted only
after normal ORCA termination, converged SCF and a finite final energy. Its
corrected Gibbs energy is computed using the **same source NumFreq result**:

```text
G(corrected) = G(gxTB NumFreq) - E(gxTB NumFreq) + E(PBE0-D4/SMD)
```

The stored gxTB electronic energy is the thermochemistry reference, not a displaced
NumFreq single point. Source database, job ID, source attempt, geometry hash,
charge/multiplicity, electronic energy and Gibbs energy are recorded with each SP.
A changed source invalidates the previous correction. Repeating `prepare` or
`submit` adds new selected sources without resetting completed calculations.

Single-point resources are separate from `[calculation]`. The defaults are
**4 CPUs, 8000 MB per worker, a 24-hour timeout, and 10 workers**. Override the
per-worker allocation in each config if needed:

```toml
[singlepoint]
cpus = 4
memory_mb = 8000
timeout_seconds = 86400
```

`singlepoint submit` prepares the small input snapshots locally, then submits the
Slurm worker array. Workers atomically claim the shared references first, then
their run's SP jobs, so separate run submissions do not duplicate the shared work.
The supplied HPC configs set `executables.mpirun_exe` to
`/software/kemi/openmpi/openmpi-4.1.1/bin/mpirun`, matching `~/bin/scripts/submit`.
They set `executables.singlepoint_orca_exe` to
`/groups/kemi/julius/orca_6_1_0/orca` from the same script. This setting applies
to DFT single points; if omitted, they use `executables.orca_exe`.
The Slurm allocation also follows that script: one node, four tasks and one CPU
per task for the default four-process calculation, with 8000 MB total memory and
requeue disabled. ORCA is launched directly with `--use-hwthread-cpus`; ORCA
starts its own MPI processes. As in the submission script, MPI binary and
library directories are prepended to `PATH` and `LD_LIBRARY_PATH`, together with
the configured ORCA directory. On each compute node, workers use the configured
MPI installation if available, then `executables.mpirun_fallback_exe` (set to
`/usr/lib64/openmpi/bin/mpirun` in the supplied configs), then an inherited
`mpirun` on `PATH`. The selected MPI installation supplies both binary and
library paths; no environment modules are loaded. If no launcher is available,
the worker exits before claiming any jobs, leaving calculations pending.
Calculations use separate temporary directories under `$SCRATCH`, or
`$SLURM_TMPDIR` if `$SCRATCH` is unset, copying results back before cleanup.
This also works when MPI was absent from the
submitting shell's `PATH`. Other clusters can set the corresponding executable
path; if omitted, the launcher is looked up as `mpirun` on `PATH`.
For individual stages and retries (stop an old worker before resetting its row):

```sh
rxflow --config hpc_ED_dF.toml singlepoint prepare
rxflow --config hpc_ED_dF.toml singlepoint worker
rxflow --config hpc_ED_dF.toml singlepoint retry --id 42
rxflow --config hpc_ED_dF.toml singlepoint retry --id 1 --references
rxflow --config hpc_ED_dF.toml singlepoint resume --workers 10
```

`resume` submits workers for existing pending rows. If gxTB minima changed or
additional states completed, rerun `prepare`/`submit` to select their new geometries.
SP retry IDs belong to the run's `singlepoints.sqlite`, or to the shared database
when `--references` is supplied. Earlier attempt outputs are preserved.

In `reaction_diagrams.ipynb`, `ENERGY_MODE = "auto"` reports DFT completion and
keeps the whole run on gxTB energies until all currently selected states have
matching successful SP results. It then switches all components, the ligand's
common reference, and both barrier directions together. The five-state plot never
combines corrected components with uncorrected components. Set `ENERGY_MODE` to
`"gxtb"` to retain the original plots, or `"dft"` to explicitly inspect partial
corrected results with missing states left unavailable. Rerun all cells to refresh.
Plots and numerical data identify the active model and retain the selected source
and single-point job IDs.

Tests exercise input/output parsing, subprocess artifact handling, queue claims,
endpoint selection, TS fallback, precursor charges and pair links, and Slurm job
dependencies with mocks. Precursor sampling is also checked with real embedding;
notebook energy accounting and pathway links have been checked against completed
cluster results from all four runs.

```sh
python -m unittest test_hpc_workflow test_hpc_precursors test_hpc_singlepoint test_reaction_diagrams test_substrate_hydrogens
```
