# Workflow reference

[Back to the workflow guide](../README.md)

## Screening and atom identities

`generate_reaction_pairs()` enumerates donor subsets and coordination isomers,
embeds one screening geometry per arrangement, and evaluates gas-phase g-xTB
single points. Arrangements within 20 kcal/mol of the lowest electronic energy
produce mapped reactant/product pairs over substrate sites and resonance forms.

The CSV provides a substrate for each ligand. All symmetry-distinct aromatic
C–H sites are considered; equivalent sites share a representative and alkyl C–H
sites are excluded. `substrate_site` is the bound carbon's one-based index in the
parsed substrate before adding hydrogens. Input atom maps do not affect symmetry.

Charges come from the molecular graph. Atom maps start at 1; Python, ORCA TS-mode
and stored donor/H/acceptor indices start at 0. Endpoint construction carries the
screened ligand stereochemistry into both endpoint Isomers while preserving
mapped atom identities and coordination seating.

## Embedding and pruning

Reaction endpoints and TS conformational sampling each request 20 conformers.
Precursors, encounters and acid/base species request 10 by default. These counts
are budgets, not a guarantee of that many surviving calculations.

Ordinary embedding steps call `prism_pruner.pruner.prune_by_rmsd` using heavy
atoms and `sampling.prune_rmsd` (default 0.5 Å). No energies or additional maximum
atom-displacement criterion are supplied, and pruning does not minimize or change
the retained coordinates. Defined stereo/coordination ensembles stay separate.
Coembedded pairs and coembedding retries are not pruned.

Endpoint sampling records the surviving counts, so preparation imports the
survivors rather than reporting deliberately pruned slots as missing seeds.
Unused precursor/reference slots have `status='skipped', stage='pruned'`.

## Endpoint optimization and thermochemistry

ORCA performs g-xTB/ALPB(water) optimization followed by a separate NumFreq.
The external executable is called as `oet_gxtb_alpb -e XTB_EXE -s h2o`: gas-phase
g-xTB plus the GFN2(ALPB water) minus GFN2(gas) correction.

For ordinary minima, the workflow rejects intended bond-length changes greater
than 1.5 Å from the starting geometry and changes in native metal coordination.
This is a change-from-reference test, not a maximum allowed metal–donor length.
Valid minima require zero imaginary frequencies. An initial imaginary-frequency
result receives one VeryTightOpt/connectivity/NumFreq recovery.

Endpoint `retained` flags use a 3 kcal/mol **electronic-energy** window. The
minimum used for DFT and diagrams is selected by successful NumFreq **Gibbs**
energy, independent of that flag. TS DFT selection has its own Gibbs-energy window.

Gibbs energies and their reference electronic energies come from ORCA NumFreq
thermochemistry. No separate solution standard-state or concentration correction
is added. The parser avoids displaced NumFreq single-point energies.

## TS guesses, validation and recovery

For each ligand/pair/seating, coembed 10 endpoint pairs and use up to the five
lowest `ts_rmsds` pairs for rxts GFN-FF scans. Each guess must contain the expected
C–H breaking and H–acceptor forming bonds. Around its selected peak, C–H must
lengthen while H–acceptor shortens. Failed guess checks permit two re-coembedding
retries with different seeds; profiles, endpoints and checks remain available.

Each accepted guess follows this sequence:

1. **Constrained preoptimization:** freeze donor, transferring H, acceptor,
   central metal and the shortest mapped acceptor-to-metal path. Cap at 50 cycles;
   a cycle-limited calculation may supply its last usable geometry.
2. **OptTS:** release constraints and use `calculation.ts_attempts` attempts
   (default two), restarting from the latest geometry with a fresh numerical Hessian.
3. **Separate NumFreq and graphRC:** require exactly one imaginary frequency.
   The mode's significant bond, angle and dihedral changes must collectively
   include donor, H and acceptor; additional atoms are allowed.

If these attempts fail, score every frame of the saved rxts scan using gas-phase
g-xTB single points. Select the highest checked interior local maximum, excluding
endpoints and monotonic shoulders. Repeat preoptimization and the OptTS budget
from that peak, with the same frequency/mode checks.

When all direct candidates fail, the group NEB fallback uses a checked scored
peak and its matching coembedded endpoints. Optimize those endpoints separately
with the normal minimum optimization/connectivity checks, then run NEB-TS and
validate its frequency/mode. Historical single-guess jobs use inline NEB recovery.

Once the initial group searches finish, reembed its lowest-G validated TS using
`rxembed.embed(..., n=20, fix=[donor, hydrogen, acceptor])`. Preserve the input
coordination/stereochemistry, rigidly align the frozen core to the parent frame,
check its coordinates/maps, prune, and run each survivor through the same
preoptimization/OptTS/NumFreq/mode sequence. Sampled TSs have no saved rxts scan
and do not start scan or NEB recovery. The group is sampled once; original and
sampled successful TSs compete for the lowest Gibbs energy.

## DFT selection and energy accounting

DFT single points use PBE0/def2-TZVP-D4/SMD(water) without geometry optimization.
For each ligand/pair/seating, select every successful original or sampled TS with
`G(g-xTB NumFreq) − G(group minimum) ≤ 3 kcal/mol`, including the boundary.
Prune these optimized TS geometries with Prism at 0.5 Å heavy-atom RMSD, in
ascending Gibbs-energy order. Keep the lowest-G representatives and preserve
distinct molecular identities. Both preparation and notebook readiness use this
same selection; this catches distinct embeddings that converge to the same TS.
For encounter states, select every successful RMSD-distinct conformer within 5 kcal/mol of the gxTB NumFreq minimum. The target is 20 encounter conformers per ligand–binding arrangement (`sampling.encounter_conformers`); increasing the target appends seeds while preserving existing results. Other states contribute their lowest-G successful conformer.

Each selected geometry has its own correction:

```text
G(corrected) = E(DFT) + G(g-xTB NumFreq) − E(g-xTB NumFreq)
```

The final TS minimum is chosen from the corrected energies. `ENERGY_MODE='auto'`
waits until the entire selected set has matching DFT results before switching
the diagram. `'gxtb'` keeps the original energies; `'dft'` permits partial
corrected results with missing states left unavailable.

Single points record the source database, job, attempt, geometry hash, charge,
multiplicity, method and source energies. Changed sources invalidate corrections.
Identical frozen elements/coordinates with matching charge, multiplicity and method can reuse
a completed electronic energy, combined with the new source's thermal correction.
XYZ comment lines and numeric formatting are ignored for reuse. Existing raw-text
hashes remain valid; historical completed rows are not rewritten.
Pending TS corrections outside the current window or removed by optimized RMSD
pruning are skipped; completed calculations are preserved. Skipped sources can
be reactivated when selected again. `singlepoint resume` always refreshes this selection.

After cancelling the old SP arrays and waiting for their removal from `squeue`,
use `singlepoint resume --reset-running --workers 10 --partition kemi1` with the
run's config. This verifies the recorded workers have stopped, creates SQLite
backups and a row manifest under `singlepoint_restarts/`, resets interrupted
running rows while preserving attempt numbers, then prepares and submits the
remaining selected work. Ongoing selected calculations are left running unless
this stopped-worker recovery is explicitly requested.

## Precursors, encounters and acid/base states

The regular submit/resume chain prepares ligand–metal protonation states, neutral
free substrates and uncoordinated substrate encounters. These use fresh graphs;
equivalent site/resonance routes share reference calculations.

For the supplied chemistry, deprotonated ligand–metal complexes have charge +1;
protonating each explicitly anionic CSV ligand site gives separate +2 complexes.
The protonation site is independent of the TS H acceptor. Encounters contain the
+1 complex and neutral substrate, sampled with rxembed noncovalent contact modes.
Encounter optimization has no geometry or fragment-separation filter; convergence
and zero imaginary frequencies are required. Native contact-mode failures are
isolated so surviving modes still supply conformers.

The separate `acid-base` workflow calculates shared acetic acid (neutral) and
acetate (−1), plus each run's protonated bound reactants (+2). All are singlets.
Shared references default to `solution_references/`; override `paths.references`
for a different protocol. Their calculation settings must match the run.

Balanced diagram states include:

| State | Gibbs-energy components |
| --- | --- |
| Protonated unbound start | Protonated ligand–metal + substrate + acetate |
| Deprotonated unbound state | Deprotonated ligand–metal + substrate + acetic acid |
| Protonated bound intermediate | Protonated bound reactant + acetate |
| Encounter, reactant, TS or product | Corresponding complex + acetic acid |

The complex already contains its substrate; do not add another free-substrate
term. These states share total charge +1 and atom inventory. Missing reference
energies stay unavailable. Acid/base reference calculations do not locate an
acetate-mediated proton-transfer TS.

## Recovery and individual stages

Run foreground sampling/preparation inside a compute allocation. The scheduled
`submit`/`resume` commands perform this work on Slurm.

```sh
rxflow --config my_run.toml sample --ligand-index 0
rxflow --config my_run.toml prepare
rxflow --config my_run.toml worker

rxflow --config my_run.toml precursors prepare
rxflow --config my_run.toml precursors worker

rxflow --config my_run.toml singlepoint plan
rxflow --config my_run.toml singlepoint prepare
rxflow --config my_run.toml singlepoint worker
```

`--ligand-index` is a zero-based CSV row; omit it to sample all rows sequentially.
Use the same conformer count/protocol when repeating precursor preparation.
Keep the seed source for saved scans and guess regeneration during recovery.

| Situation | Action |
| --- | --- |
| Missing/invalid reaction seeds (`stage='sampling'`) | Repair the seeds and repeat `prepare` |
| Failed or interrupted calculation | Stop its old worker, use the matching queue's `retry --id ID`, then resume |
| Failed TS calculations or TS sampling tasks | Inspect/reset with `retry-ts --dry-run` / `retry-ts` |
| Verify saved successful TS modes | `verify-ts --dry-run`, then `verify-ts` to queue rejected modes |
| Recheck saved endpoint geometries | `filter-endpoints --dry-run`, then `filter-endpoints` |
| Encounter retry | `precursors retry --encounters --id ID` |
| Protonated bound-state retry | `acid-base retry --bound --id ID` |
| Shared reference DFT retry | `singlepoint retry --references --id ID` |
| New minima/reference results after DFT preparation | Repeat `singlepoint submit` to refresh selection |

All actions above take `rxflow --config my_run.toml` as their prefix.
Retry IDs belong to the corresponding database. Mode audits and endpoint filters
archive prior rows/database snapshots. They preserve passing calculations.
`ts-recovery-status` reports the latest TS reset batch.

Slurm `DependencyNeverSatisfied` means an upstream dependency cannot succeed.
Inspect `sacct` and the upstream log, correct the failure, then cancel the obsolete
blocked chain and submit a replacement. `resume` creates new dependencies; it
does not repair the old Slurm chain or automatically reset interrupted rows.

## Storage, resources and provenance

Keep databases on shared storage that supports SQLite file locking. Queues use
rollback journaling rather than WAL across nodes. Schema upgrades are performed
atomically and rechecked under the migration lock. Calculation settings are
frozen per row; use fresh result paths for a different chemistry/protocol.

`workflow.sqlite` records binding permutation, donor maps, substrate/site,
resonance, stages, attempts and energies in Eh. Binding permutations remain stable
across preparation. TS sampling tasks use `ts_strategy='sampling', conformer=-1`;
sampled calculations use `'sampled'` and conformers −2, −3, etc. `ts_parent_id`
and `ts_sampling_json` preserve the parent, seed, constraints and pruning report.

Reaction artifacts live under `ligand_XX/{reactant,product,ts_guess}/`.
Precursor/reference links are in their separate databases. DFT inputs, frozen XYZs
and outputs are saved in `singlepoint/` beside the source geometry. Shared
reference DFT results are in `solution_references/singlepoints.sqlite`.

Reaction scratch uses `$TMPDIR`. DFT scratch uses `$SCRATCH`, then
`$SLURM_TMPDIR`, then the default temporary directory. Completed artifacts are
copied back; scratch wavefunctions and integrals are not retained.

```toml
[singlepoint]
cpus = 4
memory_mb = 8000
timeout_seconds = 86400
```

DFT resources are separate from `[calculation]`. `executables.singlepoint_orca_exe`
overrides the normal ORCA executable for DFT. MPI lookup tries `mpirun_exe`, then
`mpirun_fallback_exe`, then `PATH`; the selected installation supplies its binary
and library paths. Missing MPI leaves jobs pending before any claim.

## Notebooks and local exploration

In `reaction_diagrams.ipynb`, set `DB_PATH` and rerun all cells. Adjacent precursor
and shared reference databases are read without modifying their calculation rows.
The notebook preserves pathway/resonance identities and exposes `diagram_rows`,
`barrier_report`, `reactant_ts_report` and `low_energy_failed_paths`.
The final **Likely active ligands** section exposes `likely_active_ligands` and
`active_pathways`. A qualifying pathway has both barriers ≤25 kcal/mol and its
lower bound/encounter reactant within 3 kcal/mol of the ligand's lowest available
bound or encounter reactant across all comparable groups, including groups
without a TS. Substrate, charge and multiplicity stay the same in each comparison.
A competing pathway in the same reactant window vetoes the entire ligand when
its product is below an otherwise qualifying pathway's product and its reverse
barrier is >25 kcal/mol, even if its forward barrier is >25 kcal/mol. The final
section lists these exclusions; `thermodynamic_traps` and
`trap_disqualified_ligands` expose the reasons and excluded ligand IDs.

Forward barriers use `G(TS) − min[G(encounter), G(bound reactant)]`; when the
encounter is missing, the bound reactant is used. Reverse barriers use
`G(TS) − G(product)`. `BARRIER_LIMIT_KCAL` defaults to 25. `TS_FAILURE_WINDOW_KCAL`
defaults to 3 for prioritizing missing TS paths near low-energy reactants; this
reporting window is distinct from the per-group DFT conformer-selection window.

For local exploration, [`src/run_workflow.py`](../src/run_workflow.py) performs gas-phase screening,
endpoint sampling and TS-guess generation; it does
not schedule the full HPC/ORCA/DFT chain. The older
[`legacy/reaction_workflow.py`](../legacy/reaction_workflow.py) remains a reference
requiring the absent `workflow_store` module. Dependency revisions
and package sources are recorded in `pyproject.toml` and `uv.lock`.
Run it from the installed environment with `python -m run_workflow --help`.
