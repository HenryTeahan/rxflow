# rxflow

A small transition-metal reaction workflow using `ligands.csv`, `rxembed` and `rxts`.

The isolated environment is `.venv` (Python 3.13). `uv sync --locked` installs the
locked dependencies and editable sibling `../rxembed` / `../rxts` checkouts.

```sh
source .venv/bin/activate
jupyter lab manual_workflow.ipynb
```

Select **Python (rxflow)**. The notebook checks its interpreter, prints the local
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

`hpc_workflow.py` adds a separate ORCA worker, using the existing local sampling
as its starting point. It imports the three initial conformers per endpoint and
the selected coembedding/rxts TS guess; it does not import local energies or run
the local optimizer on the HPC. Preparation can be repeated as more local guesses
finish, without resetting existing jobs.

```sh
.venv/bin/python hpc_workflow.py prepare
.venv/bin/python hpc_workflow.py status
```

Configure paths and calculation settings in `hpc.toml`. Copy the code and entire
`hpc_results/` directory to the cluster, install the Python environment, and load
the cluster's ORCA/MPI modules. The prepared database no longer needs `results/`.
The current `uv.lock` also needs the sibling `rxembed` and `rxts` checkouts.
Calculation settings are frozen per row during preparation; use a new results
directory when changing the protocol.
Set the actual executable paths, then submit from the project directory:

```sh
export ORCA_EXE=/path/to/orca
export OET_GXTB_ALPB_EXE=/path/to/oet_gxtb_alpb
export XTB_EXE=/path/to/xtb
sbatch hpc_worker.sbatch
```

The example array runs up to 100 independent one-core workers. Match its memory,
CPU and time requests to `hpc.toml` when changing them. Scratch uses Python's
temporary directory (`TMPDIR` on clusters); set it to node-local storage.
The existing local runner and notebook remain unchanged.

Each worker atomically claims a row in `hpc_results/workflow.sqlite`. Its single
`jobs` table holds one row per endpoint conformer and one TS per ligand/pair/seating,
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
frequencies. A TS waits for its endpoint jobs, constrains the transferring H and
its donor/acceptor during relaxation, tries OptTS twice, and uses NEB-TS if both
attempts fail. Separate NumFreq must then give exactly one imaginary frequency.
NEB uses the lowest valid optimized endpoint on each side.

All ORCA stages use `oet_gxtb_alpb -e XTB_EXE -s h2o`: gas-phase g-xTB plus the
GFN2(ALPB water) minus GFN2(gas) correction. Gibbs energies come exclusively from
ORCA NumFreq thermochemistry; no solution
standard-state correction is added. The parser takes the reference electronic
energy from thermochemistry, avoiding displaced NumFreq single-point energies.
See the [ORCA external-program tutorial](https://www.faccts.de/docs/orca/6.1/tutorials/workflows/extopt.html)
and [thermochemistry manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html).

Failed rows are left for inspection. After stopping any previous worker for a
row, retry explicitly; attempt-specific filenames preserve earlier artifacts:

```sh
.venv/bin/python hpc_workflow.py retry --id 42
.venv/bin/python hpc_workflow.py worker
```

ORCA is not installed locally. Tests exercise input/output parsing, subprocess
artifact handling, queue claims, endpoint selection and TS fallback with mocks;
the first real ORCA calculation still needs validation on the cluster.

```sh
python -m unittest discover -s tests
```
