# Portable ORCA/UMA calculations

`rxflow uma` starts a separate method run from the lowest PBE0-corrected
conformer for every plotted reactant, TS, product, and encounter arrangement.
It preserves all ligand/site/pair/permutation identities, including pathways
outside the system-comparison sink filter. Minima from incomplete pathways are
also exported so the global reactant/encounter reference is not lost.

The source databases are read only. UMA performs its own optimization and
NumFreq calculations. The UMA Gibbs energy is the UMA electronic energy plus
its **own** ORCA NumFreq thermal correction; it is never combined with a gxTB
thermal correction or a PBE0 electronic energy.

## Prepare on the current cluster

Add an `[uma]` section to an existing system config, or adapt
[`uma.example.toml`](../uma.example.toml). The executable, precise model label,
arguments, and CPU/CUDA selection are configurable so the eventual UMA wrapper
instructions can be applied without changing the workflow. Executables and
model weights are only required on the cluster that runs the calculations.
The wrapper owns the solvent/environment treatment; gxTB ALPB and PBE0 SMD
settings are not implicitly passed to UMA.

```sh
rxflow --config uma.example.toml uma plan
rxflow --config uma.example.toml uma prepare
```

By default the output is `Ni_hpc_results/uma/` for the example. Preparation copies
the selected XYZs, mapped molecular graphs, TS atom indices, available recovery
scans and matching endpoints, and balanced reference links. The generated
`uma_manifest.json` records the frozen selection and geometry hashes.
`run_uma.toml` is a portable configuration whose paths refer only to this folder.
No calculations are submitted by `plan` or `prepare`.

The default also selects the lowest precursor, substrate, protonated-bound and
acid/base reference molecules used by the balanced diagrams, keeping the whole
diagram in one energy model. Set `include_balanced_references = false` to run
only the four core state classes; unavailable balanced levels then stay blank.

Preparation is idempotent and does not refresh a frozen run from subsequently
changed source minima. Use a new UMA results directory for a different model,
different chemistry arguments, or a new source selection. This protects completed
work and keeps CPU/GPU transfers independent of the source databases.

## Run on another cluster

Copy the **whole prepared folder**, rather than individual ORCA outputs:

```sh
rsync -a Ni_hpc_results/uma/ gpu-cluster:/scratch/project/Ni_TM_uma/
```

Install the rxflow revision containing this workflow on the destination cluster,
with the normal project dependencies (including graphRC and rxts). Install/load
ORCA and the chosen UMA external wrapper and model there. rxflow does not need
to install PyTorch/UMA into its own environment: the external wrapper can use
its separate environment. Configure the exact UMA arguments as instructed for
that wrapper. If a server/client wrapper is used, its server must be running
in the worker's allocation; a wrapper script can manage that setup.

```sh
export ORCA_EXE=/destination/path/to/orca
export ORCA_UMA_EXTERNAL_TOOLS=/destination/path/to/uma_wrapper
cd /scratch/project/Ni_TM_uma
rxflow --config run_uma.toml uma submit --partition YOUR_GPU_PARTITION --workers 4
rxflow --config run_uma.toml uma status
```

`uma submit --dry-run` prints the exact Slurm command without submitting it.
Each array worker requests one GPU by default and drains the portable queues.
For CPU testing, set `uma.gpus_per_task = 0` **and** select CPU in the wrapper.
Wrapper arguments may contain environment variables, e.g. `${UMA_DEVICE}`;
these are expanded on the destination worker, so the same frozen protocol can
be tested on CPU and then run on CUDA without embedding cluster-specific paths.
Model weights should be available in the wrapper's destination cache.

To run within an allocation rather than submit an array:

```sh
rxflow --config run_uma.toml uma worker
```

Workers never read the original source databases, CSV or seed tree. They reuse
the existing optimization/restart, separate NumFreq, and TS mode-checking code.
TSs follow constrained preoptimization → OptTS → NumFreq/graphRC validation,
with the usual group NEB fallback and TS conformer sampling. Archived rxts scans,
when available, are rescored with **UMA**, not gxTB. When a selected TS has no
archived scan, its OptTS and matching-endpoint NEB fallback remain available.

Minima must pass the connectivity/coordination and bond-length checks and have
zero imaginary frequencies. UMA encounter optimizations also pass these checks.
TSs require exactly one imaginary frequency involving the recorded transfer
atoms; changes outside the two intended H-transfer bonds are rejected.

Failed jobs are retained, with no automatic substitution of another method:

```sh
rxflow --config run_uma.toml uma retry --id 123 --queue reaction
rxflow --config run_uma.toml uma retry --id 45 --queue encounter
rxflow --config run_uma.toml uma resume --partition YOUR_GPU_PARTITION --workers 4
rxflow --config run_uma.toml uma verify-ts --dry-run
```

Stop existing workers before recovery or copying a live run. `resume` drains the
existing queues and preserves completed results; it does not reselect or import
anything from the original cluster. Logs, optimized XYZs, inputs, outputs,
frequency files, queues and provenance all remain inside the prepared folder.

## Return results and select UMA in the notebooks

After the destination workers finish:

```sh
rsync -a gpu-cluster:/scratch/project/Ni_TM_uma/ Ni_hpc_results/uma/
```

In either `reaction_diagrams.ipynb` or `system_energy_comparison.ipynb`, set:

```python
ENERGY_MODE = "uma"
```

The default lookup is the `uma` folder inside each original system result folder.
`UMA_RESULTS_ROOT` in the reaction notebook or `UMA_RESULTS_NAME` in the comparison
notebook can point to a different returned folder. The original DFT default and
plot layouts remain available. Each plot uses only the chosen method; incomplete
UMA states remain blank. The balanced comparisons, TS imaginary frequencies,
and `get_files(...)` use UMA's own selected structures when UMA is selected.
