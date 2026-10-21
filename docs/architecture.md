Architectural guide for building this workflow:

Runtime modules live in src/ and are installed as top-level modules through
pyproject.toml. Keep their import names stable for Slurm workers and notebooks.
Tests and fixtures live in tests/, batch templates in scripts/, and local
investigations in experiments/. Use uv sync --locked after layout changes.

INPUT:
	.csv file
	ligands + central metal + denticity + metal geometry
 	specific columns: ligand_smiles, metal, metal_geom, denticity, ID

STEP 1:
	embed ligand + central metal: sample all coordination modes -> Creates smiles
	Follow the methodology in reaction_smarts.ipynb:
	1. Build all central-metal to ligand permutations (denticity important here)
	1.2 Embed these using the metal geometry, just one conformer.
	2. Score using gxTB singlepoints (subprocess XTB_EXE, file, --gxtb --charge 1 --sp; skip solvent for now)
	3. Select those within 20 kcal/mol of the minimum
	4. Add the substrate in both the alpha and beta position to those within the lowest 20 kcal/mol
	5. Create product from reactant using the Chem.AllChem ReactionFromSmarts.
	6. Use Chem.ResonanceMolSupplier to create all resonance modes and keep the unique mols.

	Step 1 should give us the following:
	Atom mapped smiles for reactant and product - each ligand may have multiple configurational isomers from 1. ligand connectivity and 2. resonance structures.

STEP 2:
	Going from smiles -> Embedded minima endpoint structures and reactive conformations.
	1. Embed smiles using rxembed: request 20 conformers of each endpoint minimum, then prune by RMSD using prism_pruner.
        1.2a ---- Note: We will now be using ORCA for the following: there needs to be specific and ORCA_EXE.

	1.2 For the surviving conformers, perform gxTB optimizations using ORCA with extparams for the gxtb binary and the h2o solvent.
	1.3 Filter intended bond-length changes greater than 1.5 angstrom from their starting values, without requiring a perceived graph change. Separately reject native metal-coordination changes, including substrate dissociation after smaller movements. Select conformers within 3 kcal/mol of the lowest passing energy. Saved reaction endpoints can be filtered with filter-endpoints; reuse passing conformers and queue a verified original embedding only for states with no passing conformers.
	2. Embed smiles using rxembed.coembed: request 10 paired conformers and select the 5 lowest ts_rmsds pairs (or all surviving pairs if fewer than 5).
	2.1 Using rxts, make a TS guess from each selected pair. Keep the lowest-TS-RMSD coembedded endpoints and TS guess for NEB fallback.
	2.2a ---- Note: We will now be using ORCA for the following: there needs to be specific and ORCA_EXE.
	2.2 Preoptimize each guess with the mapped transfer atoms held, then run up to two direct OptTS attempts. If optimization or TS validation fails, score its saved rxts scan with gas-phase g-xTB, select the highest valid interior local maximum, preoptimize that guess, and try two more OptTS attempts.
	2.3 Use numerical Hessians and the TS mode with the right atom indexes.

STEP 3:
	Going from TS guess -> Full endpoint + TS reaction scheme.
	1. After creating the TS guess geometries, we want to use a TS optimizing routine.
	2. Run separate NumFreq on every converged OptTS; require exactly 1 imaginary frequency. Keep a database row for each candidate and retain the lowest Gibbs-energy successful TS per ligand/pair/seating.
	3. Only if all direct candidates fail, use NEB with the lowest-energy checked g-xTB interior peak and its matching coembedded endpoints, followed by separate NumFreq and imaginary-mode validation. Preserve all scans, guesses, scores and prior outcomes for comparison.
	4. Once the initial searches finish, reembed the lowest-G validated TS per ligand/pair/seating with rxembed.embed(..., n=20, fix=[donor, hydrogen, acceptor]). Prune by heavy-atom RMSD (prism_pruner, 0.5 Å default), then run the same constrained preoptimization, released OptTS and NumFreq/mode validation. Retain the lowest-G successful TS across original and sampled candidates. Sample each group only once, on workers, with persistent parent/seed/count metadata and retryable tasks.
	5. Apply RMSD-only prism pruning after every ordinary embedding, including precursors, references and encounters. Do not prune coembedded pairs. Fewer surviving conformers are normal, not missing-seed failures.


GUIDANCE: This should be a very minimal but clearly written code structure. Do not abstract functions too much. Use the code supplied by me as much as possible.
The regular submit/resume workflow schedules precursor and uncoordinated
substrate encounter preparation on Slurm, followed by precursor workers.
Encounter sampling covers each retained ligand–metal binding arrangement and
shares equivalent site/resonance routes. A dependent singlepoint submission job
selects minima after both reaction and precursor worker arrays terminate.
DFT selection takes the lowest-G minimum per state and all validated TS conformers
within 3 kcal/mol of the lowest g-xTB NumFreq Gibbs energy per ligand/pair/seating.
Both original and sampled TSs compete in that window. Each receives its own
thermal correction; the lowest corrected-G TS is then used in the diagrams.
Resume must perform no molecular embedding in the submitting terminal.
Encounter optimizations have no geometry filter or fragment-separation limit;
only optimization convergence and zero imaginary frequencies are required.
Preparation requeues historical encounter geometry rejections from their saved
converged XYZs at NumFreq when possible, with an archived database/row snapshot.
This is the first of a larger module, so make it modular - i.e. i'd like calls like create_permutations and run_gxtb and so on.
Use descriptive names such as generate_reaction_pairs, sample_reaction_endpoints and generate_ts_guess.
Write results under results/ligand_01/ using the CSV ID, with reactant/, product/ and ts_guess/ folders.
XYZ, INP and log files go directly in these folders; identify pairs, seatings and conformers in filenames.

GUIDANCE FILES:
Projects/workflow/FURTHER_TESTS: src and scripts contains graph connectivity functions and structural code for inspiration - we do not yet want to create a database, but write the code to be portable to a sqlite3 database structure like that
Projects/rxembed: Contains the embedding and coembedding code.
Projects/rxts: Contains the TS guess creating code.
Projects/drug-degradation: Contains reaction network structural code for inspiration. I would like to eventually make a .toml where I can give multiple sequential smarts reactions / steps that this workflow can tie together and find transition states between, like a reaction A->B->C->D and so on.
Projects/inspiration: For the ORCA approach and orchestration and solvent with ORCA use.

RUN SETTINGS: Run all ligand workflows in parallel, one worker and one CPU thread per ligand. Request 20 conformers per reaction endpoint, 20 conformers for TS sampling, and 10 coembedded reactive pairs.

