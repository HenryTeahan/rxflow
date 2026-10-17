Architectural guide for building this workflow:

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
	1. Embed smiles using rxembed: we want to embed 3 conformers of each endpoints minima. 
        1.2a ---- Note: We will now be using ORCA for the following: there needs to be specific and ORCA_EXE.

	1.2 For these 3 conformers, we want to perform gxTB optimizations: We now will be using an orca calculation - using extparams for the gxtb binary, and extparam for the h2o solvent.
	1.3 Select the conformers within 3 kcal/mol of the lowest. Only discard an optimization for connectivity if the graph changes AND the associated atom-pair distance changes by more than 2 angstrom from its starting value.
	2. Embed smiles using rxembed.coembed: request 10 paired conformers and select the pair with the lowest ts_rmsds value.
	2.1 Using rxts, create a transitionstateguesser object, and make one TS guess from that selected pair.
	2.2a ---- Note: We will now be using ORCA for the following: there needs to be specific and ORCA_EXE.
	2.2 Constrained relax this guess, keeping the TS atoms constrained to their TS guess positions.
	2.3 For each constrained relaxed TS guess, we want to write an input file, containing the guess structure, numhess = true, calcfreq=true and the TS mode with the right atom indexes.

STEP 3:
	Going from TS guess -> Full endpoint + TS reaction scheme.
	1. After creating the TS guess geometries, we want to use a TS optimizing routine.
	2. In inspiration/workflow_orca.py, a draft for this is theorized. First try optts until it converges or we give up -> alternatively, NEB fallback.
	3. Finally, once converged, we want to run a numfreq. Only success criteria is 1 imaginary frequency.


GUIDANCE: This should be a very minimal but clearly written code structure. Do not abstract functions too much. Use the code supplied by me as much as possible.
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

RUN SETTINGS: Run all ligand workflows in parallel, one worker and one CPU thread per ligand. Use 3 endpoint optimization conformers and 10 coembedded reactive pairs.

