"""Check substrate-site enumeration and mapped H transfer for the three new CSVs."""

from collections import Counter
from pathlib import Path
import unittest

from rdkit import Chem, rdBase

from rxflow import add_substrate, create_permutations, create_products, read_ligands
from rxflow_endpoints import hydrogen_transfer, mapped_molecule

FIXTURES = Path(__file__).resolve().parent / 'fixtures'
CASES = {
    "ED_dF_ligands.csv": dict(aromatic_h=3, other_h=3, sites={"site_008", "site_009"}),
    "ED_pyr_ligands.csv": dict(aromatic_h=4, other_h=3, sites={"site_006", "site_007"}),
    "Ni_ligands.csv": dict(aromatic_h=3, other_h=9, sites={"site_001", "site_002"}),
}


def aromatic_h_atoms(molecule):
    return [atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 1
            and atom.GetNeighbors()[0].GetAtomicNum() == 6
            and atom.GetNeighbors()[0].GetIsAromatic()]


def marked_site_signature(substrate, carbon_index):
    """Compare whole marked graphs independently of the production symmetry ranks."""
    marked = Chem.Mol(substrate)
    for atom in marked.GetAtoms():
        atom.SetAtomMapNum(0)
    marked.GetAtomWithIdx(carbon_index).SetIsotope(13)
    return Chem.MolToSmiles(marked, canonical=True)


def hydrogen_bonds(molecule):
    return {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
            for bond in molecule.GetBonds()
            if 1 in (bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum())}


class SubstrateHydrogenTests(unittest.TestCase):
    def test_all_csv_rows_have_the_expected_hydrogens_and_symmetry_classes(self):
        for filename, expected in CASES.items():
            rows = read_ligands(FIXTURES / filename)
            self.assertEqual(len(rows), 13)
            for row in rows:
                with self.subTest(csv=filename, ligand=row["ID"]):
                    substrate = Chem.AddHs(Chem.MolFromSmiles(row["substrate"]))
                    aromatic_h = aromatic_h_atoms(substrate)
                    other_h = [a for a in substrate.GetAtoms() if a.GetAtomicNum() == 1
                               and a.GetIdx() not in {h.GetIdx() for h in aromatic_h}]
                    self.assertEqual(len(aromatic_h), expected["aromatic_h"])
                    self.assertEqual(len(other_h), expected["other_h"])
                    signatures = {marked_site_signature(substrate, h.GetNeighbors()[0].GetIdx())
                                  for h in aromatic_h}
                    self.assertEqual(len(signatures), 2)

    def test_sites_and_h_transfer_across_every_ligand_donor_subset(self):
        totals = Counter()
        # Some chemically invalid resonance products are rejected inside create_products.
        with rdBase.BlockLogs():
            for filename, expected in CASES.items():
                for row in read_ligands(FIXTURES / filename):
                    substrate = Chem.AddHs(Chem.MolFromSmiles(row["substrate"]))
                    expected_signatures = {
                        marked_site_signature(substrate, h.GetNeighbors()[0].GetIdx())
                        for h in aromatic_h_atoms(substrate)
                    }
                    for donors, complex_mol in create_permutations(
                            row["ligand_smiles"], row["metal"], row["denticity"]):
                        with self.subTest(csv=filename, ligand=row["ID"], donors=donors):
                            totals["coordination_graphs"] += 1
                            offset = complex_mol.GetNumAtoms()
                            pairs = list(add_substrate(complex_mol, row["substrate"]))
                            self.assertEqual({site for site, _ in pairs}, expected["sites"])
                            signatures = set()
                            for site, reactant in pairs:
                                totals["substrate_sites"] += 1
                                substrate_index = int(site.removeprefix("site_")) - 1
                                signatures.add(marked_site_signature(substrate, substrate_index))
                                carbon = offset + substrate_index
                                donor_atom = reactant.GetAtomWithIdx(carbon)
                                self.assertEqual(donor_atom.GetAtomicNum(), 6)
                                self.assertTrue(donor_atom.GetIsAromatic())
                                self.assertTrue(any(a.GetAtomicNum() == 1 for a in donor_atom.GetNeighbors()))
                                coordination = reactant.GetBondBetweenAtoms(carbon, offset - 1)
                                self.assertIsNotNone(coordination)
                                self.assertEqual(coordination.GetBondType(), Chem.BondType.DATIVE)
                                products = list(create_products(reactant, row["metal"]))
                                self.assertTrue(products, "No valid H-transfer product for this substrate site")
                                for resonance, expected_acceptor, product in products:
                                    totals["mapped_h_transfers"] += 1
                                    self.assertEqual([a.GetAtomicNum() for a in reactant.GetAtoms()],
                                                     [a.GetAtomicNum() for a in product.GetAtoms()])
                                    indices = hydrogen_transfer(reactant, product)
                                    donor, hydrogen, acceptor = indices
                                    self.assertEqual(donor, carbon)
                                    self.assertEqual(acceptor, expected_acceptor)
                                    self.assertEqual(reactant.GetAtomWithIdx(hydrogen).GetAtomicNum(), 1)
                                    self.assertGreaterEqual(hydrogen, offset)
                                    self.assertLess(acceptor, offset - 1)
                                    self.assertIn(reactant.GetAtomWithIdx(acceptor).GetAtomicNum(), (7, 8, 16))
                                    self.assertEqual(hydrogen_bonds(reactant) - hydrogen_bonds(product),
                                                     {tuple(sorted((donor, hydrogen)))})
                                    self.assertEqual(hydrogen_bonds(product) - hydrogen_bonds(reactant),
                                                     {tuple(sorted((acceptor, hydrogen)))})
                                    # Exercise the same mapped-SMILES restoration used for HPC inputs.
                                    restored_r = mapped_molecule(Chem.MolToSmiles(reactant))
                                    restored_p = mapped_molecule(Chem.MolToSmiles(product))
                                    self.assertEqual(hydrogen_transfer(restored_r, restored_p), indices)
                            self.assertEqual(signatures, expected_signatures)
        print("Graph-level hydrogen checks:", dict(totals))

    def test_input_atom_maps_do_not_create_extra_substrate_sites(self):
        for filename in CASES:
            row = read_ligands(FIXTURES / filename)[0]
            _, complex_mol = next(create_permutations(row["ligand_smiles"], row["metal"], row["denticity"]))
            original = Chem.MolFromSmiles(row["substrate"])
            for atom in original.GetAtoms():
                atom.SetAtomMapNum(100 + atom.GetIdx())
            tagged_smiles = Chem.MolToSmiles(original)
            tagged_substrate = Chem.AddHs(Chem.MolFromSmiles(tagged_smiles))
            offset = complex_mol.GetNumAtoms()
            pairs = list(add_substrate(complex_mol, tagged_smiles))
            self.assertEqual(len(pairs), 2)
            signatures = {marked_site_signature(tagged_substrate, int(site.removeprefix("site_")) - 1)
                          for site, _ in pairs}
            expected = {marked_site_signature(tagged_substrate, h.GetNeighbors()[0].GetIdx())
                        for h in aromatic_h_atoms(tagged_substrate)}
            self.assertEqual(signatures, expected)
            for _, reactant in pairs:
                self.assertEqual([a.GetAtomMapNum() for a in reactant.GetAtoms()],
                                 list(range(1, reactant.GetNumAtoms() + 1)))


if __name__ == "__main__":
    unittest.main()
