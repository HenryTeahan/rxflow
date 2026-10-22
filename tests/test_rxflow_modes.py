import unittest

from rdkit import Chem

from rxflow import retain_coordination_mode


class RetainCoordinationModeTests(unittest.TestCase):
    def setUp(self):
        self.ligand = Chem.AddHs(Chem.MolFromSmiles(
            "[O-]NC(=O)COC1=NC(=CC=C1)C1=CC=CC=N1"))

    def test_keeps_tridentate_nitrogen_mode_above_general_cutoff(self):
        self.assertTrue(retain_coordination_mode([1, 7, 17], self.ligand, 23.54))

    def test_other_high_energy_modes_remain_pruned(self):
        self.assertFalse(retain_coordination_mode([0, 7, 17], self.ligand, 37.59))

    def test_general_cutoff_and_missing_energy_behavior_are_unchanged(self):
        self.assertTrue(retain_coordination_mode([0, 7, 17], self.ligand, 19.99))
        self.assertTrue(retain_coordination_mode([0, 7, 17], self.ligand, None))


if __name__ == "__main__":
    unittest.main()
