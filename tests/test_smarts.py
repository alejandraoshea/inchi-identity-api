import unittest
from rdkit import Chem
from backend.lipid.lipid_structure_detector import LipidHeadValidator


class TestLipidHeadValidator(unittest.TestCase):

    def setUp(self):
        self.validator = LipidHeadValidator()

        # Glycosylglycerol examples
        self.core = Chem.MolFromSmiles(
            "O[C@H]1[C@H](OC[C@]([H])(O)CO)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O"
        )
        self.good = Chem.MolFromSmiles(
            "O[C@H]1[C@H](OC[C@]([H])(O)COC(CCCCCC)=O)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O"
        )
        self.wrong = Chem.MolFromSmiles(
            "O[C@H]1[C@H](OC[C@]([H])(O)CO)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](OC(CCCCCC)=O)[C@H]2O)[C@H](O)[C@@H]1O"
        )

    def test_glycosylglycerol_core(self):
        self.assertFalse(self.validator.matches_any_valid_head(self.core))

    def test_glycosylglycerol_good(self):
        self.assertTrue(self.validator.matches_any_valid_head(self.good))

    def test_glycosylglycerol_wrong(self):
        self.assertFalse(self.validator.matches_any_valid_head(self.wrong))

    def test_detailed_validation(self):
        result = self.validator.validate_structure(self.good, verbose=True)
        self.assertIsInstance(result, dict)
        self.assertTrue(result.get("is_valid", False))

    def test_identify_lipid_classes(self):
        test_cases = {
            "1,2-DG": "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC",
            "1,3-DG": "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC",
            "TG": "CCCCCCCC(=O)OCC(COC(=O)CCCCCCCC)OC(=O)CCCCCCCC",
        }

        for name, smiles in test_cases.items():
            with self.subTest(name=name):
                mol = Chem.MolFromSmiles(smiles)
                lipid_class = self.validator.identify_lipid_class(mol)
                self.assertIsNotNone(lipid_class)


if __name__ == "__main__":
    unittest.main()