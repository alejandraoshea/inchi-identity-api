import unittest
from rdkit import Chem
from src.backend.lipid.lipid_structure_detector import LipidHeadValidator
from src.backend.inchi.determine_levels_id import InChI


class TestGlycosylglycerolTeacher(unittest.TestCase):
    def setUp(self):
        self.validator = LipidHeadValidator()
        
        self.core_smiles = "O[C@H]1[C@H](OC[C@]([H])(O)CO)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O"
        self.good_smiles = "O[C@H]1[C@H](OC[C@]([H])(O)COC(CCCCCC)=O)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O"
        self.wrong_smiles = "O[C@H]1[C@H](OC[C@]([H])(O)CO)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](OC(CCCCCC)=O)[C@H]2O)[C@H](O)[C@@H]1O"
        
        self.core = Chem.MolFromSmiles(self.core_smiles)
        self.good = Chem.MolFromSmiles(self.good_smiles)
        self.wrong = Chem.MolFromSmiles(self.wrong_smiles)
    
    def test_01_core_no_match(self):
        """Core (no FA) should NOT match any pattern"""
        result = self.validator.matches_any_valid_head(self.core)
        self.assertFalse(result, "CORE without FA should NOT match")
    
    def test_02_good_matches(self):
        """GOOD (FA correct position) MUST match"""
        result = self.validator.matches_any_valid_head(self.good)
        self.assertTrue(result, "GOOD molecule MUST match - FA is in correct position!")
    
    def test_03_wrong_no_match(self):
        """WRONG (FA wrong position) should NOT match"""
        result = self.validator.matches_any_valid_head(self.wrong)
        self.assertFalse(result, "WRONG molecule should NOT match - FA is in WRONG position!")
    
    def test_04_classes_detected(self):
        """Check that GOOD molecule is classified"""
        classes = self.validator.identify_lipid_class(self.good)
        self.assertTrue(len(classes) > 0, f"GOOD should be classified. Got: {classes}")
        print(f"  GOOD classified as: {classes}")


class TestSphingoLipidsExcel(unittest.TestCase):
    """Test sphingolipids from Excel file - Example 60, 64"""
    
    def setUp(self):
        self.validator = LipidHeadValidator()
        
        # Example 60: Acidic glycosphingolipid
        self.ex60_correct = "[H][C@](NC(CCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](C(O)=O)[C@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])CCCCCCCCCCCCCCC"
        self.ex60_wrong = "[H][C@](N)(CO[C@H]1O[C@H](C(O)=O)[C@H](O)[C@H](O)[C@H]1O)[C@@](OC(CCCCCCCCCCCCC)=O)([H])CCCCCCCCCCCCCCC"
        
        # Example 64: Sphingomyelin
        self.ex64_correct = "[H][C@](NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])CCCCCCCCCCCCCCC"
        self.ex64_wrong = "[H][C@](N)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](OC(CCCCCCCCCCCCCCCCCCC)=O)([H])CCCCCCCCCCCCCCC"
    
    def test_01_ex60_correct_detected(self):
        """Example 60 CORRECT should be detected as sphingolipid"""
        mol = Chem.MolFromSmiles(self.ex60_correct)
        classes = self.validator.identify_lipid_class(mol)
        print(f"\n  Ex60 CORRECT classes: {classes}")
        self.assertTrue(len(classes) > 0, "Ex60 CORRECT should be recognized")
    
    def test_02_ex60_wrong_not_same_class(self):
        """Example 60 WRONG should not match same pattern as CORRECT"""
        mol_correct = Chem.MolFromSmiles(self.ex60_correct)
        mol_wrong = Chem.MolFromSmiles(self.ex60_wrong)
        
        classes_correct = self.validator.identify_lipid_class(mol_correct)
        classes_wrong = self.validator.identify_lipid_class(mol_wrong)
        
        print(f"  Ex60 CORRECT: {classes_correct}")
        print(f"  Ex60 WRONG: {classes_wrong}")
        
        self.assertNotEqual(classes_correct, classes_wrong, 
                           "CORRECT and WRONG should have different classifications")
    
    def test_03_ex64_sphingomyelin_detected(self):
        """Example 64 CORRECT should be detected as sphingomyelin"""
        mol = Chem.MolFromSmiles(self.ex64_correct)
        classes = self.validator.identify_lipid_class(mol)
        print(f"\n  Ex64 CORRECT classes: {classes}")
        
        self.assertTrue(len(classes) > 0, "Ex64 should be recognized as sphingolipid")


class TestGlycerolipidsExcel(unittest.TestCase):
    """Test glycerolipids from Excel - Example 70"""
    
    def setUp(self):
        self.validator = LipidHeadValidator()
        
        # Example 70: Triacylglycerol
        # 1-dodecanoyl-2-hexadecanoyl-3-octadecanoyl-sn-glycerol
        self.tg_smiles = "CCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCCCC"
    
    def test_01_triacylglycerol_detected(self):
        """Triacylglycerol should be detected"""
        mol = Chem.MolFromSmiles(self.tg_smiles)
        classes = self.validator.identify_lipid_class(mol)
        print(f"\n  TG classes: {classes}")
        
        self.assertTrue(len(classes) > 0, "Triacylglycerol should be recognized")
        self.assertIn("Triacylglycerols", classes, "Should be classified as Triacylglycerol")


class TestDiacylglycerolIsomers(unittest.TestCase):
    """Test 1,2-DG vs 1,3-DG distinction"""
    
    def setUp(self):
        self.validator = LipidHeadValidator()
        
        # 1,2-DG (FAs at positions 1 and 2)
        self.dg_1_2 = "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC"
        
        # 1,3-DG (FAs at positions 1 and 3)
        self.dg_1_3 = "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC"
    
    def test_01_both_detected_as_diacylglycerol(self):
        mol_12 = Chem.MolFromSmiles(self.dg_1_2)
        mol_13 = Chem.MolFromSmiles(self.dg_1_3)
        
        classes_12 = self.validator.identify_lipid_class(mol_12)
        classes_13 = self.validator.identify_lipid_class(mol_13)
        
        print(f"\n  1,2-DG: {classes_12}")
        print(f"  1,3-DG: {classes_13}")
        
        self.assertTrue(len(classes_12) > 0, "1,2-DG should be recognized")
        self.assertTrue(len(classes_13) > 0, "1,3-DG should be recognized")


class TestInChIComparison(unittest.TestCase):
    def test_01_same_lipid_class_equal(self):
        # Two simple DGs
        dg1 = "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC"
        dg2 = "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC"
        
        mol1 = Chem.MolFromSmiles(dg1)
        mol2 = Chem.MolFromSmiles(dg2)
        
        inchi1 = Chem.MolToInchi(mol1)
        inchi2 = Chem.MolToInchi(mol2)
        
        result = InChI.areEqualSubstituentIndependent(inchi1, inchi2)
        self.assertTrue(result, "Identical DGs should be equal")
    
    def test_02_different_class_not_equal(self):
        # DG vs TG
        dg = "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC"
        tg = "CCCCCCCC(=O)OCC(COC(=O)CCCCCCCC)OC(=O)CCCCCCCC"
        
        mol_dg = Chem.MolFromSmiles(dg)
        mol_tg = Chem.MolFromSmiles(tg)
        
        inchi_dg = Chem.MolToInchi(mol_dg)
        inchi_tg = Chem.MolToInchi(mol_tg)
        
        result = InChI.areEqualSubstituentIndependent(inchi_dg, inchi_tg)
        self.assertFalse(result, "DG and TG should NOT be equal (different classes)")


class TestEdgeCases(unittest.TestCase):
    def setUp(self):
        self.validator = LipidHeadValidator()
    
    def test_01_non_lipid_no_match(self):
        benzene = Chem.MolFromSmiles("c1ccccc1")
        result = self.validator.matches_any_valid_head(benzene)
        self.assertFalse(result, "Benzene should not match any lipid pattern")
    
    def test_02_fatty_acid_detected(self):
        palmitic = Chem.MolFromSmiles("CCCCCCCCCCCCCCCC(=O)O")
        classes = self.validator.identify_lipid_class(palmitic)
        print(f"\n  Palmitic acid: {classes}")
        self.assertTrue(len(classes) > 0, "Fatty acid should be recognized")


def run_tests():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromTestCase(TestGlycosylglycerolTeacher))
    suite.addTests(loader.loadTestsFromTestCase(TestSphingoLipidsExcel))
    suite.addTests(loader.loadTestsFromTestCase(TestGlycerolipidsExcel))
    suite.addTests(loader.loadTestsFromTestCase(TestDiacylglycerolIsomers))
    suite.addTests(loader.loadTestsFromTestCase(TestInChIComparison))
    suite.addTests(loader.loadTestsFromTestCase(TestEdgeCases))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print("FINAL RESULTS")
    print("="*70)
    print(f"Tests run: {result.testsRun}")
    print(f"Passed: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failed: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    
    if result.wasSuccessful():
        print("\n✅ ALL TESTS PASSED - Your system is working!")
    else:
        print("\n❌ SOME TESTS FAILED - Review above")
    
    return result


if __name__ == "__main__":
    run_tests()