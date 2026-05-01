import unittest
from src.backend.inchi.inchi_parser import InChIParser
from rdkit import Chem

class TestInChI(unittest.TestCase):
    def test_getMainLayer(self):
        self.assertEqual(InChIParser.getMainLayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"), "C6H6")
        self.assertEqual(InChIParser.getMainLayer("InChI=1S/C2H6O/i1+1/h1H2"), "C2H6O")
        self.assertEqual(InChIParser.getMainLayer("InChI=1/C3H6O/c1-3(2)4/h1-2H3/f/h1H"), "C3H6O")

    def test_getAtomConnectionsSublayer(self):
        self.assertEqual(InChIParser.getAtomConnectionsSublayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"), "c1-2-4-6-5-3-1")
        self.assertTrue(InChIParser.getAtomConnectionsSublayer("InChI=1S/C2H6O.Na/c1-2-3;/h3H,2H2,1H3;/q;+1/rC2H6O.Na/c1-2-3;/h3H,2H2,1H3;/q;+1").startswith("c"))
        self.assertEqual(InChIParser.getAtomConnectionsSublayer("InChI=1S/C2H6O.Na/c1-2-3;/h3H,2H2,1H3;/q;+1/rC2H6O.Na/c1-2-3;/h3H,2H2,1H3;/q;+1"), "c1-2-3;")

    def test_getHydrogenAtomsSublayer(self):
        self.assertEqual(InChIParser.getHydrogenAtomsSublayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"), "h1-6H")
        self.assertEqual(InChIParser.getHydrogenAtomsSublayer("InChI=1S/C2H6O/i1+1/h1H2"), "h1H2")
        self.assertEqual(InChIParser.getHydrogenAtomsSublayer("InChI=1S/C2H6O.Na/c1-2-3;/h3H,2H2,1H3;/q;+1/rC2H6O.Na/c1-2-3;/h3H,2H2,1H3;/q;+1"), "h3H,2H2,1H3;")

    def test_getChargeSublayer(self):
        self.assertEqual(InChIParser.getChargeSublayer("InChI=1S/H3O+/h1H2/q+1"), "q+1")
        self.assertIsNone(InChIParser.getChargeSublayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"))

    def test_getProtonSublayer(self):
        self.assertEqual(InChIParser.getProtonSublayer("InChI=1S/C2H4O/q-1/p+1"), "p+1")
        self.assertIsNone(InChIParser.getProtonSublayer("InChI=1S/H3O+/h1H2/q+1"))

    def test_getDoubleBondsSublayer(self):
        self.assertEqual(InChIParser.getDoubleBondsSublayer("InChI=1S/C2H4/b1-2/h1-2H2"), "b1-2")
        self.assertEqual(InChIParser.getDoubleBondsSublayer("InChI=1S/C3H8/i2+2/b2-3/t2-/m2-/s1"), "b2-3")

    def test_getTetrahedralStereoSublayer(self):
        self.assertEqual(InChIParser.getTetrahedralStereoSublayer("InChI=1S/C4H10O/t3-/h3H"), "t3-")
        self.assertEqual(InChIParser.getTetrahedralStereoSublayer("InChI=1S/C3H8/i2+2/b2-3/t2-/m2-/s1"), "t2-")

    def test_getTypeStereoInfoSublayer(self):
        self.assertEqual(InChIParser.getTypeStereoInfoSublayer("InChI=1S/C3H8/i2+2/b2-3/t2-/m2-/s1"), "s1")
        self.assertIsNone(InChIParser.getTypeStereoInfoSublayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"))

    def test_getIsotopicLayer(self):
        self.assertEqual(InChIParser.getIsotopicLayer("InChI=1S/C2H6O/i1+1/h1H2"), "i1+1")
        self.assertEqual(InChIParser.getIsotopicLayer("InChI=1S/C3H8/i2+2/b2-3/t2-/m2-/s1"), "i2+2")
        self.assertIsNone(InChIParser.getIsotopicLayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"))

    def test_getIsotopicHydrogenSublayer(self):
        self.assertEqual(InChIParser.getIsotopicHydrogenSublayer("InChI=1S/CH4/i/h1H"), "h1H")
        self.assertEqual(InChIParser.getIsotopicHydrogenSublayer("InChI=1S/C2H6O/i1+1/h1H2"), "h1H2")

    def test_getIsotopicStereoSublayer(self):
        self.assertEqual(InChIParser.getIsotopicStereoSublayer("InChI=1S/C3H8/i2+2/b2-3/t2-/m2-/s1"), "b2-3/t2-/m2-/s1")
        self.assertEqual(InChIParser.getIsotopicStereoSublayer("InChI=1S/C3H6/i1/t1-"), "t1-")
        self.assertIsNone(InChIParser.getIsotopicStereoSublayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"))
    
    def test_getFixedHLayer(self):
        self.assertEqual(InChIParser.getFixedHLayer("InChI=1/C3H6O/c1-3(2)4/h1-2H3/f/h1H"), "f")
        self.assertIsNone(InChIParser.getFixedHLayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"))

    def test_getReconnectedLayer(self):
        self.assertTrue(InChIParser.getReconnectedLayer("InChI=1S/C2H6O.Na/c1-2-3;/h3H,2H2,1H3;/q;+1/rC2H6O.Na/c1-2-3;/h3H,2H2,1H3;/q;+1").startswith("r"))
        self.assertIsNone(InChIParser.getReconnectedLayer("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"))

    """def test_carboxylate_neutralization_connectivity(self):
        inchi = "InChI=1S/C2H3O2/c1-2(3)4/h1H3/q-1"

        mol = Chem.MolFromInchi(inchi)
        self.assertIsNotNone(mol)

        neutral = InChIParser.neutralize_molecule(mol)
        neutral_inchi = Chem.MolToInchi(neutral)

        self.assertEqual(
            InChIParser.getAtomConnectionsSublayer(inchi),
            InChIParser.getAtomConnectionsSublayer(neutral_inchi),
        )"""

if __name__ == "__main__":
    unittest.main()