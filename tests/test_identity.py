import unittest
from src.backend.inchi.determine_levels_id import InChI
from rdkit import Chem

class TestInChI(unittest.TestCase):
    def test_complete_identity(self):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertTrue(InChI.isCompleteIdentity(inchi, inchi))

    def test_are_equal_no_isotopes(self):
        inchi1="InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1"
        inchi2="InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/i1D"
        self.assertTrue(InChI.areEqualNoIsotopes(inchi1, inchi2))

    def test_are_equal_no_isotopes2(self):
        inchi1="InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/i1D"
        inchi2="InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/i9+1"
        self.assertTrue(InChI.areEqualNoIsotopes(inchi1, inchi2))

    def test_equal_salts_anion(self):
        inchi1 = "InChI=1S/C23H45NO4.ClH/c1-5-6-7-8-9-10-11-12-13-14-15-16-17-18-23(27)28-21(19-22(25)26)20-24(2,3)4;/h21H,5-20H2,1-4H3;1H"
        inchi2 = "InChI=1S/C23H45NO4.BrH/c1-5-6-7-8-9-10-11-12-13-14-15-16-17-18-23(27)28-21(19-22(25)26)20-24(2,3)4;/h21H,5-20H2,1-4H3;1H"
        self.assertTrue(InChI.areEqualDisolvedSalts(inchi1, inchi2))

    def test_equal_salts_cation(self):
        inchi1= "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        inchi2="InChI=1S/C2H4O2.K/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        self.assertTrue(InChI.areEqualDisolvedSalts(inchi1, inchi2))

    def test_equals_without_charges1(self):
        inchi1="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"
        inchi2="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChI.areEqualNoCharges(inchi1, inchi2))

    def test_equals_without_charges2(self):
        inchi1="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1"
        inchi2="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChI.areEqualNoCharges(inchi1, inchi2))
    
    def test_equals_without_charges3(self):
        inchi1="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1"
        inchi2="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p+1/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChI.areEqualNoCharges(inchi1, inchi2))

    def test_equals_without_charges4(self):
        inchi1="InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1"
        inchi2="InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3"
        self.assertTrue(InChI.areEqualNoCharges(inchi1, inchi2))

    def test_equal_no_position_double_bond(self):
        inchi1="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h9-10H,2-8,11-17H2,1H3,(H,19,20)/b10-9-"
        inchi2="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h8-9H,2-7,10-17H2,1H3,(H,19,20)/b9-8+"
        self.assertTrue(InChI.areEqualNoPositionDoubleBond(inchi1, inchi2))

    def test_equal_no_position_double_bond2(self):
        inchi1="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h8-9H,2-7,10-17H2,1H3,(H,19,20)/b9-8+"
        inchi2="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h7-8H,2-6,9-17H2,1H3,(H,19,20)/b8-7+"
        self.assertTrue(InChI.areEqualNoPositionDoubleBond(inchi1, inchi2))

    def test_equal_no_position_double_bond3(self):
        #acido cis-oleico
        inchi1="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h9-10H,2-8,11-17H2,1H3,(H,19,20)/b10-9-"
        #acido trans-oleico
        inchi2="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h9-10H,2-8,11-17H2,1H3,(H,19,20)/b10-9+"
        self.assertTrue(InChI.areEqualNoPositionDoubleBond(inchi1, inchi2))

    def test_are_equal_tautomers(self):
        #Tautomería oxoenólica
        inchi1="InChI=1S/C5H10O/c1-3-4-5(2)6/h3-4H2,1-2H3"
        inchi2="InChI=1S/C5H10O/c1-3-4-5(2)6/h4,6H,3H2,1-2H3/b5-4-"
        self.assertTrue(InChI.areEqualTautomers(inchi1, inchi2))
    
    def test_are_equal_tautomers2(self):
        #Tautomería imina-enamina
        inchi1="InChI=1S/C5H11N/c1-3-4-5(2)6/h6H,3-4H2,1-2H3"
        inchi2="InChI=1S/C5H11N/c1-3-4-5(2)6/h4H,3,6H2,1-2H3/b5-4-"
        self.assertTrue(InChI.areEqualTautomers(inchi1, inchi2))

    def test_are_equal_substituents(self):
        inchi1="InChI=1S/C26H46NO8P/c1-6-8-10-12-13-14-15-17-19-26(29)35-24(22-32-25(28)18-16-11-9-7-2)23-34-36(30,31)33-21-20-27(3,4)5/h6-9,11,16,24H,10,12-15,17-23H2,1-5H3/b8-6-,9-7+,16-11+/t24-/m1/s1"
        inchi2="InChI=1S/C26H46NO8P/c1-5-6-7-8-9-10-11-16-19-26(30)33-22-24(17-14-12-13-15-18-25(28)29)23-35-36(31,32)34-21-20-27(2,3)4/h5-6,12-15,24H,7-11,16-23H2,1-4H3,(H-,28,29,31,32)/b6-5-,14-12+,15-13+/t24-/m1/s1"
        self.assertTrue(InChI.areEqualSubstituentIndependent(inchi1, inchi2))


if __name__ == "__main__":
    unittest.main()
