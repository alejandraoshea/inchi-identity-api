import unittest
from inchi.determine_levels_id import InChi

class TestInChI(unittest.TestCase):
    def test_complete_identity(self):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertTrue(InChi.isCompleteIdentity(inchi, inchi))

    def test_are_equal_diluted_salts(self):
        inchi_a = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        inchi_b = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertTrue(InChi.areEqualDisolvedSalts(inchi_a, inchi_b))

    def test_equal_salts_anion(self):
        inchi1 = "InChI=1S/C23H45NO4.ClH/c1-5-6-7-8-9-10-11-12-13-14-15-16-17-18-23(27)28-21(19-22(25)26)20-24(2,3)4;/h21H,5-20H2,1-4H3;1H"
        inchi2 = "InChI=1S/C23H45NO4.BrH/c1-5-6-7-8-9-10-11-12-13-14-15-16-17-18-23(27)28-21(19-22(25)26)20-24(2,3)4;/h21H,5-20H2,1-4H3;1H"
        self.assertTrue(InChi.areEqualDisolvedSalts(inchi1, inchi2))

    def test_equal_salts_cation(self):
        inchi1= "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        inchi2="InChI=1S/C2H4O2.K/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        self.assertTrue(InChi.areEqualDisolvedSalts(inchi1, inchi2))


    """def test_equals_without_charges(self):
        inchi1="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"
        inchi2="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        inchi3="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1"
        inchi4="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        inchi5="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p+1/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi2))
        self.assertTrue(InChi.areEqualNoCharges(inchi3, inchi4))
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi5))
        #CHECK after having modified the code
        """
    
    #TODO: ask why this should be the same!!
    """
    def test_equals_without_charges1(self):
        inchi1="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"
        inchi2="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi2))
        #CHECK

    def test_equals_without_charges2(self):
        inchi1="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1"
        inchi2="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi2))
    def test_equals_without_charges3(self):
        inchi1="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1"
        inchi2="InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p+1/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi2))"""
 
    def test_equals_without_charges4(self):
        inchi1="InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1"
        inchi2="InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3 "
        self.assertTrue(InChi.areEqualNoCharges(inchi1, inchi2))

    def test_equal_no_position_double_bond_position(self):
        inchi1="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h9-10H,2-8,11-17H2,1H3,(H,19,20)/b10-9-"
        inchi2="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h8-9H,2-7,10-17H2,1H3,(H,19,20)/b9-8+"
        inchi3="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h7-8H,2-6,9-17H2,1H3,(H,19,20)/b8-7+"
        self.assertTrue(InChi.areEqualNoPositionDoubleBond(inchi1, inchi2))
        self.assertTrue(InChi.areEqualNoPositionDoubleBond(inchi3, inchi2))

    def test_equal_no_position_double_bond(self):
        inchi1="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h9-10H,2-8,11-17H2,1H3,(H,19,20)/b10-9-"
        inchi2="InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h9-10H,2-8,11-17H2,1H3,(H,19,20)/b10-9+"
        self.assertTrue(InChi.areEqualNoPositionDoubleBond(inchi1, inchi2))

    def test_are_equal_tautomers(self):
        inchi1="InChI=1S/C5H10O/c1-3-4-5(2)6/h3-4H2,1-2H3"
        inchi2="InChI=1S/C5H10O/c1-3-4-5(2)6/h4,6H,3H2,1-2H3/b5-4-"
        self.assertTrue(InChi.areEqualTautomers(inchi1, inchi2))

    def test_are_equal_no_isotopes(self):
        inchi1="InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1"
        inchi2="InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/i1D"
        inchi3="InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/i9+1"
        self.assertTrue(InChi.areEqualNoIsotopes(inchi1, inchi2))
        self.assertTrue(InChi.areEqualNoIsotopes(inchi3, inchi2))

if __name__ == "__main__":
    unittest.main()
