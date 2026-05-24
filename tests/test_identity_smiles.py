import unittest
from src.backend.inchi.determine_levels_id import InChI
from rdkit import Chem


class TestNormalizeInput(unittest.TestCase):
    def test_valid_smiles_becomes_inchi(self):
        result = InChI.normalize_input("CCO")
        self.assertTrue(result.startswith("InChI="), f"Expected InChI, got: {result}")

    def test_inchi_passthrough(self):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertEqual(InChI.normalize_input(inchi), inchi)

    def test_invalid_smiles_returns_input(self):
        result = InChI.normalize_input("NOT_A_SMILES!!!")
        self.assertEqual(result, "NOT_A_SMILES!!!")

    def test_empty_string(self):
        self.assertEqual(InChI.normalize_input(""), "")

    def test_none_passthrough(self):
        self.assertIsNone(InChI.normalize_input(None))

    def test_canonical_smiles_equivalents(self):
        self.assertEqual(InChI.normalize_input("CCO"), InChI.normalize_input("OCC"))


class TestCompleteIdentitySmiles(unittest.TestCase):
    def test_same_smiles(self):
        self.assertTrue(InChI.isCompleteIdentity("CCO", "CCO"))

    def test_different_smiles_same_molecule(self):
        self.assertTrue(InChI.isCompleteIdentity("CCO", "OCC"))

    def test_different_molecules(self):
        self.assertFalse(InChI.isCompleteIdentity("CCO", "CCCO"))

    def test_smiles_vs_inchi_same_molecule(self):
        smi   = "CCO"
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertTrue(InChI.isCompleteIdentity(smi, inchi))

    def test_smiles_vs_inchi_different_molecule(self):
        smi   = "CCCO"  # propanol
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"  # ethanol
        self.assertFalse(InChI.isCompleteIdentity(smi, inchi))


class TestNoIsotopesSmiles(unittest.TestCase):
    def test_smiles_matches_isotope_inchi(self):
        inchi_with_isotope = "InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/i1D"
        smi_no_isotope = "N[C@@H](Cc1ccccc1)C(=O)O"
        self.assertTrue(InChI.areEqualNoIsotopes(smi_no_isotope, inchi_with_isotope))

    def test_two_smiles_same_molecule(self):
        self.assertTrue(InChI.areEqualNoIsotopes("CCO", "OCC"))

    def test_different_molecules(self):
        self.assertFalse(InChI.areEqualNoIsotopes("CCO", "CCCO"))


class TestDissolvedSaltsSmiles(unittest.TestCase):
    def test_same_smiles(self):
        self.assertTrue(InChI.areEqualDisolvedSalts("CCO", "CCO"))

    def test_different_halide_counter_ions(self):
        smiles_hcl = "C[N+](C)(C)CC(=O)O.Cl"
        smiles_hbr = "C[N+](C)(C)CC(=O)O.Br"
        self.assertTrue(InChI.areEqualDisolvedSalts(smiles_hcl, smiles_hbr))

    def test_different_metal_cations(self):
        na_acetate = "CC(=O)[O-].[Na+]"
        k_acetate  = "CC(=O)[O-].[K+]"
        self.assertTrue(InChI.areEqualDisolvedSalts(na_acetate, k_acetate))

    def test_smiles_vs_inchi_salt(self):
        na_smi   = "CC(=O)[O-].[Na+]"
        k_inchi  = "InChI=1S/C2H4O2.K/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1"
        self.assertTrue(InChI.areEqualDisolvedSalts(na_smi, k_inchi))

    def test_different_organic_fragments(self):
        self.assertFalse(InChI.areEqualDisolvedSalts("CCO", "CCCO"))


class TestNoChargesSmiles(unittest.TestCase):
    def test_carboxylate_vs_carboxylic_acid(self):
        self.assertTrue(InChI.areEqualNoCharges("CC(=O)[O-]", "CC(=O)O"))

    def test_zwitterion_vs_protonated(self):
        zwitterion  = "C[N+](C)(C)CC(=O)[O-]"
        protonated  = "C[N+](C)(C)CC(=O)O"
        self.assertTrue(InChI.areEqualNoCharges(zwitterion, protonated))

    def test_different_molecules(self):
        self.assertFalse(InChI.areEqualNoCharges("CC(=O)O", "CCC(=O)O"))

    def test_smiles_vs_charged_inchi(self):
        inchi_charged = "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-1/t4-,6-,7-,10-/m1/s1"
        inchi_neutral  = "InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"
        self.assertTrue(InChI.areEqualNoCharges(inchi_charged, inchi_neutral))


class TestNoStereoSmiles(unittest.TestCase):
    def test_l_and_d_phenylalanine(self):
        l_phe = "N[C@@H](Cc1ccccc1)C(=O)O"
        d_phe = "N[C@H](Cc1ccccc1)C(=O)O"
        self.assertTrue(InChI.areEqualNoStereo(l_phe, d_phe))

    def test_cis_trans_2_butene(self):
        cis   = r"C/C=C\C"
        trans = "C/C=C/C"
        self.assertTrue(InChI.areEqualNoStereo(cis, trans))

    def test_different_molecules(self):
        self.assertFalse(InChI.areEqualNoStereo("CCO", "CCCO"))

    def test_same_smiles(self):
        self.assertTrue(InChI.areEqualNoStereo("CCO", "CCO"))

    def test_same_molecule_different_smiles_no_stereo(self):
        self.assertTrue(InChI.areEqualNoStereo("CCO", "OCC"))


class TestNoPositionDoubleBondSmiles(unittest.TestCase):
    def test_cis_trans_oleic_smiles(self):
        cis_oleic   = "CCCCCCCCC=CCCCCCCCC(=O)O"
        trans_oleic = "CCCCCCCC/C=C/CCCCCCCC(=O)O"
        self.assertTrue(InChI.areEqualNoPositionDoubleBond(cis_oleic, trans_oleic))

    def test_shifted_double_bond(self):
        delta9  = "CCCCCCCCC=CCCCCCCCC(=O)O"
        delta11 = "CCCCCCC=CCCCCCCCCCC(=O)O"
        self.assertTrue(InChI.areEqualNoPositionDoubleBond(delta9, delta11))

    def test_different_chain_length(self):
        c18 = "CCCCCCCCC=CCCCCCCCC(=O)O"
        c16 = "CCCCCCC=CCCCCCCCC(=O)O"
        self.assertFalse(InChI.areEqualNoPositionDoubleBond(c18, c16))


class TestMixedInputTypes(unittest.TestCase):

    ETHANOL_SMILES = "CCO"
    ETHANOL_INCHI  = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"

    def test_complete_identity_mixed(self):
        self.assertTrue(InChI.isCompleteIdentity(self.ETHANOL_SMILES, self.ETHANOL_INCHI))

    def test_no_isotopes_mixed(self):
        self.assertTrue(InChI.areEqualNoIsotopes(self.ETHANOL_SMILES, self.ETHANOL_INCHI))

    def test_dissolved_salts_mixed(self):
        self.assertTrue(InChI.areEqualDisolvedSalts(self.ETHANOL_SMILES, self.ETHANOL_INCHI))

    def test_no_charges_mixed(self):
        self.assertTrue(InChI.areEqualNoCharges(self.ETHANOL_SMILES, self.ETHANOL_INCHI))

    def test_no_stereo_mixed(self):
        self.assertTrue(InChI.areEqualNoStereo(self.ETHANOL_SMILES, self.ETHANOL_INCHI))


if __name__ == "__main__":
    unittest.main()