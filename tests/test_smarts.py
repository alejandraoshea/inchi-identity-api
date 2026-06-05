import unittest
from rdkit import Chem
from src.backend.lipid.lipid_analysis import LipidHeadValidator, LipidAnalysis
from src.backend.inchi.inchi_parser import InChIParser

class TestEx7__2_NeutralGlycosphingolipids(unittest.TestCase):
    """Galβ1-4GlcNAcα1-3Galβ1-4GlcNAcβ1-3Galβ1-4Glcβ-Cer(d18:1/16:0)"""
    NEG_CHAIN = "[H][C@](N)(CO[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O[C@H]3[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]4[C@H](O)[C@@H](O[C@@H]5[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]6[C@H](O)[C@@H](O)[C@@H](O)[C@H](O6)CO)[C@H](O5)CO)[C@@H](O)[C@H](O4)CO)[C@H](O3)CO)[C@@H](O)[C@H](O2)CO)[C@H](O1)CO)[C@@](OC(CCCCCCCCCCCCCCC)=O)([H])/C=C/CCCCCCCCCCCCC"
    VAR_SUGAR = "[H][C@](NC(CCCCCCCCCCCCCCC)=O)(CO[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O[C@H]3[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]4[C@H](O)[C@@H](O[C@H]5[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]6[C@H](O)[C@@H](O[C@H]7[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]8[C@H](O)[C@@H](O[C@@H]9[C@H](O)[C@@H](O)[C@@H](O)[C@H](O9)CO)[C@@H](O)[C@H](O8)CO)[C@H](O7)CO)[C@@H](O)[C@H](O6)CO[C@H]%10[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]%11[C@H](O)[C@@H](O[C@@H]%12[C@H](O)[C@@H](O)[C@@H](O)[C@H](O%12)CO)[C@@H](O)[C@H](O%11)CO)[C@H](O%10)CO)[C@H](O5)CO)[C@@H](O)[C@H](O4)CO[C@H]%13[C@H](NC(C)=O)[C@@H](O)[C@H](O)[C@H](O%13)CO)[C@H](O3)CO)[C@@H](O)[C@H](O2)CO)[C@H](O1)CO)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx60__2_AcidicGlycosphingolipids(unittest.TestCase):
    """N-(tetradecanoyl)-1-α-galacturonosyl-sphinganine"""
    NEG_CHAIN = "[H][C@](N)(CO[C@H]1O[C@H](C(O)=O)[C@H](O)[C@H](O)[C@H]1O)[C@@](OC(CCCCCCCCCCCCC)=O)([H])CCCCCCCCCCCCCCC"
    VAR_SUGAR = "[H][C@](NC(CCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](C(O)=O)[C@@H](O[C@H]2O[C@H](CO[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H](O)[C@H]2N)[C@H](O)[C@H]1O)[C@@](O)([H])CCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx66__2_Phosphosphingolipids(unittest.TestCase):
    """N-(2-hydroxyhexacosanoyl)-4R-hydroxyeicosasphinganine-1-phospho-(1'-my"""
    NEG_CHAIN = "[H][C@](N)(COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)[C@@](OC(C(O)CCCCCCCCCCCCCCCCCCCCCCCC)=O)([H])CCCCCCCCCCCCCCCCC"
    VAR_SUGAR = "[H][C@](NC(C(O)CCCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O)(O)=O)[C@@](O)([H])CCCCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx67__2_AcidicGlycosphingolipids(unittest.TestCase):
    """N-(15Z-tetracosenoyl)-1-β-(3'-sulfo)-glucosyl-sphing-4-enine"""
    VAR_SUGAR = "[H][C@](NC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)(CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx69__2_Unknown(unittest.TestCase):
    """N-(hexacosanoyl)-4R-hydroxysphinganine-1-phospho-(1'-[2-amino-2-deoxy-"""
    NEG_CHAIN = "CCCCCCCCCCCCCC[C@@]([H])(O)[C@](OC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O)([H])[C@](N)([H])COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2N)(O)=O"
    VAR_SUGAR = "CCCCCCCCCCCCCC[C@@]([H])(O)[C@](O)([H])[C@](NC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O)([H])COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@H](C)[C@@H](O)[C@H](O)[C@H]2N)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx3__3_BetaineLipids(unittest.TestCase):
    """1-hexadecanoyl-2-(6Z,9Z,12Z-octadecatrienoyl)-sn-glycero-3-O-(N,N,N-tr"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COCCC(C([O-])=O)[N+](C)(C)C)OC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COCCC(C([O-])=O)[N+](C)(C)C)OC(CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx76__2_OtherSphingolipids(unittest.TestCase):
    """N-(13-methyltetradecanoyl)-15-methylhexadecasphinganine-1-phospho-β-D-"""
    NEG_CHAIN = "NC(COP(O)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)=O)C(OC(CCCCCCCCCCCC(C)C)=O)CCCCCCCCCCCC(C)C"
    VAR_SUGAR = "OP(OCC(NC(CCCCCCCCCCCC(C)C)=O)C(O)CCCCCCCCCCCC(C)C)(O[C@@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx77__2_Ceramides(unittest.TestCase):
    """(2S,3S,26R,27S)-2,27-diamino-3,26,28-trihydroxyoctacosan-11-one"""
    POS = "OC[C@H](N)[C@H](O)CCCCCCCCCCCCCCC(CCCCCCC[C@H](O)[C@@H](N)C)=O"
    NEG_STEREO = "OC[C@H](N)[C@@H](O)CCCCCCCCCCCCCCC(CCCCCCC[C@H](O)[C@@H](N)C)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx78__2_Ceramides(unittest.TestCase):
    """N-(9Z-octadecenoyl)-4E,6E-tetradecasphingadienine"""
    POS = "[H][C@](NC(CCCCCCC/C=C\\CCCCCCCC)=O)(CO)[C@@](O)([H])/C=C/C=C/CCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx80__2_BasicSpNeutralAndAcidic(unittest.TestCase):
    """1-β-galactosyl-sphing-4-enine"""
    POS = "[H][C@](O)([C@](N)([H])CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)/C=C/CCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx81__02_AmphotericGlycosphingolipids(unittest.TestCase):
    """2S,3R,4E)-2-amino-3-hydroxyoctadec-4-en-1-yl β-D-galactopyranoside 6-("""
    POS = "O[C@H]([C@H]1O)[C@H](OC[C@H](N)[C@H](O)/C=C/CCCCCCCCCCCCC)O[C@H](COS(O)(=O)=O)[C@@H]1O"
    NEG_STEREO = "O[C@H]([C@H]1O)[C@H](OC[C@H](N)[C@@H](O)/C=C/CCCCCCCCCCCCC)O[C@H](COS(O)(=O)=O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx3__3_Sterols(unittest.TestCase):
    """campest-5-en-3β-yl octadecanoate"""
    POS = "C[C@]12CC[C@H](OC(CCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](C)C(C)C"
    NEG_CHAIN = "O[C@H](C1)CC[C@@]2(C)C1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](CC(CCCCCCCCCCCCCCCCC)=O)([H])CC[C@@H](C)C(C)C"
    NEG_STEREO = "C[C@]12CC[C@@H](OC(CCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](C)C(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx4__3_Sterols(unittest.TestCase):
    """Stigmast-5-en-3β-yl octadecanoate"""
    POS = "C[C@]12CC[C@H](OC(CCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](CC)C(C)C"
    NEG_CHAIN = "O[C@H](C1)CC[C@@]2(C)C1=C(C(CCCCCCCCCCCCCCCCC)=O)C[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](CC)C(C)C"
    NEG_STEREO = "C[C@]12CC[C@@H](OC(CCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])CC[C@@H](CC)C(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx5__3_Sterols(unittest.TestCase):
    """Stigmast-5,22E-dien-3β-yl docosanoate"""
    POS = "C[C@]12CC[C@H](OC(CCCCCCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])/C=C/[C@@H](CC)C(C)C"
    NEG_CHAIN = "O[C@H](C1)CC[C@@]2(C)C1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])/C=C/[C@@H](CC)C(C(CCCCCCCCCCCCCCCCCCCCC)=O)C"
    NEG_STEREO = "C[C@]12CC[C@@H](OC(CCCCCCCCCCCCCCCCCCCCC)=O)CC1=CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@@]3([H])CC[C@]4([H])[C@@](C)([H])/C=C/[C@@H](CC)C(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx6__3_Sterols(unittest.TestCase):
    """6-oxo-5α-campestan-3β,22R,23R-triol 3β-yl tetradecanoate"""
    POS = "[H][C@@]12[C@]([C@](CC[C@H](OC(CCCCCCCCCCCCC)=O)C3)(C)[C@@]3([H])C(C2)=O)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)[C@@H](O)[C@H](O)[C@@H](C)C(C)C"
    NEG_CHAIN = "[H][C@@]12[C@]([C@](CC[C@H](O)C3)(C)[C@@]3([H])C(C2)=O)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)[C@@H](OC(CCCCCCCCCCCCC)=O)[C@H](O)[C@@H](C)C(C)C"
    NEG_STEREO = "[H][C@@]12[C@]([C@](CC[C@@H](OC(CCCCCCCCCCCCC)=O)C3)(C)[C@@]3([H])C(C2)=O)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)[C@@H](O)[C@H](O)[C@@H](C)C(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx7__3_Sterols(unittest.TestCase):
    """3β,5β,10β,14β-tetrahydroxy-19-norbufa-20,22-dienolide-3β-yl-14-hydroxy"""
    POS = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54O"
    NEG_CHAIN = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)CC[C@@]54O"
    NEG_STEREO = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx8__3_Sterols(unittest.TestCase):
    """3β,5β,14β-trihydroxy-19oxo-bufa-20,22-dienolide-3β-yl-14-hydroxy-tetra"""
    POS = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54C=O"
    NEG_CHAIN = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)CC[C@@]54C=O"
    NEG_STEREO = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54C=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx9__3_Sterols(unittest.TestCase):
    """3β,5β,14β-trihydroxy-bufa-20,22-dienolide-3β-yl-14-hydroxy-tetradecano"""
    POS = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54C"
    NEG_CHAIN = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)CC[C@@]54C"
    NEG_STEREO = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@H](OC(CCCCCCCCCCCCCO)=O)CC[C@@]54C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx10__3_Sterols(unittest.TestCase):
    """1β,3β,5β,14β-tetrahydroxy-bufa-20,22-dienolide-3β-yl-14-hydroxy-tetrad"""
    POS = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@@H](OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)[C@@]54C"
    NEG_CHAIN = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)C[C@@H](O)[C@@]54C"
    NEG_STEREO = "C[C@@]1([C@]2(O)CC[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5(O)C[C@H](OC(CCCCCCCCCCCCCO)=O)C[C@@H](O)[C@@]54C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx11__3_Sterols(unittest.TestCase):
    """24-methylene,26,27-dimethylcholest-5-en-3β-yl 18-bromooctadeca-9E,17E-"""
    POS = "C[C@@]1([C@@]2([H])CC[C@]1([H])[C@@](C)([H])CCC(C(CC)CC)=C)CC[C@@]3([H])[C@@]2([H])CC=C4C[C@@H](OC(CCCC#CC#C/C=C/CCCCC#C/C=C/Br)=O)CC[C@@]43C"
    NEG_CHAIN = "C[C@@]1([C@@]2([H])CC[C@]1([H])[C@@](CC(CCCC#CC#C/C=C/CCCCC#C/C=C/Br)=O)([H])CCC(C(CC)CC)=C)CC[C@@]3([H])[C@@]2([H])CC=C4C[C@@H](O)CC[C@@]43C"
    NEG_STEREO = "C[C@@]1([C@@]2([H])CC[C@]1([H])[C@@](C)([H])CCC(C(CC)CC)=C)CC[C@@]3([H])[C@@]2([H])CC=C4C[C@H](OC(CCCC#CC#C/C=C/CCCCC#C/C=C/Br)=O)CC[C@@]43C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx12__3_Sterols(unittest.TestCase):
    """3β,5β-dihydroxy-14β,15β-epoxy-bufa-20,22-dienolide-3β-yl-16-hydroxy-9Z"""
    POS = "O=C(O[C@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@]35[C@H](O5)C[C@]4([H])C(C=C6)=COC6=O)(O)C1)CCCCCCC/C=C\\CCCCCCO"
    NEG_CHAIN = "C[C@@]1([C@]23[C@H](O3)C[C@]1([H])C(C=C4)=COC4=O)CC[C@@]5([H])[C@@]2([H])CC[C@]6(OC(CCCCCCC/C=C\\CCCCCCO)=O)C[C@@H](O)CC[C@@]65C"
    NEG_STEREO = "O=C(O[C@@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@]35[C@H](O5)C[C@]4([H])C(C=C6)=COC6=O)(O)C1)CCCCCCC/C=C\\CCCCCCO"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx13__3_Sterols(unittest.TestCase):
    """12-oxo-3β,11α,14β-trihydroxy-5β-bufa-20,22-dienolide-3β-yl-14-hydroxy-"""
    POS = "O=C(O[C@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])[C@H](O)C([C@@]4(C)[C@]3(O)CC[C@]4([H])C(C=C5)=COC5=O)=O)([H])C1)CCCCCCCCCCCCCO"
    NEG_CHAIN = "C[C@@]1([C@]2(OC(CCCCCCCCCCCCCO)=O)CC[C@]1([H])C(C=C3)=COC3=O)C([C@@H](O)[C@@]4([H])[C@@]2([H])CC[C@]5([H])C[C@@H](O)CC[C@@]54C)=O"
    NEG_STEREO = "O=C(O[C@@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])[C@H](O)C([C@@]4(C)[C@]3(O)CC[C@]4([H])C(C=C5)=COC5=O)=O)([H])C1)CCCCCCCCCCCCCO"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx14__3_Sterols(unittest.TestCase):
    """3β,14β-dihydroxy-16β-acetoxy-5β-bufa-20,22-dienolide-3β-yl-hexadecanoa"""
    POS = "C[C@@]1([C@]2(O)C[C@H](OC(C)=O)[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5([H])C[C@@H](OC(CCCCCCCCCCCCCCC)=O)CC[C@@]54C"
    NEG_CHAIN = "O[C@H]1CC[C@@]2(C)[C@@](CC[C@]3([H])[C@]2([H])CC[C@@]4(C)[C@]3(OC(CCCCCCCCCCCCCCC)=O)C[C@H](OC(C)=O)[C@]4([H])C(C=C5)=COC5=O)([H])C1"
    NEG_STEREO = "C[C@@]1([C@]2(O)C[C@H](OC(C)=O)[C@]1([H])C(C=C3)=COC3=O)CC[C@@]4([H])[C@@]2([H])CC[C@]5([H])C[C@H](OC(CCCCCCCCCCCCCCC)=O)CC[C@@]54C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx15__3_Isoprenoids(unittest.TestCase):
    """epi-Austrobuxusin H"""
    POS = "O=C(O1)[C@]([C@H](C)OC(CCCCCCCCC)=O)([H])C[C@@]1([C@H]2[C@@H]3O2)[C@@]4(C)[C@]3(O)[C@@H]5[C@H](C(C)=C)[C@@H](OC5=O)[C@H]4O"
    NEG_CHAIN = "O=C(O1)[C@]([C@H](C)O)([H])C[C@@]1([C@H]2[C@@H]3O2)[C@@]4(C)[C@]3(OC(CCCCCCCCC)=O)[C@@H]5[C@H](C(C)=C)[C@@H](OC5=O)[C@H]4O"
    NEG_STEREO = "O=C(O1)[C@]([C@H](C)OC(CCCCCCCCC)=O)([H])C[C@@]1([C@H]2[C@@H]3O2)[C@@]4(C)[C@@]3(O)[C@@H]5[C@H](C(C)=C)[C@@H](OC5=O)[C@H]4O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx16__3_Isoprenoids(unittest.TestCase):
    """Ustusoic acid B"""
    POS = "C[C@@]12[C@]([C@H](OC(/C=C/C=C/C=C/C)=O)C=C(C=O)[C@]2(O)C(O)=O)([H])C(C)(C)CCC1"
    NEG_CHAIN = "C[C@@]12[C@]([C@H](O)C=C(C=O)[C@]2(O)C(OC(/C=C/C=C/C=C/C)=O)=O)([H])C(C)(C)CCC1"
    NEG_STEREO = "C[C@@]12[C@]([C@@H](OC(/C=C/C=C/C=C/C)=O)C=C(C=O)[C@]2(O)C(O)=O)([H])C(C)(C)CCC1"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx17__3_Isoprenoids(unittest.TestCase):
    """1-O-[(2E,4E,6E,8E)-3,7-dimethyl-9-(2,6,6-trimethylcyclohex-1-en-1-yl)n"""
    POS = "CC(/C=C/C1=C(C)CCCC1(C)C)=C\\C=C\\C(C)=C\\C(O[C@@H]2O[C@H](C(O)=O)[C@@H](O)[C@H](O)[C@H]2O)=O"
    VAR_SUGAR = "CC(/C=C/C1=C(C)CCCC1(C)C)=C\\C=C\\C(C)=C\\C(O[C@@H]2O[C@@H](C)[C@@H](O)[C@H](O)[C@H]2O)=O"
    NEG_STEREO = "CC(/C=C/C1=C(C)CCCC1(C)C)=C\\C=C\\C(C)=C\\C(O[C@H]2O[C@H](C(O)=O)[C@@H](O)[C@H](O)[C@H]2O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx18__3_Isoprenoids(unittest.TestCase):
    """3,7-dimethyl-9-(2,6,6-trimethylcyclohex-1-en-1-yl)nona-2E,4Z,6E,8E-tet"""
    POS = "CC(/C=C\\C=C(C)\\C=C\\C1=C(CCCC1(C)C)C)=C\\COC(CCCCCCCCCCCCCCC)=O"
    NEG_CHAIN = "CC(/C=C\\C=C(C)\\C=C\\C1=C(C(CCCCCCCCCCCCCCC)=O)CCCC1(C)C)=C\\CO.C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))


class TestEx19__3_QuinonesAndHydroquinones(unittest.TestCase):
    """N-(5-[(2R)-6-Hydroxy-2,5,7,8-tetramethyl-3,4-dihydrochromen-2-yl]-prop"""
    POS = "CC1=C(C)C(O)=C(C)C2=C1O[C@](C)(CCC(NCC(O)=O)=O)CC2"
    NEG_CHAIN = "CC1=C(C)C(OCC(O)=O)=C(C)C2=C1O[C@](C)(CCC(N)=O)CC2"
    NEG_STEREO = "CC1=C(C)C(O)=C(C)C2=C1O[C@@](C)(CCC(NCC(O)=O)=O)CC2"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx20__3_Acylaminosugars(unittest.TestCase):
    """UDP-3-(3R-hydroxy-tetradecanoyl)-αD-glucosamine"""
    POS = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@H]3CO)O)OC(C[C@@H](CCCCCCCCCCC)O)=O)N)=O"
    NEG_CHAIN = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](OC(C[C@@H](CCCCCCCCCCC)O)=O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@H]3CO)O)O)N)=O"
    NEG_STEREO = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@@H]([C@@H]([C@H]3CO)O)OC(C[C@@H](CCCCCCCCCCC)O)=O)N)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx21__3_Acylaminosugars(unittest.TestCase):
    """2,3-bis-(3R-hydroxy-tetradecanoyl)-αD-glucosamine-1-phosphate"""
    POS = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(O)=O)O[C@@H]1CO"
    NEG_CHAIN = "N[C@H]1[C@@H](OP(O)(OC(C[C@H](O)CCCCCCCCCCC)=O)=O)O[C@H](COC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(O)=O)O[C@H]1CO"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx22__3_Acylaminosugars(unittest.TestCase):
    """Lipid IVA"""
    POS = "OP(O[C@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]2NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](O)[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]1NC(C[C@H](O)CCCCCCCCCCC)=O)(O)=O"
    NEG_CHAIN = "N[C@H]1[C@@H](OP(OC(C[C@H](O)CCCCCCCCCCC)=O)(OC(C[C@H](O)CCCCCCCCCCC)=O)=O)O[C@H](CO[C@@H]2O[C@H](COC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(OC(C[C@H](O)CCCCCCCCCCC)=O)=O)[C@H](O)[C@H]2N)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "OP(O[C@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]2NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](O)[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]1NC(C[C@H](O)CCCCCCCCCCC)=O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx23__3_Acylaminosugars(unittest.TestCase):
    """UDP-2,3-bis-(3R-hydroxy-tetradecanoyl)-αD-glucosamine"""
    POS = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@H]3CO)O)OC(C[C@@H](CCCCCCCCCCC)O)=O)NC(C[C@@H](CCCCCCCCCCC)O)=O)=O"
    NEG_CHAIN = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1OC(C[C@@H](CCCCCCCCCCC)O)=O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@H]3COC(C[C@@H](CCCCCCCCCCC)O)=O)O)O)N)=O"
    NEG_STEREO = "OP(OP(OC[C@H]1O[C@@H](N(C(N2)=O)C=CC2=O)[C@H](O)[C@@H]1O)(O)=O)(O[C@@H](O3)[C@@H]([C@H]([C@@H]([C@@H]3CO)O)OC(C[C@@H](CCCCCCCCCCC)O)=O)NC(C[C@@H](CCCCCCCCCCC)O)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx24__3_Acylaminosugars(unittest.TestCase):
    """Lipid A -disaccharide-1-phosphate"""
    POS = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(O)=O)O[C@@H]1CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]2NC(C[C@H](O)CCCCCCCCCCC)=O"
    NEG_CHAIN = "N[C@H]1[C@@H](OP(OC(C[C@H](O)CCCCCCCCCCC)=O)(OC(C[C@H](O)CCCCCCCCCCC)=O)=O)O[C@H](CO[C@@H]2O[C@H](COC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H](O)[C@H]2N)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](OP(O)(O)=O)O[C@@H]1CO[C@@H]2O[C@H](CO)[C@@H](O)[C@@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@H]2NC(C[C@H](O)CCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx25__3_AcylaminosugarGlycans(unittest.TestCase):
    """Dodecanoyl-Kdo2-Lipid IVA"""
    POS = "O[C@H]([C@@H](CO[C@@H]([C@@H]1NC(C[C@H](OC(CCCCCCCCCCC)=O)CCCCCCCCCCC)=O)O[C@H](CO[C@]2(C(O)=O)C[C@@H](O[C@@]3(C(O)=O)O[C@H]([C@H](O)CO)[C@H](O)[C@H](O)C3)[C@@H](O)[C@@H]([C@H](O)CO)O2)[C@@H](OP(O)(O)=O)[C@@H]1OC(C[C@H](O)CCCCCCCCCCC)=O)O[C@H](OP(O)(O)=O)[C@@H]4NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H]4OC(C[C@H](O)CCCCCCCCCCC)=O"
    NEG_CHAIN = "O[C@H]1[C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H](CO[C@@H]([C@@H]2N)O[C@H](CO[C@]3(C(O)=O)C[C@@H](O[C@@]4(C(O)=O)O[C@H]([C@H](O)COC(C[C@H](O)CCCCCCCCCCC)=O)[C@H](O)[C@H](O)C4)[C@@H](O)[C@@H]([C@H](OC(C[C@H](O)CCCCCCCCCCC)=O)COC(C[C@H](OC(CCCCCCCCCCC)=O)CCCCCCCCCCC)=O)O3)[C@@H](OP(O)(O)=O)[C@@H]2O)O[C@H](OP(O)(O)=O)[C@@H]1N"
    NEG_STEREO = "O[C@H]([C@@H](CO[C@@H]([C@@H]1NC(C[C@H](OC(CCCCCCCCCCC)=O)CCCCCCCCCCC)=O)O[C@H](CO[C@]2(C(O)=O)C[C@@H](O[C@@]3(C(O)=O)O[C@H]([C@H](O)CO)[C@H](O)[C@H](O)C3)[C@@H](O)[C@@H]([C@H](O)CO)O2)[C@@H](OP(O)(O)=O)[C@@H]1OC(C[C@H](O)CCCCCCCCCCC)=O)O[C@@H](OP(O)(O)=O)[C@@H]4NC(C[C@H](O)CCCCCCCCCCC)=O)[C@@H]4OC(C[C@H](O)CCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx26__3_Acyltrehaloses(unittest.TestCase):
    """2-O-hexadecanoyl-3-O-(2R,4S,6S-trimethyl-3R-hydroxy-tricosanoyl)-α,α-t"""
    POS = "O[C@H]([C@@H](CO)O[C@H](O[C@@H](O1)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](OC([C@H](C)[C@H](O)[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]1CO)[C@@H]2O)[C@@H]2O"
    NEG_CHAIN = "O[C@H]1[C@H](OC[C@H]([C@H](O)[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCC)C=O)[C@@H](CO)O[C@H](O[C@@H](O2)[C@H](O)[C@@H](O)[C@H](O)[C@H]2COC(CCCCCCCCCCCCCCC)=O)[C@@H]1O"
    NEG_STEREO = "O[C@H]([C@@H](CO)O[C@H](O[C@@H](O1)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](OC([C@H](C)[C@H](O)[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]1CO)[C@@H]2O)[C@H]2O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx27__3_Acyltrehaloses(unittest.TestCase):
    """2-O-hexadecanoyl-3-O-(2,4S,6S-trimethyl-2E-docosenoyl)-6-O-(2,4S,6S-tr"""
    POS = "O[C@H]1[C@H](OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCCCC)=O)[C@@H](CO)O[C@H](O[C@@H](O2)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]2COC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCC)=O)[C@@H]1OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "O[C@H]1[C@H](OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCCCC)=O)[C@H](CO)O[C@H](O[C@@H](O2)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]2COC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCC)=O)[C@@H]1OC(/C(C)=C/[C@@H](C)C[C@@H](C)CCCCCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx28__3_Acyltrehaloses(unittest.TestCase):
    """2-O-hexadecanoyl,3-O-(2S,4S,6S,8R,10R,12R,14R-heptamethyl)-15-hydroxy-"""
    POS = "O=C(O[C@H]1[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(=O)(O)=O)O[C@H](CO)[C@@H](O)[C@@H]1OC([C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C(O)CCCCCCCCCCCCCCC)=O)CCCCCCCCCCCCCCC"
    NEG_STEREO = "O=C(O[C@H]1[C@@H](O[C@H]2O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]2OS(=O)(O)=O)O[C@H](CO)[C@@H](O)[C@@H]1OC([C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C[C@@H](C)C(O)CCCCCCCCCCCCCCC)=O)CCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx29__3_Acyltrehaloses(unittest.TestCase):
    """6-O-hexadecanoyl-α-D-glucopyranosyl 6-O-hexadecanoyl-α-D-glucopyranosi"""
    POS = "O[C@H]1[C@@H](O[C@H]2O[C@H](COC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@H](O)[C@H]2O)O[C@H](COC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@H]1[C@@H](O[C@H]2O[C@H](COC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@@H](O)[C@H]2O)O[C@H](COC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx30__3__29_Acyltrehaloses(unittest.TestCase):
    """Acyltrehaloses"""
    POS = "O[C@H]1[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@H]2O)O[C@H](CO)[C@@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx31__3_Acyltrehaloses(unittest.TestCase):
    """Acyltrehaloses"""
    POS = "O[C@H]1[C@@H](O[C@H]2O[C@H](COC(CCCCCCCCCCCCCCCCC)=O)[C@@H](OC(CCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]2O)O[C@H](CO)[C@@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@H]1[C@@H](O[C@H]2O[C@H](COC(CCCCCCCCCCCCCCCCC)=O)[C@@H](OC(CCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H]2O)O[C@H](CO)[C@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx32__3_OtherAcylSugars(unittest.TestCase):
    """2,3-di-O-hexanoyl-α-glucopyranose"""
    POS = "O[C@H]1O[C@H](CO)[C@@H](O)[C@H](OC(CCCCC)=O)[C@H]1OC(CCCCC)=O"
    NEG_STEREO = "O[C@H]1O[C@H](CO)[C@@H](O)[C@@H](OC(CCCCC)=O)[C@H]1OC(CCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx33__3_OtherAcylSugars(unittest.TestCase):
    """Butyl 4'-O-hexadecanoyl-neohesperidoside"""
    POS = "O[C@@H]1[C@@H](O[C@@H]2O[C@@H](C)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@H]2O)[C@H](OCCCC)O[C@H](CO)[C@H]1O"
    NEG_STEREO = "O[C@@H]1[C@@H](O[C@H]2O[C@@H](C)[C@H](OC(CCCCCCCCCCCCCCC)=O)[C@@H](O)[C@H]2O)[C@H](OCCCC)O[C@H](CO)[C@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx34__3_OtherAcylSugars(unittest.TestCase):
    """Butyl 3'-O-acetyl-2'-O-butanoyl-4,6,4'-tri-O-(2-methylpropanoyl)-neohe"""
    POS = "O[C@@H]1[C@@H](O[C@@H]2O[C@@H](C)[C@H](OC(C(C)C)=O)[C@@H](OC(C)=O)[C@H]2OC(CCC)=O)[C@H](OCCCC)O[C@H](COC(C(C)C)=O)[C@H]1OC(C(C)C)=O"
    NEG_STEREO = "O[C@H]1[C@@H](O[C@@H]2O[C@@H](C)[C@H](OC(C(C)C)=O)[C@@H](OC(C)=O)[C@H]2OC(CCC)=O)[C@H](OCCCC)O[C@H](COC(C(C)C)=O)[C@H]1OC(C(C)C)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx35__3_OtherAcylSugars(unittest.TestCase):
    """(11S)-jalapinolic 11-O-α-l-rhamnopyranosyl-(1-3)-O-[4-O-n-dodecanoyl-α"""
    POS = "CCCCC[C@H](O[C@@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O2)CCCCCCCCCC(O[C@H]3[C@@H](O)[C@H]2O[C@@H](C)[C@@H]3O[C@@H]4O[C@@H](C)[C@H](O[C@@H]5O[C@@H](C)[C@H](OC(CCCCCCC)=O)[C@@H](O)[C@H]5O)[C@@H](O[C@@H]6O[C@@H](C)[C@H](O)[C@@H](O)[C@H]6O)[C@H]4OC(CCCCCCCCCCC)=O)=O"
    NEG_STEREO = "CCCCC[C@H](O[C@@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O2)CCCCCCCCCC(O[C@H]3[C@@H](O)[C@H]2O[C@@H](C)[C@@H]3O[C@H]4O[C@@H](C)[C@H](O[C@@H]5O[C@@H](C)[C@H](OC(CCCCCCC)=O)[C@@H](O)[C@H]5O)[C@@H](O[C@@H]6O[C@@H](C)[C@H](O)[C@@H](O)[C@H]6O)[C@H]4OC(CCCCCCCCCCC)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx36__3_OtherAcylSugars(unittest.TestCase):
    """(6R,9S)-9-O-β-D-glucopyranosyloxy-6'-O- ([''Z,1''Z,1''Z]-triene)-octad"""
    POS = "O=C(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)OC[C@H]1O[C@@H](O[C@@H](C)/C=C/[C@@]2(O)C(C)(C)CC(C=C2C)=O)[C@H](O)[C@@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx37__3_OtherAcylSugars(unittest.TestCase):
    """2-O-(6-O-(9Z,12Z,15Z-octadecatrienoyl)-α-D-glucopyranosyl)-β-D-fructof"""
    POS = "O[C@@H]1[C@@](O[C@H]([C@@H]([C@H]2O)O)O[C@H](COC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)[C@H]2O)(CO)O[C@H](CO)[C@H]1O"
    NEG_STEREO = "O[C@@H]1[C@@](O[C@H]([C@@H]([C@@H]2O)O)O[C@H](COC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)[C@H]2O)(CO)O[C@H](CO)[C@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx38__3_OtherAcylSugars(unittest.TestCase):
    """6-O-isopropoyl-6'-O-(14-methyl-2Z,4E-hexadecadienoyl)-maltose"""
    POS = "O[C@H]1[C@@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2COC(C(C)C)=O)O[C@H](COC(/C=C\\C=C\\CCCCCCCCC(CC)C)=O)[C@@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx39__3_OtherAcylSugars(unittest.TestCase):
    """(2-O-octadecanoyl-3-O-isobutyroyl)-2R-(α-D-glucopyranosyloxy)-3-hydrox"""
    POS = "OC[C@H](C(O)=O)O[C@@H](O1)[C@@H]([C@H]([C@@H]([C@H]1CO)O)OC(C(C)C)=O)OC(CCCCCCCCCCCCCCCCC)=O"
    NEG_CHAIN = "OC[C@H](C(OC(C(C)C)=O)=O)O[C@@H](O1)[C@@H]([C@H]([C@@H]([C@H]1COC(CCCCCCCCCCCCCCCCC)=O)O)O)O"
    NEG_STEREO = "OC[C@H](C(O)=O)O[C@@H](O1)[C@@H]([C@H]([C@H]([C@H]1CO)O)OC(C(C)C)=O)OC(CCCCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx40__3_OtherGlycerophospholipids(unittest.TestCase):
    """1,-2-diacyl-sn-glycero-3-phospho-O-[N-(2-hydroxyethyl)glycine]"""
    POS = "O=C(C)OC[C@]([H])(COP(OCCNCC(O)=O)(O)=O)OC(C)=O"
    NEG_STEREO = "O=C(C)OCC([H])(COP(OCCNCC(O)=O)(O)=O)OC(C)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx41__3_OtherGlycerophospholipids(unittest.TestCase):
    """1,2-dihexadecanoyl-sn-glycero-3-phosphosulfocholine"""
    POS = "[H][C@](OC(CCCCCCCCCCCCCCC)=O)(COP(OCC[S+](C)C)([O-])=O)COC(CCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "[H]C(OC(CCCCCCCCCCCCCCC)=O)(COP(OCC[S+](C)C)([O-])=O)COC(CCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx42__3_OtherGlycerophospholipids(unittest.TestCase):
    """PE 16:0/18:1(9Z)-15-isoLG pyrrole"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCCN1C(C)=C(C/C=C\\CCCC(O)=O)C(/C=C/C(O)CCCCC)=C1)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OCCN1C(C)=C(C/C=C\\CCCC(O)=O)C(/C=C/C(OC(CCCCCCCCCCCCCCC)=O)CCCCC)=C1)(OC(CCCCCCC/C=C\\CCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OCCN1C(C)=C(C/C=C\\CCCC(O)=O)C(/C=C/C(O)CCCCC)=C1)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx42_OH__3_OtherGlycerophospholipids(unittest.TestCase):
    """PE 16:0/18:1(9Z)-15-isoLG hydroxylactam"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCCN(C1=O)C(C)(O)C(C/C=C\\CCCC(O)=O)=C1/C=C/C(O)CCCCC)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx43__3_OtherGlycerophospholipids(unittest.TestCase):
    """2,3-di-O-phytanyl-sn-glycero-1-phospho-(3'-sn-glycerol-1'-sulfate)"""
    POS = "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@](OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)([H])COP(OC[C@@](O)([H])COS(O)(=O)=O)(O)=O"
    NEG_STEREO = "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@](OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)([H])COP(OCC(O)([H])COS(O)(=O)=O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx44__3_OtherGlycerophospholipids(unittest.TestCase):
    """1-O-hexadecyl-sn-glycero-3-phosphoric acid methyl ester"""
    POS = "[H][C@](O)(COP(OC)(O)=O)COCCCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx45__3_OtherGlycerophospholipids(unittest.TestCase):
    """Cholesteryl-6-O-(1,2-ditetradecanoyl-sn-glycero-3-phospho)-α-D-glucopy"""
    POS = "O=C(CCCCCCCCCCCCC)OC[C@]([H])(COP(O)(OC[C@H]1O[C@H](O[C@@H]2CC3=CC[C@]4([H])[C@@](CC[C@@]5(C)[C@@]4([H])CC[C@@H]5[C@@H](CCCC(C)C)C)([H])[C@@]3(C)CC2)[C@H](O)[C@@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OC(CCCCCCCCCCCCC)=O)(OC[C@H]1O[C@H](O[C@@H]2CC3=CC[C@]4([H])[C@@](CC[C@@]5(C)[C@@]4([H])CC[C@@H]5[C@@H](CCCC(C)C)C)([H])[C@@]3(C)CC2)[C@H](OC(CCCCCCCCCCCCC)=O)[C@@H](O)[C@@H]1O)=O)O"
    VAR_SUGAR = "O=C(CCCCCCCCCCCCC)OC[C@]([H])(COP(O)(OC[C@H]1O[C@H](O[C@@H]2CC3=CC[C@]4([H])[C@@](CC[C@@]5(C)[C@@]4([H])CC[C@@H]5[C@@H](CCCC(C)C)C)([H])[C@@]3(C)CC2)[C@H](O)[C@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCC)OCC([H])(COP(O)(OC[C@H]1O[C@H](O[C@@H]2CC3=CC[C@]4([H])[C@@](CC[C@@]5(C)[C@@]4([H])CC[C@@H]5[C@@H](CCCC(C)C)C)([H])[C@@]3(C)CC2)[C@H](O)[C@@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx46__3_OtherGlycerophospholipids(unittest.TestCase):
    """1',3'-Bis-[1-(9Z-octadecenoyl)-2-hexadecanoyl-sn-glycero-3-phospho]-di"""
    POS = "O=P(O)(OC[C@@](OC(CCCCCCCCCCCCCCC)=O)([H])COC(CCCCCCC/C=C\\CCCCCCCC)=O)OCCNCCOP(O)(OC[C@@](OC(CCCCCCCCCCCCCCC)=O)([H])COC(CCCCCCC/C=C\\CCCCCCCC)=O)=O"
    NEG_STEREO = "O=P(O)(OCC(OC(CCCCCCCCCCCCCCC)=O)([H])COC(CCCCCCC/C=C\\CCCCCCCC)=O)OCCNCCOP(O)(OC[C@@](OC(CCCCCCCCCCCCCCC)=O)([H])COC(CCCCCCC/C=C\\CCCCCCCC)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx47__3_Glycerophosphocholines(unittest.TestCase):
    """1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphocholine"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCC[N+](C)(C)C)([O-])=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OCC[N+](C)(C)C)([O-])=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx48__3_Glycerophosphoethanolamines(unittest.TestCase):
    """1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphoethanolamine"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCCN)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OCCNC(CCCCCCCCCCCCCCC)=O)(OC(CCCCCCC/C=C\\CCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OCCN)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx49__3_Glycerophosphoserines(unittest.TestCase):
    """1-dodecanoyl-2-tridecanoyl-sn-glycero-3-phosphoserine"""
    POS = "O=C(CCCCCCCCCCC)OC[C@]([H])(COP(OC[C@](C(O)=O)([H])N)(O)=O)OC(CCCCCCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OC[C@](C(OC(CCCCCCCCCCC)=O)=O)([H])N)(OC(CCCCCCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCC)OCC([H])(COP(OC[C@](C(O)=O)([H])N)(O)=O)OC(CCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx50__3_Glycerophosphoglycerols(unittest.TestCase):
    """1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phospho-(1'-rac-glycer"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(OCC(O)CO)(O)=O"
    NEG_CHAIN = "OC[C@@]([H])(O)COP(OCC(O)COC(CCCCCCC/C=C\\CCCCCCCC)=O)(OC(CCCCCCCCCCCCCCC)=O)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(OCC(O)CO)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx51__3_Glycerophosphoglycerophosphates(unittest.TestCase):
    """1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phospho-(1'-sn-glycero"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OC[C@](O)([H])COP(O)(O)=O)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OC[C@](OC(CCCCCCCCCCCCCCC)=O)([H])COP(O)(OC(CCCCCCC/C=C\\CCCCCCCC)=O)=O)(O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OC[C@](O)([H])COP(O)(O)=O)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx52__3_Glycerophosphoinositols(unittest.TestCase):
    """1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phospho-(1'-myo-inosit"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O"
    NEG_CHAIN = "OC[C@@]([H])(O)COP(O[C@@H]1[C@H](OC(CCCCCCC/C=C\\CCCCCCCC)=O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(OC(CCCCCCCCCCCCCCC)=O)=O"
    VAR_SUGAR = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx52__3___mono_Glycerophosphoinositols(unittest.TestCase):
    """Glycerophosphoinositols"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCC/C=C\\CCCCCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx52__3___bi_Glycerophosphoinositols(unittest.TestCase):
    """Glycerophosphoinositols"""
    POS = "O=C(CCCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx52__3___tri_Glycerophosphoinositols(unittest.TestCase):
    """Glycerophosphoinositols"""
    POS = "O=C(CCCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)COP(O[C@@H]1[C@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx53__3_Glycerophosphates(unittest.TestCase):
    """1-dodecanoyl-2-tridecanoyl-sn-glycero-3-phosphate"""
    POS = "O=C(CCCCCCCCCCC)OC[C@]([H])(COP(O)(O)=O)OC(CCCCCCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OC(CCCCCCCCCCCC)=O)(OC(CCCCCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCC)OCC([H])(COP(O)(O)=O)OC(CCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx54__3_Glyceropyrophosphates(unittest.TestCase):
    """1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-pyrophosphate"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OP(O)(O)=O)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(OP(OC(CCCCCCC/C=C\\CCCCCCCC)=O)(O)=O)(OC(CCCCCCCCCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OP(O)(O)=O)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx55__3_Glycerophosphoglycerophosphoglycerols(unittest.TestCase):
    """1',3'-Bis[1,2-Di-(9Z-12Z-octadecadienoyl)-sn-glycero-3-phospho]-sn-gly"""
    POS = "O=P(O)(OC[C@@](OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)([H])COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)OCC([H])(O)COP(O)(OC[C@@](OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)([H])COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)=O"
    NEG_STEREO = "O=P(O)(OCC(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)([H])COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)OC[C@@]([H])(O)COP(O)(OC[C@@](OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)([H])COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx56__3_Cdpglycerols(unittest.TestCase):
    """1,2-Didodecanoyl-sn-glycero-3-cytidine-5'-diphosphate"""
    POS = "O=C(CCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCCCCCC)=O)COP(OP(OC[C@H]1O[C@@H](N2C=CC(N)=NC2=O)[C@H](O)[C@@H]1O)(O)=O)(O)=O"
    NEG_CHAIN = "OC[C@@]([H])(O)COP(OP(OC[C@H]1O[C@@H](N2C=CC(NC(CCCCCCCCCCC)=O)=NC2=O)[C@H](O)[C@@H]1OC(CCCCCCCCCCC)=O)(O)=O)(O)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCC)OCC([H])(OC(CCCCCCCCCCC)=O)COP(OP(OC[C@H]1O[C@@H](N2C=CC(N)=NC2=O)[C@H](O)[C@@H]1O)(O)=O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx57__3_Glycosylglycerophospholipids(unittest.TestCase):
    """1-octadecanoyl-2-(5Z,8Z,11Z,14Z-eicosatetraenoyl)-sn-glycero-3-phospho"""
    POS = "O=C(CCCCCCCCCCCCCCCCC)OC[C@]([H])(COP(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](OC(CCCCCCCCCCCCCCCCC)=O)[C@H]1O)(OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCCCC)OCC([H])(COP(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx58__3_Glycosylglycerophospholipids(unittest.TestCase):
    """1-octadecanoyl-2-eicosanoyl-sn-glycero-3-phospho-(1'-β-D-6'-O-acetyl-g"""
    POS = "O=C(CCCCCCCCCCCCCCCCC)OC[C@]([H])(COP(O[C@@H]1O[C@H](COC(C)=O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)OC(CCCCCCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCCCC)OCC([H])(COP(O[C@@H]1O[C@H](COC(C)=O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)OC(CCCCCCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx59__3_Glycosylglycerophospholipids(unittest.TestCase):
    """1-(11Z-octadecenoyl)-2-(hexadecenoyl)-sn-glycero-3-phospho-2'-α-D-6-gl"""
    POS = "[H][C@@](CO)(O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1N)COP(OC[C@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCC/C=C\\CCCCCC)=O)(O)=O"
    NEG_STEREO = "[H][C@@](CO)(O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1N)COP(OCC([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCC/C=C\\CCCCCC)=O)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx60__3_Glycosylglycerophospholipids(unittest.TestCase):
    """1-tetradecanoyl-2-hexadecanoyl-sn-glycero-3-phospho-(1'-glycerol-3'-)5"""
    POS = "O=C(CCCCCCCCCCCCC)OC[C@]([H])(COP(OCC(O)CO[C@@H]1O[C@H](C[As](C)(C)=O)[C@@H](O)[C@H]1O)(O)=O)OC(CCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "OC[C@]([H])(COP(OCC(OC(CCCCCCCCCCCCC)=O)CO[C@@H]1O[C@H](C[As](C)(C)=O)[C@@H](O)[C@H]1O)(OC(CCCCCCCCCCCCCCC)=O)=O)O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx61__3_Glycerophosphoinositolglycans(unittest.TestCase):
    """EtN-P-6Manα1-2Manα1-6Manα1-4GlcNα1-6-PI(14:0/14:0)"""
    POS = "O=C(CCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCCCCCCCC)=O)COP(O)(O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O[C@H]2O[C@H](CO)[C@@H](O[C@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@H]5O[C@H](COP(OCCN)(O)=O)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O)[C@@H]3O)[C@H](O)[C@H]2N)=O"
    NEG_CHAIN = "OC[C@@]([H])(O)COP(O)(O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O[C@H]2O[C@H](COC(CCCCCCCCCCCCC)=O)[C@@H](O[C@H]3O[C@H](CO[C@H]4O[C@H](COC(CCCCCCCCCCCCC)=O)[C@@H](O)[C@H](O)[C@@H]4O[C@H]5O[C@H](COP(OCCN)(O)=O)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O)[C@@H]3O)[C@H](O)[C@H]2N)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))


class TestEx62__3_Glycerophosphoinositolglycans(unittest.TestCase):
    """2'-O-(α-D-Manp)-(1-hexadecanoyl-2-tetradecanoyl-sn-glycero-3-phospho-1"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(O)(O[C@@H]1[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"
    VAR_SUGAR = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(O)(O[C@@H]1[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)OC(CCCCCCCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(O)(O[C@@H]1[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)=O)OC(CCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx63__3_DiglycerolTetraetherPhospholipidsCaldarc(unittest.TestCase):
    """sn-caldarchaeo-1-phosphoethanolamine"""
    POS = "NCCOP(O)(OC[C@]1([H])COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@](COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO1)([H])CO)=O"
    NEG_STEREO = "NCCOP(O)(OC[C@@]1([H])COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@](COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO1)([H])CO)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx64__3_GlycerolnonitolTetraetherPhospholipids(unittest.TestCase):
    """sn-caldito-1-phosphoethanolamine"""
    POS = "OC([C@](OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H]1C)([H])COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@](COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC1)([H])COP(OCCN)(O)=O)C(CO)(O)C(O)C(O)C(O)CO"
    NEG_STEREO = "OC([C@@](OCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@@H]1C)([H])COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCO[C@](COCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CCC[C@H](C)CC1)([H])COP(OCCN)(O)=O)C(CO)(O)C(O)C(O)C(O)CO"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx65__3_OxidizedGlycerophospholipids(unittest.TestCase):
    """1-O-(1Z-hexadecenyl)-2-(12S-hydroxy-5Z,8Z,10E,14Z-eicosatetraenoyl)-sn"""
    POS = "[H][C@](OC(CCC/C=C\\C/C=C\\C=C\\[C@@H](O)C/C=C\\CCCCC)=O)(COP(OCCN)(O)=O)CO/C=C\\CCCCCCCCCCCCCC"
    NEG_CHAIN = "[H][C@](O)(COP(OCCN)(OC(CCC/C=C\\C/C=C\\C=C\\[C@@H](O)C/C=C\\CCCCC)=O)=O)CO/C=C\\CCCCCCCCCCCCCC"
    NEG_STEREO = "[H]C(OC(CCC/C=C\\C/C=C\\C=C\\[C@@H](O)C/C=C\\CCCCC)=O)(COP(OCCN)(O)=O)CO/C=C\\CCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx66__3_GlycerophosphoethanolamineGlycans(unittest.TestCase):
    """N-(1-deoxyfructosyl)-1-hexadecanoyl-2-(4Z,7Z,10Z,13Z,16Z,19Z-docosahex"""
    POS = "O[C@]1(CNCCOP(OC[C@]([H])(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCCCCCCCCCC)=O)(O)=O)OC[C@@H](O)[C@@H](O)[C@@H]1O"
    VAR_SUGAR = "O[C@]1(CNCCOP(OC[C@]([H])(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCCCCCCCCCC)=O)(O)=O)OC[C@H](O)[C@H](O)[C@@H]1O"
    NEG_STEREO = "O[C@]1(CNCCOP(OCC([H])(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCCCCCCCCCC)=O)(O)=O)OC[C@@H](O)[C@@H](O)[C@@H]1O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx67__3_Dihydroxyacetonephospates(unittest.TestCase):
    """1-hexadecanoyl-glycerone 3-phosphate"""
    POS = "O=C(CCCCCCCCCCCCCCC)OCC(COP(O)(O)=O)=O"
    NEG_CHAIN = "OCC(COP(O)(OC(CCCCCCCCCCCCCCC)=O)=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))


class TestEx68__3_Glycerophosphoethanols(unittest.TestCase):
    """1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphoethanol"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@]([H])(COP(OCC)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(COP(OCC)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx69__3_Glycerophosphothreonines(unittest.TestCase):
    """1-octadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphothreonine"""
    POS = "O=C(CCCCCCCCCCCCCCCCC)OC[C@]([H])(COP(O[C@H](C)[C@](C(O)=O)([H])N)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"
    NEG_CHAIN = "OC[C@]([H])(COP(O[C@H](C)[C@](C(O)=O)([H])NC(CCCCCCC/C=C\\CCCCCCCC)=O)(OC(CCCCCCCCCCCCCCCCC)=O)=O)O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCCCC)OCC([H])(COP(O[C@H](C)[C@](C(O)=O)([H])N)(O)=O)OC(CCCCCCC/C=C\\CCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx71__3_OtherGlycerolipids(unittest.TestCase):
    """1-(9Z,1Z-octadecadienoyl)-2-(10Z,13Z,16Z,19Z-docosatetraenoyl)-3-O-[hy"""
    POS = "O=C(CCCCCCC/C=C\\C/C=C\\CCCCC)OCC([H])(COCC([H])(C([O-])=O)C[N+](C)(C)C)OC(CCCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx72__3_OtherGlycerolipids(unittest.TestCase):
    """1-O-(1'S,2'S,3'R,4'R,5'S-tetrahydroxycyclopentyl)-2-(9-methylpentadeca"""
    POS = "[H][C@@](COCCCCCCCCCC(C)CCCCCC)(OC(CCCCCCCC(C)CCCCCC)=O)COC1[C@@H](O)[C@H](O)[C@H]([C@@H]1O)O"
    NEG_CHAIN = "[H][C@@](COCCCCCCCCCC(C)CCCCCC)(O)CO[C@@H]1[C@@H](O)[C@H](O)[C@H]([C@@H]1O)OC(CCCCCCCC(C)CCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))


class TestEx73__3_OtherGlycerolipids(unittest.TestCase):
    """1-hexadecanoyl-2-(10-methyl-octadecanoyl)-3-O-(2S,5-diaminohexanoyl)-s"""
    POS = "O=C(CCCCCCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCCCC(C)CCCCCCCC)=O)COC([C@H](N)CCCCN)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCCCCCC)OCC([H])(OC(CCCCCCCCC(C)CCCCCCCC)=O)COC([C@H](N)CCCCN)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx74__3_Glycerolipids(unittest.TestCase):
    """1-dodecanoyl-sn-glycerol"""
    POS = "OC[C@]([H])(O)COC(CCCCCCCCCCC)=O"
    NEG_STEREO = "OCC([H])(O)COC(CCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx75__3_Glycerolipids(unittest.TestCase):
    """1-acyl-3-O-β-D-galactosyl-sn-glycerol"""
    POS = "O=C(C)OC[C@@]([H])(O)CO[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1"
    NEG_CHAIN = "OC[C@@]([H])(O)CO[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](COC(C)=O)O1"
    NEG_STEREO = "O=C(C)OCC([H])(O)CO[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx76__3_Glycerolipids(unittest.TestCase):
    """1-hexadecanyl-2-((2'-α-glucosyl)-β-glucosyl)-3-β-xylosyl-sn-glycerol"""
    POS = "[H][C@](COCCCCCCCCCCCCCCCC)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)CO[C@H]3O[C@@H](CO)[C@@H](O)[C@@H]3O"
    NEG_CHAIN = "[H][C@](COCCCCCCCCCCCCCCCC)(O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O[C@H]2O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]2O)CO[C@H]3O[C@@H](CO)[C@@H](O)[C@@H]3O"
    NEG_STEREO = "[H]C(COCCCCCCCCCCCCCCCC)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)CO[C@H]3O[C@@H](CO)[C@@H](O)[C@@H]3O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx77__3_Glycerolipids(unittest.TestCase):
    """1-(9Z-hexadecenoyl)-3-(6'-sulfo-α-D-quinovosyl)-sn-glycerol"""
    POS = "[H][C@](O)(CO[C@H]1O[C@H](CS(=O)(O)=O)[C@@H](O)[C@H](O)[C@H]1O)COC(CCCCCCC/C=C\\CCCCCC)=O"
    NEG_CHAIN = "[H]C(O)(CO[C@H]1O[C@H](CS(=O)(O)=O)[C@@H](O)[C@H](O)[C@H]1O)COC(CCCCCCC/C=C\\CCCCCC)=O"
    VAR_SUGAR = "[H][C@](O)(CO[C@H]1O[C@@H](CS(=O)(O)=O)[C@@H](O)[C@@H](O)[C@H]1O)COC(CCCCCCC/C=C\\CCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))


class TestEx78__3_Betaine(unittest.TestCase):
    """1-(9Z-octadecenoyl)-2-(10Z,13Z,16Z,19Z-docosatetraenoyl)-sn-glycero-3-"""
    POS = "O=C(CCCCCCC/C=C\\CCCCCCCC)OC[C@]([H])(COCC(C([O-])=O)C[N+](C)(C)C)OC(CCCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O"
    NEG_STEREO = "O=C(CCCCCCC/C=C\\CCCCCCCC)OCC([H])(COCC(C([O-])=O)C[N+](C)(C)C)OC(CCCCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx79__3_Phosphosphingolipids(unittest.TestCase):
    """N-(eicosanoyl)-hexadecasphing-4-enine-1-phosphocholine"""
    POS = "[H][C@](NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])/C=C/CCCCCCCCCCC"
    NEG_STEREO = "[H]C(NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])/C=C/CCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx80__3_Phosphosphingolipids(unittest.TestCase):
    """N-(2R-hydroxy--15-methyl-3E-hexadecenoyl)-9-methyl-eicosasphinga-4E,8E"""
    POS = "[H][C@](NC([C@H](O)/C=C/CCCCCCCCCCC(C)C)=O)(COP(OCCN)(O)=O)[C@@](O)([H])/C=C/CC/C=C(C)/CCCCCCCCCCC"
    NEG_CHAIN = "[H][C@](N)(COP(OCCN)(OC([C@H](O)/C=C/CCCCCCCCCCC(C)C)=O)=O)[C@@](O)([H])/C=C/CC/C=C(C)/CCCCCCCCCCC"
    NEG_STEREO = "[H]C(NC([C@H](O)/C=C/CCCCCCCCCCC(C)C)=O)(COP(OCCN)(O)=O)[C@@](O)([H])/C=C/CC/C=C(C)/CCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx81__3_Phosphosphingolipids(unittest.TestCase):
    """N-(docosanoyl)-sphing-4-enine-1-phospho-(1'-myo-inositol)"""
    POS = "[H][C@](NC(CCCCCCCCCCCCCCCCCCCCC)=O)(COP(O)(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"
    VAR_SUGAR = "[H][C@](NC(CCCCCCCCCCCCCCCCCCCCC)=O)(COP(O)(O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"
    NEG_STEREO = "[H]C(NC(CCCCCCCCCCCCCCCCCCCCC)=O)(COP(O)(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx82__3_NeutralGlycosphingolipids(unittest.TestCase):
    """N-(34S-methyl-5Z,9Z,12Z,15Z,18Z,21Z-hexatriacontahexaenoyl)-1-sulfo-β-"""
    POS = "O=C(N[C@@]([H])(COS(O[C@@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O)(=O)=O)[C@@](O)([H])/C=C/CCCCCCCCCC[C@H](C)CC)CCC/C=C\\CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCCCCCCCC[C@@H](C)CC"
    NEG_CHAIN = "[H][C@](N)(COS(O[C@H]1O[C@H](C)[C@H](OC(CCC/C=C\\CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCCCCCCCC[C@@H](C)CC)=O)[C@H](O)[C@H]1O)(=O)=O)[C@@](O)([H])/C=C/CCCCCCCCCC[C@H](C)CC"
    VAR_SUGAR = "O=C(N[C@@]([H])(COS(O[C@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O)(=O)=O)[C@@](O)([H])/C=C/CCCCCCCCCC[C@H](C)CC)CCC/C=C\\CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCCCCCCCC[C@@H](C)CC"
    NEG_STEREO = "O=C(NC([H])(COS(O[C@@H]1O[C@H](C)[C@H](O)[C@H](O)[C@H]1O)(=O)=O)[C@@](O)([H])/C=C/CCCCCCCCCC[C@H](C)CC)CCC/C=C\\CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCCCCCCCC[C@@H](C)CC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx83__3_NeutralGlycosphingolipids(unittest.TestCase):
    """1-O-melibiosoyl-(N-(2R-hydroxy-heneicosanoyl)-4R-hydroxy-17-methyl-sph"""
    POS = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCC)=O)(CO[C@@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@H](O)/C=C/CCCCCCCCCC(C)C"
    NEG_CHAIN = "[H][C@](N)(CO[C@@H]1O[C@H](CO[C@H]2O[C@H](COC([C@H](O)CCCCCCCCCCCCCCCCCCC)=O)[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@H](O)/C=C/CCCCCCCCCCC.C"
    VAR_SUGAR = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@H](O)/C=C/CCCCCCCCCC(C)C"
    NEG_STEREO = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCC)=O)(CO[C@@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)C(O)([H])[C@H](O)/C=C/CCCCCCCCCC(C)C"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx84__3_AcidicGlycosphingolipids(unittest.TestCase):
    """1-O-(6'-phosphoethanolaminy-β-D-glucopyranosyl)-(N-(2R-hydroxy-tetraco"""
    POS = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O)(CO[C@@H]1O[C@H](COP(O)(OCCN)=O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@@H](CCCCCCCCCCC(C)C)O"
    NEG_CHAIN = "[H][C@](N)(CO[C@@H]1O[C@H](COP(O)(OCCN)=O)[C@@H](O)[C@H](O)[C@H]1OC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O)[C@@](O)([H])[C@@H](CCCCCCCCCCC(C)C)O"
    VAR_SUGAR = "[H][C@](NC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](COP(O)(OCCN)=O)[C@@H](O)[C@H](O)C1)[C@@](O)([H])[C@@H](CCCCCCCCCCC(C)C)O"
    NEG_STEREO = "[H]C(NC([C@H](O)CCCCCCCCCCCCCCCCCCCCCC)=O)(CO[C@@H]1O[C@H](COP(O)(OCCN)=O)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])[C@@H](CCCCCCCCCCC(C)C)O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_chain(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_CHAIN)
        ref = LipidHeadValidator.get_inchi(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_inchi=ref))

    def test_sugar_variant_still_matches(self):
        m = LipidAnalysis.parse_smiles(self.VAR_SUGAR)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


class TestEx85__3_OtherPolyketides(unittest.TestCase):
    """1S-((4S-acetoxy-5R-methyl-3-methylene-6-phenylhexyl)-6-(E)-4S,6S-dimet"""
    POS = "O[C@@H]1[C@](O[C@H](C(O)=O)[C@]2(C(O)=O)O)(CCC([C@H]([C@H](C)CC3=CC=CC=C3)OC(C)=O)=C)O[C@]2(C(O)=O)[C@@H]1OC(/C=C/[C@@H](C)C[C@@H](C)CC)=O"
    NEG_STEREO = "O[C@@H]1[C@](O[C@H](C(O)=O)[C@@]2(C(O)=O)O)(CCC([C@H]([C@H](C)CC3=CC=CC=C3)OC(C)=O)=C)O[C@]2(C(O)=O)[C@@H]1OC(/C=C/[C@@H](C)C[C@@H](C)CC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = LipidAnalysis.parse_smiles(self.POS)
        self.assertTrue(self.v.matches_any_valid_head(m))

    def test_negative_stereo(self):
        pos = LipidAnalysis.parse_smiles(self.POS)
        neg = LipidAnalysis.parse_smiles(self.NEG_STEREO)
        ref = InChIParser.get_stereo_layer(pos)
        self.assertFalse(self.v.matches_any_valid_head(neg, reference_stereo_inchi=ref))


if __name__ == "__main__":
    unittest.main(verbosity=2)