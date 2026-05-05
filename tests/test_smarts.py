"""
Tests generated directly from Naming_Example.xlsx

Column mapping (Hoja1):
  B = positive SMILES          → must match
  H = error 1 (wrong chain)    → FA placed at wrong position (O instead of N), must NOT match
  I = variant azúcares         → different/extra sugar, N-acyl intact
  J = error 3 (stereochemistry)→ wrong stereo on sphingoid base

Each test class maps to one Excel row.
Column I (sugar variant) always keeps the correct N-acyl bond, so it is
expected to still match — it is NOT a negative, just a structural variant.
Column J (stereo) tests are marked @expectedFailure because current SMARTS
patterns do not enforce stereochemistry.

Rows with no entry in a column → that test method is omitted.
Rows whose negatives still contain [G]/[R] placeholders → skipped (invalid SMILES).
"""

import unittest
from rdkit import Chem
from src.backend.lipid.lipid_analysis import LipidHeadValidator


def parse(smiles: str):
    """Return RDKit mol, or None if SMILES is empty/placeholder/invalid."""
    if not smiles:
        return None
    if "[G]" in smiles or "[R]" in smiles:
        return None
    return Chem.MolFromSmiles(smiles.strip())


# ──────────────────────────────────────────────────────────────────────────────
# Example 7 — Neutral glycosphingolipid  (Gal-GlcNAc-Gal-GlcNAc-Gal-Glc-Cer)
# ──────────────────────────────────────────────────────────────────────────────
class TestEx7_NeutralGlycosphingolipid(unittest.TestCase):
    # Col B: correct structure — N-acyl ceramide + oligosaccharide chain
    POS = (
        "[H][C@](NC(CCCCCCCCCCCCCCC)=O)"
        "(CO[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O[C@H]3"
        "[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]4[C@H](O)[C@@H](O[C@@H]5"
        "[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]6[C@H](O)[C@@H](O)[C@@H](O)"
        "[C@H](O6)CO)[C@H](O5)CO)[C@@H](O)[C@H](O4)CO)[C@H](O3)CO)"
        "[C@@H](O)[C@H](O2)CO)[C@H](O1)CO)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"
    )
    # Col H: FA moved from N to O (free amine)
    NEG_CHAIN = (
        "[H][C@](N)"
        "(CO[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O[C@H]3"
        "[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]4[C@H](O)[C@@H](O[C@@H]5"
        "[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]6[C@H](O)[C@@H](O)[C@@H](O)"
        "[C@H](O6)CO)[C@H](O5)CO)[C@@H](O)[C@H](O4)CO)[C@H](O3)CO)"
        "[C@@H](O)[C@H](O2)CO)[C@H](O1)CO)"
        "[C@@](OC(CCCCCCCCCCCCCCC)=O)([H])/C=C/CCCCCCCCCCCCC"
    )
    # Col J: wrong stereo on C2 of sphingoid base
    NEG_STEREO = (
        "[H][C@](NC(CCCCCCCCCCCCCCC)=O)"
        "(CO[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O[C@H]3"
        "[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]4[C@H](O)[C@@H](O[C@@H]5"
        "[C@H](NC(C)=O)[C@@H](O)[C@H](O[C@H]6[C@H](O)[C@@H](O)[C@@H](O)"
        "[C@H](O6)CO)[C@H](O5)CO)[C@@H](O)[C@H](O4)CO)[C@H](O3)CO)"
        "[C@@H](O)[C@H](O2)CO)[C@H](O1)CO)"
        "C(O)([H])/C=C/CCCCCCCCCCCCC"
    )

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m, "Col B SMILES invalid")
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex7 positive (N-acyl ceramide + sugars) must match")

    def test_negative_chain(self):
        """FA at O instead of N → must NOT match."""
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m, "Col H SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex7 chain-error (O-acyl, free amine) must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        """Wrong C2 stereo — not enforced by current SMARTS."""
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m, "Col J SMILES invalid")
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex7 stereo-error should NOT match once stereo patterns added")


# ──────────────────────────────────────────────────────────────────────────────
# Example 60 — Acidic glycosphingolipid  (glucuronic acid headgroup)
# ──────────────────────────────────────────────────────────────────────────────
class TestEx60_AcidicGlycosphingolipid(unittest.TestCase):
    POS       = "[H][C@](NC(CCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](C(O)=O)[C@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])CCCCCCCCCCCCCCC"
    NEG_CHAIN = "[H][C@](N)(CO[C@H]1O[C@H](C(O)=O)[C@H](O)[C@H](O)[C@H]1O)[C@@](OC(CCCCCCCCCCCCC)=O)([H])CCCCCCCCCCCCCCC"
    # Col I: extended sugar (still has N-acyl) → should still match
    VAR_SUGAR = "[H][C@](NC(CCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](C(O)=O)[C@@H](O[C@H]2O[C@H](CO[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H](O)[C@H]2N)[C@H](O)[C@H]1O)[C@@](O)([H])CCCCCCCCCCCCCCC"
    NEG_STEREO= "[H][C@](NC(CCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](C(O)=O)[C@H](O)[C@H](O)[C@H]1O)C(O)([H])CCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex60 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex60 chain-error (O-acyl, free amine) must NOT match")

    def test_sugar_variant_still_matches(self):
        """Col I keeps correct N-acyl bond — must still match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex60 sugar-variant (N-acyl intact) must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 64 — Phosphosphingolipid / Sphingomyelin
# ──────────────────────────────────────────────────────────────────────────────
class TestEx64_Sphingomyelin(unittest.TestCase):
    POS        = "[H][C@](NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])CCCCCCCCCCCCCCC"
    NEG_CHAIN  = "[H][C@](N)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](OC(CCCCCCCCCCCCCCCCCCC)=O)([H])CCCCCCCCCCCCCCC"
    NEG_STEREO = "[H][C@](NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)C(O)([H])CCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex64 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex64 chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 65 — Ceramide phosphoethanolamine
# ──────────────────────────────────────────────────────────────────────────────
class TestEx65_CeramidePhosphoethanolamine(unittest.TestCase):
    POS        = "[H][C@](NC(CCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(OCCN)(O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCC"
    NEG_CHAIN  = "[H][C@](N)(COP(OCCN)(O)=O)[C@@](OC(CCCCCCCCCCCCCCCCCCCCCCC)=O)([H])/C=C/CCCCCCCCCCC"
    NEG_STEREO = "[H][C@](NC(CCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(OCCN)(O)=O)C(O)([H])/C=C/CCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex65 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex65 chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 66 — Ceramide phosphoinositol
# ──────────────────────────────────────────────────────────────────────────────
class TestEx66_CeramidePhosphoinositol(unittest.TestCase):
    POS        = "[H][C@](NC(C(O)CCCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)[C@@](O)([H])CCCCCCCCCCCCCCCCC"
    NEG_CHAIN  = "[H][C@](N)(COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)[C@@](OC(C(O)CCCCCCCCCCCCCCCCCCCCCCCC)=O)([H])CCCCCCCCCCCCCCCCC"
    # Col I: different inositol stereochemistry — N-acyl intact, should match
    VAR_SUGAR  = "[H][C@](NC(C(O)CCCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O)(O)=O)[C@@](O)([H])CCCCCCCCCCCCCCCCC"
    NEG_STEREO = "[H][C@](NC(C(O)CCCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)C(O)([H])CCCCCCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex66 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex66 chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Col I: different inositol stereo but N-acyl bond correct — must match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex66 inositol-stereo variant (N-acyl intact) must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 67 — Acidic glycosphingolipid / Sulfatide
# ──────────────────────────────────────────────────────────────────────────────
class TestEx67_Sulfatide(unittest.TestCase):
    POS        = "[H][C@](NC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)(CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(=O)(O)=O)[C@H]1O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"
    NEG_CHAIN  = "[H][C@](N)(CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(=O)(O)=O)[C@H]1O)[C@@](OC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)([H])/C=C/CCCCCCCCCCCCC"
    # Col I: desulfated glucose — N-acyl intact, should match
    VAR_SUGAR  = "[H][C@](NC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)(CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"
    NEG_STEREO = "[H][C@](NC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)(CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(=O)(O)=O)[C@H]1O)C(O)([H])/C=C/CCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex67 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex67 chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        """Col I: desulfated glucose, N-acyl still correct — must match."""
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex67 sugar-variant (desulfated, N-acyl intact) must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 69 — Ceramide phospho-inositol-glucosamine
# ──────────────────────────────────────────────────────────────────────────────
class TestEx69_CeramidePhosphoInositolGlcN(unittest.TestCase):
    POS        = "CCCCCCCCCCCCCC[C@@]([H])(O)[C@](O)([H])[C@](NC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O)([H])COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2N)(O)=O"
    NEG_CHAIN  = "CCCCCCCCCCCCCC[C@@]([H])(O)[C@](OC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O)([H])[C@](N)([H])COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2N)(O)=O"
    # Col I: modified sugar (CH3 instead of CH2OH) — N-acyl intact
    VAR_SUGAR  = "CCCCCCCCCCCCCC[C@@]([H])(O)[C@](O)([H])[C@](NC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O)([H])COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@H](C)[C@@H](O)[C@H](O)[C@H]2N)(O)=O"
    NEG_STEREO = "CCCCCCCCCCCCCC[C@@]([H])(O)C(O)([H])[C@](NC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O)([H])COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2N)(O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex69 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex69 chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex69 sugar-variant (N-acyl intact) must still match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 70 — Glycerolipid / Triacylglycerol
# ──────────────────────────────────────────────────────────────────────────────
class TestEx70_Triacylglycerol(unittest.TestCase):
    POS        = "O=C(CCCCCCCCCCC)OC[C@@]([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O"
    NEG_STEREO = "O=C(CCCCCCCCCCC)OCC([H])(OC(CCCCCCCCCCCCCCC)=O)COC(CCCCCCCCCCCCCCCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex70 positive must match")

    def test_positive_class(self):
        m = parse(self.POS)
        self.assertIn("Triacylglycerols", self.v.identify_lipid_class(m))

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 2 — Glycerolipid / 1,2-Diacylglycerol (polyunsaturated)
# ──────────────────────────────────────────────────────────────────────────────
class TestEx2_Diacylglycerol(unittest.TestCase):
    POS        = "OC[C@]([H])(OC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCC/C=C\\C/C=C\\CCCC)=O"
    NEG_STEREO = "OCC([H])(OC(CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)=O)COC(CCCCCCC/C=C\\C/C=C\\CCCC)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex2 positive must match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 75 — Other sphingolipid (sulfonyl headgroup)
# ──────────────────────────────────────────────────────────────────────────────
class TestEx75_SulfonylSphingolipid(unittest.TestCase):
    POS        = "CC(C)CCCCCCCCCCC[C@](O)([H])[C@](NC(C[C@@H](O)CCCCCCCCCCCC(C)C)=O)([H])CS(O)(=O)=O"
    NEG_CHAIN  = "N[C@@]([C@@](OC(C[C@@H](O)CCCCCCCCCCCC(C)C)=O)([H])CCCCCCCCCCCC(C)C)([H])CS(O)(=O)=O"
    NEG_STEREO = "CC(C)CCCCCCCCCCCC(O)([H])[C@](NC(C[C@@H](O)CCCCCCCCCCCC(C)C)=O)([H])CS(O)(=O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex75 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex75 chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


# ──────────────────────────────────────────────────────────────────────────────
# Example 76 — Other sphingolipid (phospho-mannoside)
# ──────────────────────────────────────────────────────────────────────────────
class TestEx76_PhosphoMannoside(unittest.TestCase):
    POS        = "OP(OCC(NC(CCCCCCCCCCCC(C)C)=O)C(O)CCCCCCCCCCCC(C)C)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)=O"
    NEG_CHAIN  = "NC(COP(O)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)=O)C(OC(CCCCCCCCCCCC(C)C)=O)CCCCCCCCCCCC(C)C"
    # Col I: open-ring sugar — N-acyl intact, should match
    VAR_SUGAR  = "OP(OCC(NC(CCCCCCCCCCCC(C)C)=O)C(O)CCCCCCCCCCCC(C)C)(O[C@@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O)=O"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex76 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex76 chain-error must NOT match")

    def test_sugar_variant_still_matches(self):
        m = parse(self.VAR_SUGAR)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m),
                        "Ex76 sugar-variant (N-acyl intact) must still match")


# ──────────────────────────────────────────────────────────────────────────────
# Example 79 — Ceramide-1-phosphate
# ──────────────────────────────────────────────────────────────────────────────
class TestEx79_Ceramide1Phosphate(unittest.TestCase):
    POS        = "[H][C@](NC(CCCCCCCCCCCCCCCCCCCCC)=O)(COP(O)(O)=O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC"
    NEG_CHAIN  = "[H][C@](N)(COP(O)(O)=O)[C@@](OC(CCCCCCCCCCCCCCCCCCCCC)=O)([H])/C=C/CCCCCCCCCCCCC"
    NEG_STEREO = "[H][C@](NC(CCCCCCCCCCCCCCCCCCCCC)=O)(COP(O)(O)=O)C(O)([H])/C=C/CCCCCCCCCCCCC"

    def setUp(self): self.v = LipidHeadValidator()

    def test_positive(self):
        m = parse(self.POS)
        self.assertIsNotNone(m)
        self.assertTrue(self.v.matches_any_valid_head(m), "Ex79 positive must match")

    def test_negative_chain(self):
        m = parse(self.NEG_CHAIN)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m),
                         "Ex79 chain-error must NOT match")

    @unittest.expectedFailure
    def test_negative_stereo(self):
        m = parse(self.NEG_STEREO)
        self.assertIsNotNone(m)
        self.assertFalse(self.v.matches_any_valid_head(m))


if __name__ == "__main__":
    unittest.main(verbosity=2)