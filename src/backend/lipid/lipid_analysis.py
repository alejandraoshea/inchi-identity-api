import sys
import requests
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from rdkit.Chem import inchi as rdInchi
from collections import Counter
from backend.lipid.lipid_tail_extraction import TailExtractor
from backend.lipid.lipid_pattern_generator import build_combined_patterns
from backend.inchi.inchi_parser import InChIParser
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class HeadgroupPattern:
    """
    Represents a lipid headgroup with positional specificity.
    Attributes:
        name:         Human-readable name
        smarts:       SMARTS pattern defining the headgroup
        lipid_class:  Classification (e.g. "Glycerophospholipids")
        fa_positions: sn positions where FAs attach (e.g. ["sn-1", "sn-2"])
        description:  Additional details
    """
    name: str
    smarts: str
    lipid_class: str
    fa_positions: List[str]
    description: str = ""


class LipidHeadValidator:
    MANUAL_PATTERNS = {
        "phosphatidylcholine": HeadgroupPattern(
            name="Phosphatidylcholine (PC)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CH2X4][NX4+]([CH3X4])([CH3X4])[CH3X4])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylcholines",
            fa_positions=["sn-1", "sn-2"],
        ),
        "phosphatidylethanolamine": HeadgroupPattern(
            name="Phosphatidylethanolamine (PE)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CH2X4][NX3H2,NX4H3+])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylethanolamines",
            fa_positions=["sn-1", "sn-2"],
        ),
        "phosphatidylserine": HeadgroupPattern(
            name="Phosphatidylserine (PS)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CHX4]([NX3H2,NX4H3+])[CX3](=[OX1])[OX2H0,OX1-,OX1H1])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylserines",
            fa_positions=["sn-1", "sn-2"],
        ),
        "phosphatidylglycerol": HeadgroupPattern(
            name="Phosphatidylglycerol (PG)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CHX4]([OX2H1])[CH2X4][OX2H1])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylglycerols",
            fa_positions=["sn-1", "sn-2"],
        ),
        "cardiolipin": HeadgroupPattern(
            name="Cardiolipin (CL)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CHX4][CH2X4])[OX2])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Cardiolipins",
            fa_positions=["sn-1", "sn-2", "sn-1'", "sn-2'"],
        ),

        "gpl_diacyl": HeadgroupPattern(
            name="Diacyl-glycerophospholipid (generic)",
            smarts="[CH2X4]([OX2][CX3](=[OX1])[#6])[CHX4]([OX2][CX3](=[OX1])[#6])[CH2X4][OX2][PX4]",
            lipid_class="Glycerophospholipids",
            fa_positions=[],
            description="Covers all diacyl-GPL variants regardless of headgroup identity",
        ),
        "gpl_plasmalogen_full": HeadgroupPattern(
            name="Plasmalogen (vinyl ether at sn-1)",
            smarts="[CH2X4]([OX2]/[CH]=[CH]/[#6])[CHX4]([OX2][CX3](=[OX1])[#6])[CH2X4][OX2][PX4]",
            lipid_class="Glycerophospholipids",
            fa_positions=[],
        ),
        "gpl_ether_chain": HeadgroupPattern(
            name="Ether-linked glycerophospholipid",
            smarts="[CX4][OX2][CX4]~[CX4]~[CX4]~[CX4]~[CX4].[CH2X4][OX2][PX4]",
            lipid_class="Glycerophospholipids",
            fa_positions=[],
        ),
        "gpl_dhap": HeadgroupPattern(
            name="DHAP-type glycerophospholipid",
            smarts="[CH2X4][OX2][CX3](=[OX1])[#6].[CH2X4][OX2][PX4]",
            lipid_class="Glycerophospholipids",
            fa_positions=[],
        ),

        "ceramide_nacyl_relaxed": HeadgroupPattern(
            name="Ceramide N-acyl (sphingoid context required)",
            smarts="[NX3H1]([CX3](=[OX1])[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6])[CX4H1][CX4][OX2H1]",
            lipid_class="Ceramides",
            fa_positions=["N-acyl"],
            description="N-amide must be on CH adjacent to CH2-OH (sphingoid base context).",
        ),
        "ceramide": HeadgroupPattern(
            name="Ceramide",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]",
            lipid_class="Ceramides",
            fa_positions=["N-acyl"],
        ),
        "sphingomyelin": HeadgroupPattern(
            name="Sphingomyelin (SM)",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]",
            lipid_class="Sphingomyelins",
            fa_positions=["N-acyl"],
        ),
        "ceramide_phosphoinositol": HeadgroupPattern(
            name="Ceramide phosphoinositol",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4].[PX4](=[OX1])[OX2][CH]1[CH]([OX2H1])[CH]([OX2H1])[CH]([OX2H1])[CH]([OX2H1])[CH]1[OX2H1]",
            lipid_class="Phosphosphingolipids",
            fa_positions=["N-acyl"],
        ),
        "sphingoid_base_free": HeadgroupPattern(
            name="Free sphingoid base",
            smarts="[NX3H2][CX4][CH2X4][OX2H1]",
            lipid_class="Sphingolipids",
            fa_positions=[],
        ),
        "ganglioside_core": HeadgroupPattern(
            name="Ganglioside core",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4].[CH]1O[CH](CO)[CH](O)[CH](O)[CH]1O",
            lipid_class="Gangliosides",
            fa_positions=["N-acyl"],
        ),
        "glycosphingolipid_galactose": HeadgroupPattern(
            name="Glycosphingolipid with Galactose",
            smarts="[H][C@](NC(CCCCCCCCCCCCCCC)=O)(CO[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC",
            lipid_class="Neutral glycosphingolipids",
            fa_positions=["N-acyl"],
        ),
        "glycosphingolipid_glucose": HeadgroupPattern(
            name="Glycosphingolipid with Glucose",
            smarts="[H][C@](NC(CCCCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC",
            lipid_class="Neutral glycosphingolipids",
            fa_positions=["N-acyl"],
        ),
        "acidic_glycosphingolipid_glucuronic": HeadgroupPattern(
            name="Acidic Glycosphingolipid",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]",
            lipid_class="Acidic glycosphingolipids",
            fa_positions=["N-acyl"],
        ),
        "sphingomyelin_detailed": HeadgroupPattern(
            name="Sphingomyelin (detailed)",
            smarts="[H][C@](NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])CCCCCCCCCCCCCCC",
            lipid_class="Phosphosphingolipids",
            fa_positions=["N-acyl"],
        ),
        "ceramide_phosphoethanolamine": HeadgroupPattern(
            name="Ceramide Phosphoethanolamine",
            smarts="[H][C@](NC(CCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(OCCN)(O)=O)[C@@](O)([H])CCCCCCCCCCCCCCC",
            lipid_class="Phosphosphingolipids",
            fa_positions=["N-acyl"],
        ),
        "ceramide_phosphoinositol_detailed": HeadgroupPattern(
            name="Ceramide Phosphoinositol (detailed)",
            smarts="[H][C@](NC(C(O)CCCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)[C@@](O)([H])CCCCCCCCCCCCCCCCC",
            lipid_class="Phosphosphingolipids",
            fa_positions=["N-acyl"],
        ),

        "saccharolipid_n_acyl": HeadgroupPattern(
            name="Saccharolipid (N-acyl aminosugar)",
            smarts="[C;R][NX3][CX3](=[OX1])[#6]~[#6]~[#6]",
            lipid_class="Saccharolipids",
            fa_positions=[],
        ),
        "saccharolipid_o_acyl": HeadgroupPattern(
            name="Saccharolipid (O-acyl sugar)",
            smarts="[C;R][OX2][CX3](=[OX1])[#6]~[#6]~[#6]~[#6]",
            lipid_class="Saccharolipids",
            fa_positions=[],
            description="Ester on sugar ring carbon. Post-match: ring must contain O (pyranose/furanose); "
                        "molecule must not have free sphingoid base.",
        ),

        "sterol_ester_generic": HeadgroupPattern(
            name="Sterol ester (ring O-acyl)",
            smarts="[C;R][OX2][CX3](=[OX1])[CX4,CX3]~[CX4]~[CX4]",
            lipid_class="Sterols",
            fa_positions=[],
            description="Post-match: requires ≥3 rings (steroid nucleus) and no phosphate.",
        ),
        "sterol_core": HeadgroupPattern(
            name="Sterol core",
            smarts="[C@]12CC[C@H]3[C@@H](CC[C@]4(C)[C@@H]3CC=C4)C1CCC2",
            lipid_class="Sterols",
            fa_positions=[],
        ),

        "glycerolipid_o_glycosyl": HeadgroupPattern(
            name="O-glycosylglycerolipid",
            smarts="[CX4][OX2][CH]1O[CH][CH][CH][CH]1",
            lipid_class="Glycerolipids",
            fa_positions=[],
        ),
        "glycerolipid_monoacyl": HeadgroupPattern(
            name="Monoacyl-glycerolipid (non-phosphate sn-3)",
            smarts="[CH2X4]([OX2][CX3](=[OX1])[#6])[CHX4]([OX2,OH1])[CH2X4][OX2;!$([OX2][PX4])]",
            lipid_class="Glycerolipids",
            fa_positions=[],
        ),
        "glycerolipid_acyl_ch": HeadgroupPattern(
            name="Acyl-glycerolipid (ester on central CH)",
            smarts="[CHX4]([OX2][CX3](=[OX1])[#6])([CH2X4])[CH2X4][OX2;!$([OX2][PX4])]",
            lipid_class="Glycerolipids",
            fa_positions=[],
        ),
        "glycosylglycerol_sn1": HeadgroupPattern(
            name="Glycosylglycerol (FA at sn-1)",
            smarts="O[C@H]1[C@H](OC[C@H](O)C[OX2][CX3](=O)[CX4,CX3]~[CX4,CX3]~[CX4,CX3])O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O",
            lipid_class="Glycoglycerolipids",
            fa_positions=["sn-1"],
        ),
        "diacylglycerol_1_2": HeadgroupPattern(
            name="1,2-Diacylglycerol",
            smarts="[OX2H1][CH2X4][CHX4]([OX2][CX3](=[OX1])[#6])[CH2X4][OX2][CX3](=[OX1])[#6]",
            lipid_class="Diacylglycerols",
            fa_positions=["sn-1", "sn-2"],
        ),
        "diacylglycerol_1_3": HeadgroupPattern(
            name="1,3-Diacylglycerol",
            smarts="[OX2][CX3](=[OX1])[#6].[CH2X4][CHX4]([OX2H1])[CH2X4][OX2][CX3](=[OX1])[#6]",
            lipid_class="Diacylglycerols",
            fa_positions=["sn-1", "sn-3"],
        ),
        "monoacylglycerol_sn1": HeadgroupPattern(
            name="1-Monoacylglycerol",
            smarts="[OX2H1][CH2X4][CHX4]([OX2H1])[CH2X4][OX2][CX3](=[OX1])[#6]",
            lipid_class="Monoacylglycerols",
            fa_positions=["sn-1"],
        ),
        "monoacylglycerol_sn2": HeadgroupPattern(
            name="2-Monoacylglycerol",
            smarts="[OX2H1][CH2X4][CHX4]([OX2][CX3](=[OX1])[#6])[CH2X4][OX2H1]",
            lipid_class="Monoacylglycerols",
            fa_positions=["sn-2"],
        ),
        "triacylglycerol": HeadgroupPattern(
            name="Triacylglycerol",
            smarts="[CH2X4]([OX2][CX3](=[OX1])[#6])[CHX4]([OX2][CX3](=[OX1])[#6])[CH2X4][OX2][CX3](=[OX1])[#6]",
            lipid_class="Triacylglycerols",
            fa_positions=["sn-1", "sn-2", "sn-3"],
        ),

        "fatty_acid": HeadgroupPattern(
            name="Fatty acid",
            smarts="[CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CX3](=O)[OX2H1]",
            lipid_class="Fatty acids",
            fa_positions=[],
        ),
        "plasmalogen": HeadgroupPattern(
            name="Plasmalogen",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][CH]=[CH])[OX2][PX4]",
            lipid_class="Ether phospholipids",
            fa_positions=["sn-1", "sn-2"],
        ),

        "prenol_lactone": HeadgroupPattern(
            name="Prenol lipid lactone",
            smarts="[O;R][CX3](=[OX1])[#6]",
            lipid_class="Prenol lipids",
            fa_positions=[],
        ),
        "prenol_allylic_ester": HeadgroupPattern(
            name="Prenol allylic ester",
            smarts="[OX2][CX3](=[OX1])[CX4,CX3][CX3]=[CX3]",
            lipid_class="Prenol lipids",
            fa_positions=[],
        ),
        "prenol_polyene_ester": HeadgroupPattern(
            name="Prenol polyene ester (retinoid)",
            smarts="[CX3]=[CX3][CX4][OX2][CX3](=[OX1])",
            lipid_class="Prenol lipids",
            fa_positions=[],
        ),
        "prenol_amide": HeadgroupPattern(
            name="Prenol lipid amide (quinone derivative)",
            smarts="[NX3][CX3](=[OX1])[CX4]~[CX4]~[CX4].[c,C][OX2H1]",
            lipid_class="Prenol lipids",
            fa_positions=[],
        ),
    }

    HEADGROUP_PATTERNS = build_combined_patterns(
        manual_patterns=MANUAL_PATTERNS, excel_path="Naming_Example.xlsx"
    )

    # Priority order: most specific → least specific.
    PATTERN_PRIORITY = [
        ["phosphatidylcholine", "phosphatidylethanolamine", "phosphatidylserine",
         "phosphatidylglycerol", "cardiolipin"],
        ["gpl_diacyl", "gpl_plasmalogen_full", "gpl_ether_chain", "gpl_dhap"],
        ["ceramide_nacyl_relaxed", "ceramide", "sphingomyelin",
         "ceramide_phosphoinositol", "sphingoid_base_free", "ganglioside_core",
         "glycosphingolipid_galactose", "glycosphingolipid_glucose",
         "acidic_glycosphingolipid_glucuronic", "sphingomyelin_detailed",
         "ceramide_phosphoethanolamine", "ceramide_phosphoinositol_detailed"],
        ["saccharolipid_n_acyl", "saccharolipid_o_acyl"],
        ["sterol_ester_generic", "sterol_core"],
        ["glycerolipid_o_glycosyl", "glycerolipid_monoacyl", "glycerolipid_acyl_ch",
         "glycosylglycerol_sn1", "diacylglycerol_1_2", "diacylglycerol_1_3",
         "monoacylglycerol_sn1", "monoacylglycerol_sn2", "triacylglycerol"],
        ["prenol_lactone", "prenol_allylic_ester", "prenol_polyene_ester",
         "prenol_amide", "fatty_acid", "plasmalogen"],
    ]

    def __init__(self):
        self.compiled_patterns: Dict[str, Tuple[HeadgroupPattern, Chem.Mol]] = {}
        self.compile_patterns()

    def compile_patterns(self):
        for pattern_id, pattern in self.HEADGROUP_PATTERNS.items():
            try:
                mol_pattern = Chem.MolFromSmarts(pattern.smarts)
                if mol_pattern is not None:
                    self.compiled_patterns[pattern_id] = (pattern, mol_pattern)
                else:
                    print(f"Warning: Failed to compile pattern '{pattern_id}': {pattern.name}",
                          file=sys.stderr)
            except Exception as e:
                print(f"Error compiling pattern '{pattern_id}': {e}", file=sys.stderr)

    def matches_pattern(self, mol: Chem.Mol, pattern_id: str) -> bool:
        if pattern_id not in self.compiled_patterns:
            return False

        pattern_info, mol_pattern = self.compiled_patterns[pattern_id]

        # Step 1: SMARTS substructure match
        if not mol.HasSubstructMatch(mol_pattern):
            return False

        # Step 2: FA chain guard
        if pattern_info.fa_positions:
            if not LipidAnalysis.has_long_carbon_chain(mol, min_len=6):
                return False

        # Step 3: Positional constraints (pattern-specific)
        if pattern_id == "saccharolipid_o_acyl":
            ring_info = mol.GetRingInfo()
            # Build atom → rings map for fast O-ring lookup
            ring_map: Dict[int, List[set]] = {}
            for ring in ring_info.AtomRings():
                for idx in ring:
                    ring_map.setdefault(idx, []).append(set(ring))

            def _in_O_ring(idx: int) -> bool:
                return any(
                    any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in r)
                    for r in ring_map.get(idx, [])
                )

            p_ring = Chem.MolFromSmarts("[C;R][OX2][CX3](=[OX1])[#6]~[#6]~[#6]~[#6]")
            p_exo  = Chem.MolFromSmarts("[CH2X4][OX2][CX3](=[OX1])[#6]~[#6]~[#6]~[#6]")

            o_ring_ok = False
            if p_ring:
                for match in mol.GetSubstructMatches(p_ring):
                    if _in_O_ring(match[0]):
                        o_ring_ok = True
                        break
            if not o_ring_ok and p_exo:
                for match in mol.GetSubstructMatches(p_exo):
                    if any(_in_O_ring(n.GetIdx())
                           for n in mol.GetAtomWithIdx(match[0]).GetNeighbors()):
                        o_ring_ok = True
                        break
            if not o_ring_ok:
                return False

            p_sph = Chem.MolFromSmarts("[NX3H2][CX4H1][CX4][OX2H1]")
            if p_sph and mol.HasSubstructMatch(p_sph):
                return False

        if pattern_id == "sterol_ester_generic":
            if mol.GetRingInfo().NumRings() < 3:
                return False
            p_phos = Chem.MolFromSmarts("[PX4]")
            if p_phos and mol.HasSubstructMatch(p_phos):
                return False

        return True

    @staticmethod
    def get_inchi(mol: Chem.Mol) -> str:
        return rdInchi.MolToInchi(mol) or ''

    def matches_any_valid_head(
        self,
        mol: Chem.Mol,
        reference_stereo_inchi: str = None,
        reference_inchi: str = None,
    ) -> bool:
        """
        Three modes depending on which reference is supplied:
        Headgroup class (no reference):
            Pure SMARTS topology check. Confirms the molecule belongs to a
            recognised lipid class. Stereo and FA position are ignored.

        Stereo-validated (reference_stereo_inchi):
            Headgroup class check AND InChI stereo layer comparison (/t /m /s).
            Rejects the molecule if its stereocentres differ from the reference.

        Complete identity (reference_inchi):
        """
        class_match = False
        for priority_group in self.PATTERN_PRIORITY:
            for pattern_id in priority_group:
                if self.matches_pattern(mol, pattern_id):
                    class_match = True
                    break
            if class_match:
                break

        if not class_match:
            manual_ids = {pid for group in self.PATTERN_PRIORITY for pid in group}
            for pattern_id in self.compiled_patterns:
                if pattern_id not in manual_ids:
                    if self.matches_pattern(mol, pattern_id):
                        class_match = True
                        break

        if not class_match:
            return False

        if reference_inchi:
            mol_inchi = LipidHeadValidator.get_inchi(mol)
            if not (mol_inchi == reference_inchi):
                return False

        if reference_stereo_inchi:
            mol_stereo = InChIParser.get_stereo_layer(mol)
            if mol_stereo != reference_stereo_inchi:
                return False

        return True

    def get_matching_patterns(self, mol: Chem.Mol) -> List[HeadgroupPattern]:
        for priority_group in self.PATTERN_PRIORITY:
            for pattern_id in priority_group:
                if self.matches_pattern(mol, pattern_id):
                    pattern_info, _ = self.compiled_patterns[pattern_id]
                    return [pattern_info]
        manual_ids = {pid for group in self.PATTERN_PRIORITY for pid in group}
        for pattern_id in self.compiled_patterns:
            if pattern_id not in manual_ids:
                if self.matches_pattern(mol, pattern_id):
                    pattern_info, _ = self.compiled_patterns[pattern_id]
                    return [pattern_info]
        return []

    def identify_lipid_class(self, mol: Chem.Mol) -> List[str]:
        """Return the lipid class of the best-matching pattern, or []."""
        matches = self.get_matching_patterns(mol)
        return sorted({m.lipid_class for m in matches}) if matches else []

    def validate_structure(self, mol: Chem.Mol, verbose: bool = False) -> Dict:
        """
        Full validation of a lipid structure.
        Returns:
            {
                "is_valid":        bool,
                "lipid_class":     str or None,
                "matched_patterns":List[HeadgroupPattern],
                "fa_positions":    List[str]
            }
        """
        matches = self.get_matching_patterns(mol)
        result = {
            "is_valid":         len(matches) > 0,
            "lipid_class":      matches[0].lipid_class if matches else None,
            "matched_patterns": matches,
            "fa_positions":     matches[0].fa_positions if matches else [],
        }
        if verbose and matches:
            print(f"Valid lipid: {result['lipid_class']} / {matches[0].name}", file=sys.stderr)
        return result

    @staticmethod
    def is_valid_lipid_structure(mol: Chem.Mol) -> bool:
        return LipidHeadValidator().matches_any_valid_head(mol)


# ─────────────────────────────────────────────────────────────────────────────

class LipidAnalysis:
    CLASSYFIRE_URL = "https://classyfire.wishartlab.com/entities.json"
    MIN_TAIL_CARBONS = 6

    HEAD_ANCHORS = {
        "carboxyl":             Chem.MolFromSmarts("C(=O)[O;H,-]"),
        "phosphate":            Chem.MolFromSmarts("P(=O)([O;H,-])[O;H,-]"),
        "phosphorylcholine":    Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C"),
        "phosphorylethanolamine": Chem.MolFromSmarts("OP(=O)([O-])OCCN"),
        "phosphorylserine":     Chem.MolFromSmarts("OP(=O)([O-])OC(N)C(=O)[O-]"),
        "phosphorylinositol":   Chem.MolFromSmarts("OP(=O)([O-])OC1C(O)C(O)C(O)C(O)C1O"),
        "sulfate":              Chem.MolFromSmarts("OS(=O)(=O)[O;H,-]"),
        "choline":              Chem.MolFromSmarts("OCC[N+](C)(C)C"),
        "ethanolamine":         Chem.MolFromSmarts("OCCN"),
        "serine":               Chem.MolFromSmarts("OC(N)C(=O)[O-]"),
    }

    @staticmethod
    def is_lipid_rdkit(mol):
        if mol is None:
            return False
        has_chain = LipidAnalysis.has_long_carbon_chain(mol, min_len=8)
        head_atoms = TailExtractor.detect_head_atoms(mol)
        ester, amide, ether = LipidAnalysis.count_lipid_linkages(mol)
        if len(head_atoms) > 0:
            tails = TailExtractor.extract_tails(mol)
            if len(tails) >= 1:
                return True
        if has_chain and mol.HasSubstructMatch(LipidAnalysis.HEAD_ANCHORS["carboxyl"]):
            return True
        if has_chain and amide >= 1:
            return True
        if has_chain and (ester >= 2 or ether >= 2):
            return True
        if mol.GetRingInfo().NumRings() >= 4 and has_chain:
            return True
        return False

    @staticmethod
    @lru_cache(maxsize=10000)
    def is_lipid_classyfire(inchi: str) -> bool:
        try:
            response = requests.get(
                LipidAnalysis.CLASSYFIRE_URL, params={"inchi": inchi}, timeout=5
            )
            if response.status_code != 200:
                return False
            data = response.json()
            if not data:
                return False
            entry = data[0]
            for field in ["kingdom", "superclass", "class", "subclass"]:
                name = entry.get(field, {}).get("name", "")
                if name and "lipid" in name.lower():
                    return True
            return False
        except Exception:
            return False

    @staticmethod
    def is_lipid(inchi: str, mol=None, use_classyfire=True) -> bool:
        if mol is None:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                return False
        if use_classyfire and LipidAnalysis.is_lipid_classyfire(inchi):
            return True
        return LipidAnalysis.is_lipid_rdkit(mol)

    @staticmethod
    def has_long_carbon_chain(mol, min_len=8):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            visited: set = set()
            stack = [(atom.GetIdx(), 0)]
            while stack:
                idx, length = stack.pop()
                if idx in visited:
                    continue
                visited.add(idx)
                a = mol.GetAtomWithIdx(idx)
                if a.GetAtomicNum() != 6:
                    continue
                length += 1
                if length >= min_len:
                    return True
                for nbr in a.GetNeighbors():
                    stack.append((nbr.GetIdx(), length))
        return False

    @staticmethod
    def count_lipid_linkages(mol):
        ester = amide = ether = 0
        for bond in mol.GetBonds():
            a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                if {a1.GetSymbol(), a2.GetSymbol()} == {"C", "O"}:
                    ester += 1
            if a1.GetSymbol() == "C" and a2.GetSymbol() == "N":
                amide += 1
            elif a2.GetSymbol() == "C" and a1.GetSymbol() == "N":
                amide += 1
            if bond.GetBondType() == Chem.BondType.SINGLE:
                if {a1.GetSymbol(), a2.GetSymbol()} == {"C", "O"}:
                    ether += 1
        return ester, amide, ether

    @staticmethod
    def tail_sig_levelB(t):
        return (t["C"], t["DB"], t["O"], t["DB_positions"], t["O_positions"])

    @staticmethod
    def tail_sig_levelC(t):
        return (t["C"], t["DB"], t["O"])

    @staticmethod
    def atom_count(mol):
        counts: Counter = Counter()
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] += 1
        return counts


class LipidComparator:
    @staticmethod
    def lipid_signature(mol):
        head_atoms = TailExtractor.detect_head_atoms(mol)
        if not head_atoms:
            return None
        tails = TailExtractor.extract_tails(mol)
        return tuple(sorted([
            (tail["C"], tail["DB"], tuple(sorted(tail["O_positions"])))
            for tail in tails
        ]))