import requests
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from collections import Counter
from backend.lipid.lipid_tail_extraction import TailExtractor
from backend.lipid.lipid_pattern_generator import build_combined_patterns
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

@dataclass
class HeadgroupPattern:
    """
    Represents a lipid headgroup with positional specificity.
    Attributes:
        name: Human-readable name (e.g., "Glycosylglycerol Type A")
        smarts: SMARTS pattern that defines the headgroup AND FA attachment point(s)
        lipid_class: Classification (e.g., "Glycoglycerolipids")
        fa_positions: List of sn positions where FAs can attach (e.g., ["sn-1"])
        description: Additional details
    """
    name: str
    smarts: str
    lipid_class: str
    fa_positions: List[str]
    description: str = ""

class LipidHeadValidator:
    """
    Validates lipid structures using positional-specific SMARTS patterns.
    This class checks whether fatty acids are attached at the CORRECT positions
    on the lipid headgroup (sn-1, sn-2, sn-3, etc.).
    """
    MANUAL_PATTERNS = {
        # glycosylglycerols
        "glycosylglycerol_sn1": HeadgroupPattern(
            name="Glycosylglycerol (FA at sn-1)",
            smarts="O[C@H]1[C@H](OC[C@H](O)C[OX2][CX3](=O)[CX4,CX3]~[CX4,CX3]~[CX4,CX3])O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O",
            lipid_class="Glycoglycerolipids",
            fa_positions=["sn-1"],
            description="Glycosylglycerol"
        ),
        
        #diacylglycerols (DG)
        "diacylglycerol_1_2": HeadgroupPattern(
            name="1,2-Diacylglycerol",
            smarts="[OX2H1][CH2X4][CHX4]([OX2][CX3](=[OX1])[#6])[CH2X4][OX2][CX3](=[OX1])[#6]",
            lipid_class="Diacylglycerols",
            fa_positions=["sn-1", "sn-2"],
            description="FAs at positions 1 and 2 of glycerol backbone"
        ),
        
        "diacylglycerol_1_3": HeadgroupPattern(
            name="1,3-Diacylglycerol",
            smarts="[OX2][CX3](=[OX1])[#6].[CH2X4][CHX4]([OX2H1])[CH2X4][OX2][CX3](=[OX1])[#6]",
            lipid_class="Diacylglycerols",
            fa_positions=["sn-1", "sn-3"],
            description="FAs at positions 1 and 3 of glycerol backbone"
        ),
        
        #monoacylglyerols (MG)
        "monoacylglycerol_sn1": HeadgroupPattern(
            name="1-Monoacylglycerol",
            smarts="[OX2H1][CH2X4][CHX4]([OX2H1])[CH2X4][OX2][CX3](=[OX1])[#6]",
            lipid_class="Monoacylglycerols",
            fa_positions=["sn-1"],
            description="FA at position 1 of glycerol"
        ),
        
        "monoacylglycerol_sn2": HeadgroupPattern(
            name="2-Monoacylglycerol",
            smarts="[OX2H1][CH2X4][CHX4]([OX2][CX3](=[OX1])[#6])[CH2X4][OX2H1]",
            lipid_class="Monoacylglycerols",
            fa_positions=["sn-2"],
            description="FA at position 2 of glycerol"
        ),
        
        #triacylglycerols (TG)
        "triacylglycerol": HeadgroupPattern(
            name="Triacylglycerol",
            smarts="[CH2X4]([OX2][CX3](=[OX1])[#6])[CHX4]([OX2][CX3](=[OX1])[#6])[CH2X4][OX2][CX3](=[OX1])[#6]",
            lipid_class="Triacylglycerols",
            fa_positions=["sn-1", "sn-2", "sn-3"],
            description="FAs at all three positions of glycerol"
        ),
        
        #phosphatidylcholines (PC)
        "phosphatidylcholine": HeadgroupPattern(
            name="Phosphatidylcholine (PC)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CH2X4][NX4+]([CH3X4])([CH3X4])[CH3X4])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylcholines",
            fa_positions=["sn-1", "sn-2"],
            description="PC headgroup with FAs at sn-1 and sn-2"
        ),
        
        #phosphatidylethanolamines (PE)
        "phosphatidylethanolamine": HeadgroupPattern(
            name="Phosphatidylethanolamine (PE)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CH2X4][NX3H2,NX4H3+])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylethanolamines",
            fa_positions=["sn-1", "sn-2"],
            description="PE headgroup with FAs at sn-1 and sn-2"
        ),
        
        #phosphatidylserines (PS)
        "phosphatidylserine": HeadgroupPattern(
            name="Phosphatidylserine (PS)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CHX4]([NX3H2,NX4H3+])[CX3](=[OX1])[OX2H0,OX1-,OX1H1])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylserines",
            fa_positions=["sn-1", "sn-2"],
            description="PS headgroup with FAs at sn-1 and sn-2"
        ),
        
        #phosphatidylglycerols (PG)
        "phosphatidylglycerol": HeadgroupPattern(
            name="Phosphatidylglycerol (PG)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CHX4]([OX2H1])[CH2X4][OX2H1])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylglycerols",
            fa_positions=["sn-1", "sn-2"],
            description="PG headgroup with FAs at sn-1 and sn-2"
        ),
        
        #cardiolipins (CL)
        "cardiolipin": HeadgroupPattern(
            name="Cardiolipin (CL)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1-,OX1H1])[OX2][CH2X4][CHX4][CH2X4])[OX2])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Cardiolipins",
            fa_positions=["sn-1", "sn-2", "sn-1'", "sn-2'"],
            description="CL with 4 FA chains"
        ),
        
        #cermides (Cer) - Sphingolipids
        "ceramide": HeadgroupPattern(
            name="Ceramide",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]",
            lipid_class="Ceramides",
            fa_positions=["N-acyl"],
            description="Sphingoid base with N-acyl FA"
        ),
        
        # sphingomyelin — was matching GlcNAc acetamide
        "sphingomyelin": HeadgroupPattern(
            name="Sphingomyelin (SM)",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]",
            lipid_class="Sphingomyelins",
            fa_positions=["N-acyl"],
            description="SM - N-acyl FA on sphingoid base (8+ carbons required)"
        ),

        # ceramide_phosphoinositol — compound pattern, fix the ceramide half
        "ceramide_phosphoinositol": HeadgroupPattern(
            name="Ceramide phosphoinositol",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4].[PX4](=[OX1])[OX2][CH]1[CH]([OX2H1])[CH]([OX2H1])[CH]([OX2H1])[CH]([OX2H1])[CH]1[OX2H1]",
            lipid_class="Phosphosphingolipids",
            fa_positions=["N-acyl"],
            description="Ceramide with inositol phosphate"
        ),

        "sterol_core": HeadgroupPattern(
            name="Sterol core",
            smarts="[C@]12CC[C@H]3[C@@H](CC[C@]4(C)[C@@H]3CC=C4)C1CCC2",
            lipid_class="Sterols",
            fa_positions=[],
            description="Steroid nucleus"
        ),

        "plasmalogen": HeadgroupPattern(
            name="Plasmalogen",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][CH]=[CH])[OX2][PX4]",
            lipid_class="Ether phospholipids",
            fa_positions=["sn-1", "sn-2"],
        ),

        "fatty_acid": HeadgroupPattern(
            name="Fatty acid",
            smarts="[CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CX3](=O)[OX2H1]",
            lipid_class="Fatty acids",
            fa_positions=[]
        ),

        # ganglioside_core — fix the ceramide half of the compound pattern
        "ganglioside_core": HeadgroupPattern(
            name="Ganglioside core",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4].[CH]1O[CH](CO)[CH](O)[CH](O)[CH]1O",
            lipid_class="Gangliosides",
            fa_positions=["N-acyl"],
            description="Ganglioside - ceramide base with sugar (8+ carbon FA required)"
        ),
        
        "glycosphingolipid_galactose": HeadgroupPattern(
            name="Glycosphingolipid with Galactose",
            smarts="[H][C@](NC(CCCCCCCCCCCCCCC)=O)(CO[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC",
            lipid_class="Neutral glycosphingolipids",
            fa_positions=["N-acyl"],
            description="Neutral glycosphingolipid with galactose (Example 7)"
        ),
        
        "glycosphingolipid_glucose": HeadgroupPattern(
            name="Glycosphingolipid with Glucose",
            smarts="[H][C@](NC(CCCCCCCCCCCCCCC)=O)(CO[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@@](O)([H])/C=C/CCCCCCCCCCCCC",
            lipid_class="Neutral glycosphingolipids",
            fa_positions=["N-acyl"],
            description="Neutral glycosphingolipid with glucose"
        ),
        
        # acidic_glycosphingolipid_glucuronic — same fix
        "acidic_glycosphingolipid_glucuronic": HeadgroupPattern(
            name="Acidic Glycosphingolipid with Glucuronic acid",
            smarts="[NX3H1][CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]",
            lipid_class="Acidic glycosphingolipids",
            fa_positions=["N-acyl"],
            description="Acidic glycosphingolipid - requires N-acyl FA with 8+ carbons"
        ),
        
        "sphingomyelin_detailed": HeadgroupPattern(
            name="Sphingomyelin (detailed)",
            smarts="[H][C@](NC(CCCCCCCCCCCCCCCCCCC)=O)(COP(OCC[N+](C)(C)C)([O-])=O)[C@@](O)([H])CCCCCCCCCCCCCCC",
            lipid_class="Phosphosphingolipids",
            fa_positions=["N-acyl"],
            description="Sphingomyelin with phosphocholine (Example 64 from Excel)"
        ),
        
        "ceramide_phosphoethanolamine": HeadgroupPattern(
            name="Ceramide Phosphoethanolamine",
            smarts="[H][C@](NC(CCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(OCCN)(O)=O)[C@@](O)([H])CCCCCCCCCCCCCCC",
            lipid_class="Phosphosphingolipids",
            fa_positions=["N-acyl"],
            description="Ceramide phosphoethanolamine (Example 65 from Excel)"
        ),
        
        "ceramide_phosphoinositol_detailed": HeadgroupPattern(
            name="Ceramide Phosphoinositol (detailed)",
            smarts="[H][C@](NC(C(O)CCCCCCCCCCCCCCCCCCCCCCCC)=O)(COP(O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)(O)=O)[C@@](O)([H])CCCCCCCCCCCCCCCCC",
            lipid_class="Phosphosphingolipids",
            fa_positions=["N-acyl"],
            description="Ceramide with inositol phosphate (Example 66 from Excel)"
        ),
    }

    HEADGROUP_PATTERNS = build_combined_patterns(manual_patterns=MANUAL_PATTERNS, excel_path="Naming_Example.xlsx")
    
    def __init__(self):
        """Initialize the validator and compile all SMARTS patterns."""
        self.compiled_patterns: Dict[str, Tuple[HeadgroupPattern, Chem.Mol]] = {}
        self.compile_patterns()
    
    def compile_patterns(self):
        """Compile all SMARTS patterns for efficient matching."""
        for pattern_id, pattern in self.HEADGROUP_PATTERNS.items():
            try:
                mol_pattern = Chem.MolFromSmarts(pattern.smarts)
                if mol_pattern is not None:
                    self.compiled_patterns[pattern_id] = (pattern, mol_pattern)
                else:
                    print(f"Warning: Failed to compile pattern '{pattern_id}': {pattern.name}")
            except Exception as e:
                print(f"Error compiling pattern '{pattern_id}': {e}")
    
    def matches_pattern(self, mol: Chem.Mol, pattern_id: str) -> bool:
        """
        Check if a molecule matches a specific headgroup pattern.
        
        Args:
            mol: RDKit molecule to test
            pattern_id: ID of the pattern (key in HEADGROUP_PATTERNS)
        
        Returns:
            True if molecule matches the pattern
        """
        if pattern_id not in self.compiled_patterns:
            return False
        
        pattern_info, mol_pattern = self.compiled_patterns[pattern_id]
        
        # STEP 1: Check SMARTS match
        if not mol.HasSubstructMatch(mol_pattern):
            return False
        
        # STEP 2: If pattern expects FA(s), verify they exist
        if pattern_info.fa_positions:
            if not LipidAnalysis.has_long_carbon_chain(mol, min_len=6):
                return False
        
        return True
    
    def matches_any_valid_head(self, mol: Chem.Mol) -> bool:
        """
        Check if molecule matches ANY valid headgroup pattern.
        Args:
            mol: RDKit molecule to test
        Returns:
            True if molecule has a recognized lipid headgroup with FA in correct position
        """
        for pattern_id in self.compiled_patterns:
            if self.matches_pattern(mol, pattern_id):
                return True
        return False
    
    def get_matching_patterns(self, mol: Chem.Mol) -> List[HeadgroupPattern]:
        """
        Get all headgroup patterns that match this molecule.
        Args:
            mol: RDKit molecule to test
        
        Returns:
            List of matching HeadgroupPattern objects
        """
        matches = []
        for pattern_id, (pattern_info, mol_pattern) in self.compiled_patterns.items():
            if self.matches_pattern(mol, pattern_id):
                matches.append(pattern_info)
        return matches
    
    def identify_lipid_class(self, mol: Chem.Mol) -> List[str]: 
        """
        Identify the lipid class based on headgroup pattern.
        Args:
            mol: RDKit molecule to classify
        Returns:
            List of lipid class names (can be multiple if molecule matches multiple patterns)
        """
        matches = self.get_matching_patterns(mol)
        return sorted(set(m.lipid_class for m in matches)) if matches else []
    
    
    def validate_structure(self, mol: Chem.Mol, verbose: bool = False) -> Dict:
        """
        Complete validation of a lipid structure.
        Args:
            mol: RDKit molecule to validate
            verbose: If True, return detailed information
        Returns:
            Dictionary with validation results:
            {
                "is_valid": bool,
                "lipid_class": str or None,
                "matched_patterns": List[HeadgroupPattern],
                "fa_positions": List[str]
            }
        """
        matches = self.get_matching_patterns(mol)
        
        result = {
            "is_valid": len(matches) > 0,
            "lipid_class": matches[0].lipid_class if matches else None,
            "matched_patterns": matches,
            "fa_positions": matches[0].fa_positions if matches else []
        }
        
        if verbose and matches:
            print(f"Valid lipid structure detected:")
            print(f"Class: {result['lipid_class']}")
            print(f"Pattern: {matches[0].name}")
            print(f"FA positions: {', '.join(result['fa_positions'])}")
        
        return result

    def is_valid_lipid_structure(mol: Chem.Mol) -> bool:
        """
        Quick check: Does this molecule have a valid lipid headgroup
        with FA attached at the CORRECT position?
        
        Args:
            mol: RDKit molecule
        
        Returns:
            True if valid lipid structure
        """
        validator = LipidHeadValidator()
        return validator.matches_any_valid_head(mol)

class LipidAnalysis:
    CLASSYFIRE_URL = "https://classyfire.wishartlab.com/entities.json"

    #min chain length to consider it a FA
    MIN_TAIL_CARBONS = 6

    HEAD_ANCHORS = {
        "carboxyl": Chem.MolFromSmarts("C(=O)[O;H,-]"),
        "phosphate": Chem.MolFromSmarts("P(=O)([O;H,-])[O;H,-]"),
        "phosphorylcholine": Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C"),
        "phosphorylethanolamine": Chem.MolFromSmarts("OP(=O)([O-])OCCN"),
        "phosphorylserine": Chem.MolFromSmarts("OP(=O)([O-])OC(N)C(=O)[O-]"),
        "phosphorylinositol": Chem.MolFromSmarts("OP(=O)([O-])OC1C(O)C(O)C(O)C(O)C1O"),
        "sulfate": Chem.MolFromSmarts("OS(=O)(=O)[O;H,-]"),
        "choline": Chem.MolFromSmarts("OCC[N+](C)(C)C"),
        "ethanolamine": Chem.MolFromSmarts("OCCN"),
        "serine": Chem.MolFromSmarts("OC(N)C(=O)[O-]")
    }


    @staticmethod
    def is_lipid_rdkit(mol):
        if mol is None:
            return False

        has_chain = LipidAnalysis.has_long_carbon_chain(mol, min_len=8)
        head_atoms = TailExtractor.detect_head_atoms(mol)
        ester, amide, ether = LipidAnalysis.count_lipid_linkages(mol)

        # Strong lipid signal
        if len(head_atoms) > 0:
            tails = TailExtractor.extract_tails(mol)
            if len(tails) >= 1:
                return True

        # Fatty acids
        if has_chain and mol.HasSubstructMatch(LipidAnalysis.HEAD_ANCHORS["carboxyl"]):
            return True

        # Sphingolipids
        if has_chain and amide >= 1:
            return True

        # Glycerolipids
        if has_chain and (ester >= 2 or ether >= 2):
            return True

        # Sterols
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() >= 4 and has_chain:
            return True

        return False
    

    @staticmethod
    @lru_cache(maxsize=10000)
    def is_lipid_classyfire(inchi: str) -> bool:
        try:
            response = requests.get(
                LipidAnalysis.CLASSYFIRE_URL,
                params={"inchi": inchi},
                timeout=5
            )

            if response.status_code != 200:
                return False

            data = response.json()
            if not data:
                return False

            entry = data[0]

            taxonomy_fields = [
                entry.get("kingdom", {}).get("name", ""),
                entry.get("superclass", {}).get("name", ""),
                entry.get("class", {}).get("name", ""),
                entry.get("subclass", {}).get("name", "")
            ]

            for field in taxonomy_fields:
                if field and "lipid" in field.lower():
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

        # STEP 1: Classyfire
        if use_classyfire:
            if LipidAnalysis.is_lipid_classyfire(inchi):
                return True

        # STEP 2: fallback RDKit
        return LipidAnalysis.is_lipid_rdkit(mol)
    

    @staticmethod
    def has_long_carbon_chain(mol, min_len=8):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue

            visited = set()
            stack = [(atom.GetIdx(), 0)]

            while stack:
                idx, length = stack.pop()
                if idx in visited:
                    continue

                visited.add(idx)
                atom_obj = mol.GetAtomWithIdx(idx)

                if atom_obj.GetAtomicNum() != 6:
                    continue

                length += 1
                if length >= min_len:
                    return True

                for nbr in atom_obj.GetNeighbors():
                    stack.append((nbr.GetIdx(), length))

        return False
    
    @staticmethod
    def count_lipid_linkages(mol):
        ester = 0
        amide = 0
        ether = 0

        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()

            # Ester: C=O connected to O-C
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                if {a1.GetSymbol(), a2.GetSymbol()} == {"C", "O"}:
                    ester += 1

            # Amide
            if a1.GetSymbol() == "C" and a2.GetSymbol() == "N":
                amide += 1
            elif a2.GetSymbol() == "C" and a1.GetSymbol() == "N":
                amide += 1

            # Ether
            if bond.GetBondType() == Chem.BondType.SINGLE:
                if {a1.GetSymbol(), a2.GetSymbol()} == {"C", "O"}:
                    ether += 1

        return ester, amide, ether

    def tail_sig_levelB(t):
        return (
            t["C"],
            t["DB"],
            t["O"],
            t["DB_positions"],
            t["O_positions"],
        )

    def tail_sig_levelC(t):
        return (
            t["C"],
            t["DB"],
            t["O"],
        )

    def atom_count(mol):
        counts = Counter()
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] += 1
        return counts
    
class LipidComparator:
    @staticmethod
    def lipid_signature(mol):
        """
        Generate a signature for lipid comparison.
        """
        head_atoms = TailExtractor.detect_head_atoms(mol)

        if not head_atoms:
            return None

        tails = TailExtractor.extract_tails(mol)
        
        return tuple(
            sorted([
                (tail["C"], tail["DB"], tuple(sorted(tail["O_positions"])))
                for tail in tails
            ])
        )