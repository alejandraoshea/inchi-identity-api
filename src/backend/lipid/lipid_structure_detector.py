from rdkit import Chem
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
    HEADGROUP_PATTERNS = {
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
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1H1-])[OX2][CH2X4][CH2X4][NX4+]([CH3X4])([CH3X4])[CH3X4])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylcholines",
            fa_positions=["sn-1", "sn-2"],
            description="PC headgroup with FAs at sn-1 and sn-2"
        ),
        
        #phosphatidylethanolamines (PE)
        "phosphatidylethanolamine": HeadgroupPattern(
            name="Phosphatidylethanolamine (PE)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1H1-])[OX2][CH2X4][CH2X4][NX3H2,NX4H3+])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylethanolamines",
            fa_positions=["sn-1", "sn-2"],
            description="PE headgroup with FAs at sn-1 and sn-2"
        ),
        
        #phosphatidylserines (PS)
        "phosphatidylserine": HeadgroupPattern(
            name="Phosphatidylserine (PS)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1H1-])[OX2][CH2X4][CHX4]([NX3H2,NX4H3+])[CX3](=[OX1])[OX2H0,OX1H1-])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylserines",
            fa_positions=["sn-1", "sn-2"],
            description="PS headgroup with FAs at sn-1 and sn-2"
        ),
        
        #phosphatidylglycerols (PG)
        "phosphatidylglycerol": HeadgroupPattern(
            name="Phosphatidylglycerol (PG)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1H1-])[OX2][CH2X4][CHX4]([OX2H1])[CH2X4][OX2H1])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Phosphatidylglycerols",
            fa_positions=["sn-1", "sn-2"],
            description="PG headgroup with FAs at sn-1 and sn-2"
        ),
        
        #cardiolipins (CL)
        "cardiolipin": HeadgroupPattern(
            name="Cardiolipin (CL)",
            smarts="[CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1H1-])[OX2][CH2X4][CHX4]([CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1H1-])[OX2][CH2X4][CHX4][CH2X4])[OX2])[OX2][CX3](=[OX1])[#6]",
            lipid_class="Cardiolipins",
            fa_positions=["sn-1", "sn-2", "sn-1'", "sn-2'"],
            description="CL with 4 FA chains"
        ),
        
        #cermides (Cer) - Sphingolipids
        "ceramide": HeadgroupPattern(
            name="Ceramide",
            smarts="[C@H](NC(=O)[#6])[C@H](O)[#6]",
            lipid_class="Ceramides",
            fa_positions=["N-acyl"],
            description="Sphingoid base with N-acyl FA"
        ),
        
        #sphingomyelins (SM)
        "sphingomyelin": HeadgroupPattern(
            name="Sphingomyelin (SM)",
            smarts="[NX3H1,NX3H2][CX3](=[OX1])[#6]",  # N-acyl amide linkage
            lipid_class="Sphingomyelins",
            fa_positions=["N-acyl"],
            description="SM - detects N-acyl FA on sphingoid base"
        ),

        "neutral_glycosphingolipid": HeadgroupPattern(
            name="Neutral glycosphingolipid",
            smarts="[CHX4]([NX3H1][CX3](=[OX1])[#6])[CHX4]=[CHX3]",  # Ceramide base
            lipid_class="Neutral glycosphingolipids",
            fa_positions=["N-acyl"],
            description="Glycosphingolipid with sugar residues"
        ),
        
        "acidic_glycosphingolipid": HeadgroupPattern(
            name="Acidic glycosphingolipid",
            smarts="[CHX4]([NX3H1][CX3](=[OX1])[#6])[CHX4]",  # Sphinganine base
            lipid_class="Acidic glycosphingolipids",
            fa_positions=["N-acyl"],
            description="Glycosphingolipid with acidic groups"
        ),
        
        "ceramide_phosphoinositol": HeadgroupPattern(
            name="Ceramide phosphoinositol",
            smarts="[NX3H1][CX3](=[OX1])[#6].[PX4](=[OX1])[OX2][CH]1[CH]([OX2H1])[CH]([OX2H1])[CH]([OX2H1])[CH]([OX2H1])[CH]1[OX2H1]",
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
            smarts="[CX3](=O)[OX2H1]",
            lipid_class="Fatty acids",
            fa_positions=[]
        ),

        "ganglioside_core": HeadgroupPattern(
            name="Ganglioside core",
            smarts="[NX3H1][CX3](=[OX1])[#6].[CH]1O[CH](CO)[CH](O)[CH](O)[CH]1O",
            lipid_class="Gangliosides",
            fa_positions=["N-acyl"],
            description="Ganglioside - requires ceramide base with sugar"
        ),
        
        #TODO: ADD MORE
    }
    
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
            True if molecule matches the pattern (FA in correct position)
        """
        if pattern_id not in self.compiled_patterns:
            return False
        
        _, mol_pattern = self.compiled_patterns[pattern_id]
        match = mol.HasSubstructMatch(mol_pattern)

        if not match:
            return False
        
        # 1. Fatty acid must REALLY exist
        if pattern_id in ["glycosylglycerol_sn1"]:
            if not self.has_fatty_acid_tail(mol):
                return False
            
        if pattern_id == "fatty_acid":
            if not self.has_fatty_acid_tail(mol):
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
    
    def has_fatty_acid_tail(self, mol: Chem.Mol, min_carbons: int = 6) -> bool:
        """Check if molecule has at least one long carbon chain (fatty acid-like)."""
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
                if length >= min_carbons:
                    return True

                for nbr in atom_obj.GetNeighbors():
                    stack.append((nbr.GetIdx(), length))

        return False
    
    def identify_lipid_class(self, mol: Chem.Mol) -> List[str]:  # Note: List[str] not Optional[str]
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
