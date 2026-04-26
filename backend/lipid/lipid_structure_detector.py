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
            smarts="O[C@H]1[C@H](OC[C@H](O)C[OX2][CX3](=O)[#6])O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O",
            lipid_class="Glycoglycerolipids",
            fa_positions=["sn-1"],
            description="Your teacher's example - FA attached at the correct glycerol position"
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
            smarts="[CH2X4][OX2H1].[CHX4]([NX3H1][CX3](=[OX1])[#6])[CHX4]=[CHX3]",
            lipid_class="Ceramides",
            fa_positions=["N-acyl"],
            description="Sphingoid base with N-acyl FA"
        ),
        
        #sphingomyelins (SM)
        "sphingomyelin": HeadgroupPattern(
            name="Sphingomyelin (SM)",
            smarts="[CH2X4][OX2][PX4](=[OX1])([OX2H0,OX1H1-])[OX2][CH2X4][CH2X4][NX4+]([CH3X4])([CH3X4])[CH3X4].[CHX4]([NX3H1][CX3](=[OX1])[#6])[CHX4]=[CHX3]",
            lipid_class="Sphingomyelins",
            fa_positions=["N-acyl"],
            description="SM with phosphocholine headgroup"
        ),
        
        #TODO: ADD MORE
    }
    
    def __init__(self):
        """Initialize the validator and compile all SMARTS patterns."""
        self._compiled_patterns: Dict[str, Tuple[HeadgroupPattern, Chem.Mol]] = {}
        self._compile_patterns()
    
    def _compile_patterns(self):
        """Compile all SMARTS patterns for efficient matching."""
        for pattern_id, pattern in self.HEADGROUP_PATTERNS.items():
            try:
                mol_pattern = Chem.MolFromSmarts(pattern.smarts)
                if mol_pattern is not None:
                    self._compiled_patterns[pattern_id] = (pattern, mol_pattern)
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
        if pattern_id not in self._compiled_patterns:
            return False
        
        _, mol_pattern = self._compiled_patterns[pattern_id]
        return mol.HasSubstructMatch(mol_pattern)
    
    def matches_any_valid_head(self, mol: Chem.Mol) -> bool:
        """
        Check if molecule matches ANY valid headgroup pattern.
        Args:
            mol: RDKit molecule to test
        Returns:
            True if molecule has a recognized lipid headgroup with FA in correct position
        """
        for pattern_id in self._compiled_patterns:
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
        for pattern_id, (pattern_info, mol_pattern) in self._compiled_patterns.items():
            if mol.HasSubstructMatch(mol_pattern):
                matches.append(pattern_info)
        return matches
    
    def identify_lipid_class(self, mol: Chem.Mol) -> Optional[str]:
        """
        Identify the lipid class based on headgroup pattern.
        Args:
            mol: RDKit molecule to classify
        Returns:
            Lipid class name (e.g., "Glycoglycerolipids") or None
        """
        matches = self.get_matching_patterns(mol)
        if matches:
            return matches[0].lipid_class
        return None
    
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


def get_lipid_class(mol: Chem.Mol) -> Optional[str]:
    """
    Identify the lipid class of a molecule.
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Lipid class name or None
    """
    validator = LipidHeadValidator()
    return validator.identify_lipid_class(mol)
