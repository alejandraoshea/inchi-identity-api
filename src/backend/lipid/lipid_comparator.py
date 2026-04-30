import requests
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from collections import Counter
from backend.lipid.lipid_tail_extraction import TailExtractor
from backend.lipid.lipid_comparator import LipidComparator

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
        
        # use already-computed values - no recomputation needed
        return tuple(
            sorted([
                (tail["C"], tail["DB"], tuple(sorted(tail["O_positions"])))
                for tail in tails
            ])
        )
    
    @staticmethod
    def equal_ignore_double_bond_position(mol1, mol2):
        sig1 = LipidComparator.lipid_signature(mol1)
        sig2 = LipidComparator.lipid_signature(mol2)

        if sig1 is None or sig2 is None:
            return Chem.MolToSmiles(mol1, canonical=True) == Chem.MolToSmiles(mol2, canonical=True)

        return sig1 == sig2