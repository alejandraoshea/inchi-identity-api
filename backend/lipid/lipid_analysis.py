import requests
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from collections import Counter
from backend.inchi.lipid_tail_extraction import TailExtractor

class LipidAnalysis:
    CLASSYFIRE_URL = "https://classyfire.wishartlab.com/entities.json"

    #min chain length to consider it a FA
    MIN_TAIL_CARBONS = 6

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

            # check multiple taxonomy levels (more robust)
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
    