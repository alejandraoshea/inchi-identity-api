import requests
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from collections import Counter

class TailExtractor:
    @staticmethod
    def detect_head_atoms(mol):
        head_atoms = set()

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in [7, 8, 15, 16]:  # N, O, P, S
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6:
                        head_atoms.add(atom.GetIdx())
                        break

        return head_atoms
    
    @staticmethod
    def extract_tails(mol):
        #Detect and exclude head atoms first
        head_atoms = TailExtractor.detect_head_atoms(mol)
        
        visited_global = set()
        visited_global.update(head_atoms)
        tails = []

        def walk_chain(start_atom_idx):
            visited = set()
            stack = [(start_atom_idx, None)]

            carbons = 0
            double_bonds = 0
            oxygens = 0
            db_positions = []
            o_positions = []

            position = 0

            while stack:
                atom_idx, parent = stack.pop()

                if atom_idx in visited or atom_idx in head_atoms:
                    continue

                visited.add(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)

                if atom.GetAtomicNum() == 6:
                    carbons += 1
                    position += 1

                    for bond in atom.GetBonds():
                        neighbor = bond.GetOtherAtom(atom)
                        neighbor_idx = neighbor.GetIdx()

                        if neighbor_idx == parent:
                            continue
                        
                        # stop at head atoms
                        if neighbor_idx in head_atoms:
                            continue

                        # continue carbon chain
                        if neighbor.GetAtomicNum() == 6:
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                double_bonds += 1
                                db_positions.append(position)

                            stack.append((neighbor_idx, atom_idx))

                        # oxygen attached (but not in head)
                        elif neighbor.GetAtomicNum() == 8:
                            oxygens += 1
                            o_positions.append(position)

            if carbons >= 6:  # filter small fragments
                return {
                    "C": carbons,
                    "DB": double_bonds,
                    "O": oxygens,
                    "DB_positions": tuple(sorted(db_positions)),
                    "O_positions": tuple(sorted(o_positions)),
                    "atoms": tuple(sorted(visited)),
                }

            return None

        # start from all carbons NOT in the headgroup
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                idx = atom.GetIdx()

                if idx in visited_global:
                    continue

                tail = walk_chain(idx)

                if tail:
                    tails.append(tail)
                    visited_global.update(tail["atoms"])

        # deduplicate by atoms
        unique = []
        seen_atoms = set()

        for t in tails:
            atom_key = t["atoms"]

            if atom_key not in seen_atoms:
                seen_atoms.add(atom_key)
                unique.append(t)

        return unique
    
    
    @staticmethod
    def walk_chain(start_atom, coming_from):
            stack = [(start_atom, None)]
            visited = set()

            carbons = 0
            double_bonds = 0
            oxygens = 0

            db_positions = []
            o_positions = []

            position = 0  # relative position along chain

            while stack:
                atom, prev = stack.pop()
                idx = atom.GetIdx()

                if idx in visited:
                    continue
                visited.add(idx)

                if atom.GetSymbol() == "C":
                    carbons += 1
                    position += 1
                elif atom.GetSymbol() == "O":
                    oxygens += 1
                    o_positions.append(position)

                for bond in atom.GetBonds():
                    nbr = bond.GetOtherAtom(atom)

                    if prev is not None and nbr.GetIdx() == prev:
                        continue

                    # detect C=C double bonds
                    if (
                        bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                        and atom.GetSymbol() == "C"
                        and nbr.GetSymbol() == "C"
                    ):
                        double_bonds += 1
                        db_positions.append(position)

                    if nbr.GetSymbol() in ["C", "O"]:
                        stack.append((nbr, atom.GetIdx()))

            return {
                "C": carbons,
                "DB": double_bonds,
                "O": oxygens,
                "DB_positions": tuple(sorted(db_positions)),
                "O_positions": tuple(sorted(o_positions)),
                "atoms": tuple(sorted(visited)),
            }
    