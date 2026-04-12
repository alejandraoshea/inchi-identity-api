import requests
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from collections import Counter

class LipidAnalysis:
    CLASSYFIRE_URL = "https://classyfire.wishartlab.com/entities.json"

    # SMARTS patterns used to detect polar lipid headgroups
    HEAD_ANCHORS = {
        "phosphate": Chem.MolFromSmarts("P(=O)(O)(O)"), #detects phospolipids
        "phosphocholine": Chem.MolFromSmarts("P(=O)(O)(O)OCC[N+](C)(C)C"),
        "phosphoethanolamine": Chem.MolFromSmarts("P(=O)(O)(O)OCCN"),
        "phosphoserine": Chem.MolFromSmarts("P(=O)(O)(O)OCC(N)C(=O)O"),
        "amine": Chem.MolFromSmarts("[NX3;H2,H1,H0;+0,+1]"), #detects nitrogen group like ethanolamine, choline, serine
        "quaternary_amine": Chem.MolFromSmarts("[NX4+]"),
        "amide": Chem.MolFromSmarts("C(=O)N"), #detects sphingolipid linkages
        "sugar_ring": Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O"), #detects glycolipid heads
        "glycerol": Chem.MolFromSmarts("OCC(O)CO"),
        "carboxyl": Chem.MolFromSmarts("C(=O)[O;H,-]"), #detects FA headgroups
        "sulfonate": Chem.MolFromSmarts("S(=O)(=O)O"),
        "sterol_core": Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCCC4"),
    }

    #min chain length to consider it a FA
    MIN_TAIL_CARBONS = 6

    @staticmethod
    def is_lipid_rdkit(mol):
        if mol is None:
            return False

        has_chain = LipidAnalysis.has_long_carbon_chain(mol, min_len=8)
        head_atoms = LipidAnalysis.detect_head_atoms(mol)
        ester, amide, ether = LipidAnalysis.count_lipid_linkages(mol)

        # Strong lipid signal
        if has_chain and len(head_atoms) > 0:
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

    # STEP 1 remove cis/trans stereochemistry
    @staticmethod
    def remove_cis_trans(mol):
        mol = Chem.Mol(mol)
        for bond in mol.GetBonds():
            if (
                bond.GetBondType() == Chem.BondType.DOUBLE
                and bond.GetBeginAtom().GetAtomicNum() == 6
                and bond.GetEndAtom().GetAtomicNum() == 6):
                    bond.SetStereo(Chem.BondStereo.STEREONONE) #rdkit recomputes stereochemistry 
                
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        return mol

    # STEP 2 detect head atoms via SMARTS anchors
    @staticmethod
    def detect_head_atoms(mol):
        """
        Find atoms belonging to the lipid headgroup
        For each SMARTS pattern: find matches and store atoms
        """
        head_atoms = set()
        for smarts in LipidAnalysis.HEAD_ANCHORS.values():
            matches = mol.GetSubstructMatches(smarts)

            for match in matches:
                head_atoms.update(match)

        return head_atoms

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
    
    # STEP 3 extract tails via graph traversal - FIXED VERSION
    @staticmethod
    def extract_tails(mol):
        # CRITICAL: Detect and exclude head atoms first
        head_atoms = LipidAnalysis.detect_head_atoms(mol)
        
        visited_global = set()
        visited_global.update(head_atoms)  # Don't start chain walks from head atoms
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
                        
                        # Stop at head atoms
                        if neighbor_idx in head_atoms:
                            continue

                        # Continue carbon chain
                        if neighbor.GetAtomicNum() == 6:
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                double_bonds += 1
                                db_positions.append(position)

                            stack.append((neighbor_idx, atom_idx))

                        # Oxygen attached (but not in head)
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

        # Start from all carbons NOT in the headgroup
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                idx = atom.GetIdx()

                if idx in visited_global:
                    continue

                tail = walk_chain(idx)

                if tail:
                    tails.append(tail)
                    visited_global.update(tail["atoms"])

        # Deduplicate by atoms
        unique = []
        seen_atoms = set()

        for t in tails:
            atom_key = t["atoms"]

            if atom_key not in seen_atoms:
                seen_atoms.add(atom_key)
                unique.append(t)

        return unique
    
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
    
    # STEP 4 compute tail signature
    @staticmethod
    def tail_signature(mol, tail_atoms):
        """
        Converts a tail into a simple numeric signature:
        CH3-(CH2)16-CH=CH-COOH into (18 carbons, 1 double bond)
        """
        carbons = 0
        double_bonds = 0

        for idx in tail_atoms:
            atom = mol.GetAtomWithIdx(idx)

            #we count C atoms
            if atom.GetAtomicNum() == 6:
                carbons += 1

            for bond in atom.GetBonds():
                #we count double bonds
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_bonds += 1

        #each bond appears twice in transversal
        double_bonds = double_bonds // 2

        #discard fragments that are too short
        if carbons < LipidAnalysis.MIN_TAIL_CARBONS:
            return None

        return (carbons, double_bonds)

    # STEP 5: build lipid signature
    @staticmethod
    def lipid_signature(mol):
        head_atoms = LipidAnalysis.detect_head_atoms(mol)

        if not head_atoms:
            return None

        tails = LipidAnalysis.extract_tails(mol, head_atoms)
        tail_sigs = []

        for tail in tails:
            sig = LipidAnalysis.tail_signature(mol, tail)
            if sig:
                tail_sigs.append(sig)

        return tuple(sorted(tail_sigs))

    # STEP 6 compare molecules
    @staticmethod
    def equal_ignore_double_bond_position(mol1, mol2):
        sig1 = LipidAnalysis.lipid_signature(mol1)
        sig2 = LipidAnalysis.lipid_signature(mol2)

        if sig1 is None or sig2 is None:
            return Chem.MolToSmiles(mol1, canonical=True) == Chem.MolToSmiles(mol2, canonical=True)

        return sig1 == sig2