import requests
from functools import lru_cache
from rdkit import Chem
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

            # Check multiple taxonomy levels (more robust)
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

    @staticmethod
    def lipid_identity_levels(mol1, mol2):
        """
        Returns lipid identity levels
        Level A: Cis/trans y resto idéntico
        Level B: Posicion cadenas
        Level C: Posición dobles enlaces y oxígeno
        Level D: Número total de carbonos, dobles enlaces y oxígenos
        """

        # --- STEP A: remove cis/trans ---
        mol1_A = Chem.Mol(mol1)
        mol2_A = Chem.Mol(mol2)

        LipidAnalysis.remove_cis_trans(mol1_A)
        LipidAnalysis.remove_cis_trans(mol2_A)

        # IMPORTANT: extract from modified mols
        tails1 = LipidAnalysis.extract_detailed_tails(mol1_A)
        tails2 = LipidAnalysis.extract_detailed_tails(mol2_A)

        # --- LEVEL A ---
        LEVELA = tails1 == tails2

        # --- LEVEL B ---
        LEVELB = sorted(tails1) == sorted(tails2)

        # --- LEVEL C ---
        sig1 = sorted([(t["C"], t["DB"], t["O"]) for t in tails1])
        sig2 = sorted([(t["C"], t["DB"], t["O"]) for t in tails2])

        LEVELC = sig1 == sig2

        # --- LEVEL D ---
        total1 = (
            sum(t["C"] for t in tails1),
            sum(t["DB"] for t in tails1),
            sum(t["O"] for t in tails1),
        )

        total2 = (
            sum(t["C"] for t in tails2),
            sum(t["DB"] for t in tails2),
            sum(t["O"] for t in tails2),
        )

        LEVELD = total1 == total2

        return {
            "LEVELA": LEVELA,
            "LEVELB": LEVELB,
            "LEVELC": LEVELC,
            "LEVELD": LEVELD
        }

    # STEP 1 remove cis/trans stereochemistry
    @staticmethod
    def remove_cis_trans(mol):
        for bond in mol.GetBonds():
            if (
                bond.GetBondType() == Chem.BondType.DOUBLE
                and bond.GetBeginAtom().GetAtomicNum() == 6
                and bond.GetEndAtom().GetAtomicNum() == 6):
                    bond.SetStereo(Chem.BondStereo.STEREONONE) #rdkit recomputes stereochemistry 
                
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

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

    # STEP 3 extract tails via graph traversal
    @staticmethod
    def extract_tails(mol, head_atoms):
        """
        Extract lipid hydrocarbon tails.
        1. Start from carbon atoms directly connected to the polar head.
        2. Traverse only through carbon atoms.
        3. Stop if we return to head atoms.
        """
        
        visited = set(head_atoms)
        tails = []

        for head_idx in head_atoms:
            head_atom = mol.GetAtomWithIdx(head_idx)

            for nbr in head_atom.GetNeighbors():
                # tails must start with carbon
                if nbr.GetAtomicNum() != 6:
                    continue

                start_tail = nbr.GetIdx()
                if start_tail in visited:
                    continue

                stack = [start_tail]
                tail_atoms = set()

                while stack:
                    current_idx = stack.pop()
                    if current_idx in visited or current_idx in head_atoms:
                        continue

                    atom = mol.GetAtomWithIdx(current_idx)
                    # stop if not carbon (prevents entering headgroup)
                    if atom.GetAtomicNum() != 6:
                        continue
                    
                    visited.add(current_idx)
                    tail_atoms.add(current_idx)

                    for next_atom in atom.GetNeighbors():
                        stack.append(next_atom.GetIdx())

                if len(tail_atoms) >= LipidAnalysis.MIN_TAIL_CARBONS:
                    tails.append(tail_atoms)

        return tails

    @staticmethod
    def extract_detailed_tails(mol):
        head_atoms = LipidAnalysis.detect_head_atoms(mol)

        if not head_atoms:
            return []

        tails = LipidAnalysis.extract_tails(mol, head_atoms)

        detailed = []

        for tail in tails:
            carbons = 0
            double_bonds = 0
            oxygens = 0

            for idx in tail:
                atom = mol.GetAtomWithIdx(idx)

                if atom.GetAtomicNum() == 6:
                    carbons += 1
                elif atom.GetAtomicNum() == 8:
                    oxygens += 1

                for bond in atom.GetBonds():
                    if (
                        bond.GetBondType() == Chem.BondType.DOUBLE and
                        bond.GetBeginAtom().GetAtomicNum() == 6 and
                        bond.GetEndAtom().GetAtomicNum() == 6
                    ):
                        double_bonds += 1

            double_bonds //= 2

            if carbons >= LipidAnalysis.MIN_TAIL_CARBONS:
                detailed.append({
                    "C": carbons,
                    "DB": double_bonds,
                    "O": oxygens
                })

        return detailed

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

    @staticmethod
    def lipid_chain_signatures(mol):
            #extracts substituets
            chains = []

            for bond in mol.GetBonds():
                #examine every bond
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()

                # detect ester C=O and we identify the carbonyl atom
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    if atom1.GetSymbol() == "C" and atom2.GetSymbol() == "O":
                        carbon = atom1
                    elif atom2.GetSymbol() == "C" and atom1.GetSymbol() == "O":
                        carbon = atom2

                    else:
                        continue

                    # follow carbon (FA) chain
                    for nbr in carbon.GetNeighbors():
                        if nbr.GetSymbol() == "C":
                            #generate the canonical SMILES representation of substituents
                            chain = Chem.MolFragmentToSmiles(
                                mol,
                                atomsToUse=[nbr.GetIdx()],
                                canonical=True,
                                isomericSmiles=False
                            )

                            chains.append(chain)

                #ether detection
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    if atom1.GetSymbol() == "O" and atom2.GetSymbol() == "C":
                        if atom1.GetDegree() > 1:
                            for nbr in atom2.GetNeighbors():
                                if nbr.GetSymbol() == "C" and nbr.GetIdx() != atom1.GetIdx():
                                    chain = Chem.MolFragmentToSmiles(
                                        mol,
                                        atomsToUse=[nbr.GetIdx()],
                                        canonical=True,
                                        isomericSmiles=False
                                    )
                                    if chain.count("C") >= 8:
                                        chains.append(chain)

                    elif atom2.GetSymbol() == "O" and atom1.GetSymbol() == "C":
                        if atom1.GetDegree() > 1:
                            for nbr in atom1.GetNeighbors():
                                if nbr.GetSymbol() == "C" and nbr.GetIdx() != atom2.GetIdx():
                                    chain = Chem.MolFragmentToSmiles(
                                        mol,
                                        atomsToUse=[nbr.GetIdx()],
                                        canonical=True,
                                        isomericSmiles=False
                                    )

                                if chain.count("C") >= 8:
                                    chains.append(chain)

            return Counter(chains)

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
    