from rdkit import Chem


class LipidAnalysis:

    # SMARTS patterns used to detect polar lipid headgroups
    HEAD_ANCHORS = {
        "phosphate": Chem.MolFromSmarts("P(=O)(O)(O)"), #detects phospolipids
        "amine": Chem.MolFromSmarts("[NX3;H2,H1,H0;+0,+1]"), #detects nitrogen group like ethanolamine, choline, serine
        "carboxyl": Chem.MolFromSmarts("C(=O)[O;H,-]"), #detects FA headgroups
        "amide": Chem.MolFromSmarts("C(=O)N"), #detects sphingolipid linkages
        "sugar": Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]1O"), #detects glycolipid heads
    }

    #min chain length to consider it a FA
    MIN_TAIL_CARBONS = 6

    # STEP 1 detect head atoms via SMARTS anchors
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

    # STEP 2 extract tails via graph traversal
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

    # STEP 3 compute tail signature
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

    # STEP 4: build lipid signature
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

    # STEP 5 remove cis/trans stereochemistry
    @staticmethod
    def remove_cis_trans(mol):
        for bond in mol.GetBonds():
            if (
                bond.GetBondType() == Chem.BondType.DOUBLE
                and bond.GetBeginAtom().GetAtomicNum() == 6
                and bond.GetEndAtom().GetAtomicNum() == 6):
                    bond.SetStereo(Chem.BondStereo.STEREONONE) #rdkit recomputes stereochemistry 
                
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # STEP 6 compare molecules
    @staticmethod
    def equal_ignore_double_bond_position(mol1, mol2):
        #FA : lipid comparison
        #other: normal structure comparison
        sig1 = LipidAnalysis.lipid_signature(mol1)
        sig2 = LipidAnalysis.lipid_signature(mol2)

        # if not a fatty-acid lipid
        if sig1 is None or sig2 is None:
            return Chem.MolToSmiles(mol1, canonical=True) == Chem.MolToSmiles(mol2, canonical=True)

        return sig1 == sig2