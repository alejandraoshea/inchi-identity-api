from rdkit import Chem


class LipidAnalysis:

    # generic polar head anchors 
    HEAD_ANCHORS = {
        "phosphate": Chem.MolFromSmarts("P(=O)(O)(O)"),
        "amine": Chem.MolFromSmarts("[NX3;H2,H1,H0;+0,+1]"),
        "carboxyl": Chem.MolFromSmarts("C(=O)[O;H,-]"),
        "amide": Chem.MolFromSmarts("C(=O)N"),
        "sugar": Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]1O"),
    }

    MIN_TAIL_CARBONS = 6

    # STEP 1 remove cis/trans stereochemistry
    @staticmethod
    def remove_cis_trans(mol):
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                bond.SetStereo(Chem.BondStereo.STEREONONE)

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # STEP 2 detect head atoms via SMARTS anchors
    @staticmethod
    def detect_head_atoms(mol):
        head_atoms = set()
        for smarts in LipidAnalysis.HEAD_ANCHORS.values():
            matches = mol.GetSubstructMatches(smarts)

            for match in matches:
                head_atoms.update(match)

        return head_atoms

    # STEP 3 extract tails via graph traversal
    @staticmethod
    def extract_tails(mol, head_atoms):
        visited = set(head_atoms)
        tails = []

        for atom_idx in head_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)

            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()

                if nbr_idx in visited:
                    continue

                stack = [nbr_idx]
                tail_atoms = set()

                while stack:
                    current = stack.pop()
                    if current in visited:
                        continue

                    visited.add(current)
                    tail_atoms.add(current)
                    current_atom = mol.GetAtomWithIdx(current)

                    for next_atom in current_atom.GetNeighbors():
                        next_idx = next_atom.GetIdx()

                        if next_idx not in visited:
                            stack.append(next_idx)

                if tail_atoms:
                    tails.append(tail_atoms)

        return tails

    # STEP 4 compute tail signature
    @staticmethod
    def tail_signature(mol, tail_atoms):
        carbons = 0
        double_bonds = 0

        for idx in tail_atoms:
            atom = mol.GetAtomWithIdx(idx)

            if atom.GetAtomicNum() == 6:
                carbons += 1

            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_bonds += 1

        double_bonds = double_bonds // 2

        if carbons < LipidAnalysis.MIN_TAIL_CARBONS:
            return None

        return (carbons, double_bonds)

    # STEP 5: build lipid signature
    @staticmethod
    def lipid_signature(mol):
        LipidAnalysis.remove_cis_trans(mol)
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
        #FA : lipid comparison
        #other: normal structure comparison
        sig1 = LipidAnalysis.lipid_signature(mol1)
        sig2 = LipidAnalysis.lipid_signature(mol2)

        # if not a fatty-acid lipid
        if sig1 is None or sig2 is None:
            return Chem.MolToSmiles(mol1, canonical=True) == Chem.MolToSmiles(mol2, canonical=True)

        return sig1 == sig2