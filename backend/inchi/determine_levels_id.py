from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdmolops
from rdkit.Chem import inchi
from rdkit.Chem.SaltRemover import SaltRemover
from backend.inchi.inchi_parser import InChiParser
from backend.inchi.inchi_layers_enum import InchiLayers
from backend.inchi.lipid_analysis import LipidAnalysis
import subprocess, os
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import Counter
from rdkit.Chem.MolStandardize import rdMolStandardize

class InChi:
    def isCompleteIdentity(inchi1: str, inchi2: str) -> bool:
        return inchi1 == inchi2

    def mol_from_inchi(inchi: str):
        try:
            mol = Chem.MolFromInchi(inchi.strip())
            if mol is None:
                raise ValueError(f"Invalid InChI: {inchi}")
            return mol
        except Exception as e:
            print(f"RDKit conversion failed for InChI: {inchi}\nError: {e}")
            return None


    def main_fragment(mol):
        try:
            frags = Chem.GetMolFrags(mol, asMols=True)
            # pick largest fragment (most atoms)
            main = max(frags, key=lambda m: m.GetNumAtoms())
            Chem.SanitizeMol(main)
            Chem.AssignStereochemistry(main, cleanIt=True, force=True)
            return main
        except Exception:
            return mol
    
    def areEqualNoIsotopes(inchi1: str, inchi2: str) -> bool:
        inchi1_isotopes = InChiParser.removeIsotopicLayers(inchi1)
        inchi2_isotopes = InChiParser.removeIsotopicLayers(inchi2)
        return inchi1_isotopes == inchi2_isotopes
    
    def areEqualDisolvedSalts(inchi1: str, inchi2: str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)

        if not mol1 or not mol2:
            return False    

        main1 = InChi.main_fragment(mol1)
        main2 = InChi.main_fragment(mol2)

        try:
            inchi_main1 = inchi.MolToInchi(main1)
            inchi_main2 = inchi.MolToInchi(main2)
        except Exception:
            return False

        return inchi_main1 == inchi_main2
    
    #helper method: detect negative charge to neutralize the mol
    def has_negative_charge(inchi: str) -> bool:
        for part in inchi.split("/"):
            if part.startswith("q-") or part.startswith("p-"):
                return True
        return False
    
    def get_charge_info(inchi: str):
        has_p_plus = False
        has_p_minus = False
        has_q_plus = False
        has_q_minus = False

        for part in inchi.split("/"):
            if part.startswith("p+"):
                has_p_plus = True
            elif part.startswith("p-"):
                has_p_minus = True
            elif part.startswith("q+"):
                has_q_plus = True
            elif part.startswith("q-"):
                has_q_minus = True

        return has_p_plus, has_p_minus, has_q_plus, has_q_minus
    
    def remove_only_p_layer(inchi: str) -> str:
        parts = []
        for p in inchi.strip().split("/"):
            if p.startswith("p"):
                continue
            if p:  # avoid empty segments
                parts.append(p.strip())
        return "/".join(parts)
  
    def neutralize_molecule(mol):
        # standardize and cleanup the molecule
        clean_mol = rdMolStandardize.Cleanup(mol)
        # remove fragments and keep largest, then neutralize
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
        uncharger = rdMolStandardize.Uncharger()
        return uncharger.uncharge(parent_clean_mol)

    def areEqualNoCharges(inchi1: str, inchi2: str) -> bool:
        # STEP 1: remove isotopes
        inchi1 = InChiParser.removeIsotopicLayers(inchi1)
        inchi2 = InChiParser.removeIsotopicLayers(inchi2)

        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        # STEP 2: remove salts
        mol1 = InChi.main_fragment(mol1)
        mol2 = InChi.main_fragment(mol2)

        # we convert back to InChI so we can inspect p/q layers
        inchi1 = Chem.MolToInchi(mol1)
        inchi2 = Chem.MolToInchi(mol2)

        inchi1 = inchi1.strip()
        inchi2 = inchi2.strip()

        # STEP 3: charge layers
        p1_plus, p1_minus, q1_plus, q1_minus = InChi.get_charge_info(inchi1)
        p2_plus, p2_minus, q2_plus, q2_minus = InChi.get_charge_info(inchi2)

        # CASE 1: p+N → remove p layer
        if (p1_plus or p2_plus) and not (p1_minus or p2_minus or q1_minus or q2_minus):
            inchi1 = InChi.remove_only_p_layer(inchi1)
            inchi2 = InChi.remove_only_p_layer(inchi2)

            return inchi1 == inchi2

        # CASE 2: q-N or p-N → neutralize molecules
        if p1_minus or p2_minus or q1_minus or q2_minus:

            mol1 = InChi.mol_from_inchi(inchi1)
            mol2 = InChi.mol_from_inchi(inchi2)

            if mol1 is None or mol2 is None:
                return False

            mol1 = InChi.neutralize_molecule(mol1)
            mol2 = InChi.neutralize_molecule(mol2)

            sig1 = Chem.MolToSmiles(mol1, canonical=True, isomericSmiles=False)
            sig2 = Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=False)

            return sig1 == sig2

        # CASE 3: q+N → we leave it
        return inchi1 == inchi2

 
    def areEqualNoStereo(inchi1: str, inchi2: str) -> bool:
        inchi1_no_stereo = InChiParser.removeStereoLayers(inchi1)
        inchi2_no_stereo = InChiParser.removeStereoLayers(inchi2)
        return inchi1_no_stereo == inchi2_no_stereo
 
    #stereochemical layer - sublayer
    def areEqualNoPositionDoubleBond(inchi1: str, inchi2: str) -> bool:
        if inchi1 == inchi2:
            return True

        # STEP 1: remove isotopes
        inchi1 = InChiParser.removeIsotopicLayers(inchi1)
        inchi2 = InChiParser.removeIsotopicLayers(inchi2)

        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        # STEP 2: remove salts
        mol1 = InChi.main_fragment(mol1)
        mol2 = InChi.main_fragment(mol2)

        # STEP 3: neutralize charges
        mol1 = InChi.neutralize_molecule(mol1)
        mol2 = InChi.neutralize_molecule(mol2)

        # STEP 4: lipid check
        if not (LipidAnalysis.is_lipid(inchi1, mol1) and
                LipidAnalysis.is_lipid(inchi2, mol2)):

            # fallback: simple structural comparison
            sig1 = Chem.MolToSmiles(mol1, canonical=True, isomericSmiles=False)
            sig2 = Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=False)
            return sig1 == sig2

        # STEP 5: remove cis/trans
        LipidAnalysis.remove_cis_trans(mol1)
        LipidAnalysis.remove_cis_trans(mol2)

        # STEP 6: extract tails
        tails1 = LipidAnalysis.extract_detailed_tails(mol1)
        tails2 = LipidAnalysis.extract_detailed_tails(mol2)

        # LEVEL A
        LEVELA = tails1 == tails2

        # LEVEL B
        LEVELB = sorted(tails1) == sorted(tails2)

        # LEVEL C
        sig1_C = sorted([(t["C"], t["DB"], t["O"]) for t in tails1])
        sig2_C = sorted([(t["C"], t["DB"], t["O"]) for t in tails2])
        LEVELC = sig1_C == sig2_C

        # LEVEL D
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

        return LEVELA or LEVELB or LEVELC or LEVELD
    
   
    def run_inchitrust(mol, inchitrust_path):
        try:
            molblock = Chem.MolToMolBlock(mol)
            process = subprocess.run(
                [inchitrust_path],
                input=molblock,
                text=True,
                capture_output=True,
                check=True
            )

            output = process.stdout.strip()
            return output

        except Exception as e:
            print(f"InChI Trust execution failed: {e}")
            return None
        
    def areEqualTautomers(inchi1: str, inchi2: str, inchitrust_path: str | None = None) -> bool:
        # if path not provided, read from environment variable
        if inchitrust_path is None:
            inchitrust_path = os.getenv("INCHITRUST_PATH")

        if not inchitrust_path:
            raise ValueError(
                "InChI Trust path not provided. Set INCHITRUST_PATH environment variable."
            )

        # STEP 1: remove isotope layers
        inchi1 = InChiParser.removeIsotopicLayers(inchi1)
        inchi2 = InChiParser.removeIsotopicLayers(inchi2)

        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        # STEP 2: remove salts
        mol1 = InChi.main_fragment(mol1)
        mol2 = InChi.main_fragment(mol2)

        # STEP 3: neutralize charges
        mol1 = InChi.neutralize_molecule(mol1)
        mol2 = InChi.neutralize_molecule(mol2)

        # STEP 4: generate canonical tautomer
        taut1 = InChi.run_inchitrust(mol1, inchitrust_path)
        taut2 = InChi.run_inchitrust(mol2, inchitrust_path)

        if taut1 is None or taut2 is None:
            return False

        return taut1 == taut2
    
    @staticmethod
    def get_substituent_signatures(mol):
        mol = InChi.main_fragment(mol)
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)

        # detect molecules with no Murcko scaffold (e.g. lipids)
        use_scaffold = scaffold.GetNumAtoms() < mol.GetNumAtoms()

        if use_scaffold:
            scaffold_smiles = Chem.MolToSmiles(
                scaffold,
                canonical=True,
                isomericSmiles=False
            )
            scaffold_atoms = {a.GetIdx() for a in scaffold.GetAtoms()}
        else:
            # no scaffold → treat whole molecule as substituent space
            scaffold_smiles = "NO_SCAFFOLD"
            scaffold_atoms = set()

        visited = set()
        substituents = []

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx in scaffold_atoms or idx in visited:
                continue

            stack = [idx]
            frag_atoms = set()

            while stack:
                current = stack.pop()
                if current in visited:
                    continue

                visited.add(current)
                frag_atoms.add(current)
                atom_obj = mol.GetAtomWithIdx(current)

                for nbr in atom_obj.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in scaffold_atoms and nbr_idx not in visited:
                        stack.append(nbr_idx)

            if frag_atoms:
                smiles = Chem.MolFragmentToSmiles(
                    mol,
                    atomsToUse=list(frag_atoms),
                    canonical=True,
                    isomericSmiles=False
                )
                substituents.append(smiles)

        return scaffold_smiles, Counter(substituents)
    
    @staticmethod
    def substituent_position_independent_signature(inchi: str):
        mol = InChi.mol_from_inchi(inchi)

        if mol is None:
            return None

        scaffold, subs = InChi.get_substituent_signatures(mol)
        return (scaffold, subs)
    
    
    @staticmethod
    def areEqualSubstituentIndependent(inchi1: str, inchi2: str) -> bool:
        # STEP 1: remove isotopes
        inchi1 = InChiParser.removeIsotopicLayers(inchi1)
        inchi2 = InChiParser.removeIsotopicLayers(inchi2)

        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        #STEP 2: remove salts
        mol1 = InChi.main_fragment(mol1)
        mol2 = InChi.main_fragment(mol2)

        #Step 3: charges
        mol1 = InChi.neutralize_molecule(mol1)
        mol2 = InChi.neutralize_molecule(mol2)

        #STEP 4: cis/trans 
        LipidAnalysis.remove_cis_trans(mol1)
        LipidAnalysis.remove_cis_trans(mol2)

        sig1 = LipidAnalysis.lipid_chain_signatures(mol1)
        sig2 = LipidAnalysis.lipid_chain_signatures(mol2)

        #case 1: molecule is a lipid
        if (LipidAnalysis.is_lipid(inchi1, mol1) and
            LipidAnalysis.is_lipid(inchi2, mol2)):
            return sig1 == sig2
        
        #case 2 - Murcko scaffold substituent comparison
        sig1 = InChi.substituent_position_independent_signature(inchi1)
        sig2 = InChi.substituent_position_independent_signature(inchi2)

        return sig1 == sig2


    @staticmethod
    def get_ids(inchi1: str, inchi2: str, config: dict) -> dict:
        criteria = config.get("identity_criteria", {})

        results = {}

        # COMPLETE IDENTITY
        if criteria["complete_identity"]["enabled"]:
            results[InchiLayers.COMPLETE_IDENTITY] = (
                InChi.isCompleteIdentity(inchi1, inchi2)
            )

        # ISOTOPES
        if criteria["isotope_independence"]["isotope_independent_identity"]:
            results[InchiLayers.ISOTOPIC_INDEPENDENCE] = (
                InChi.areEqualNoIsotopes(inchi1, inchi2)
            )

        # SALTS
        if criteria["salt_independence"]["desalted_identity"]:
            results[InchiLayers.SALTS_INDEPENDENCE] = (
                InChi.areEqualDisolvedSalts(inchi1, inchi2)
            )

        # CHARGES
        if criteria["charge_independence"]["charge_independent_identity"]:
            results[InchiLayers.CHARGES_INDEPENDENCE] = (
                InChi.areEqualNoCharges(inchi1, inchi2)
            )

        #ISOMER
        # DOUBLE BOND POSITION
        if criteria["isomer_independence"]["double_bond_position_independent_identity"]:
            results[InchiLayers.DOUBLE_BONDS_INDEPENDENCE] = (
                InChi.areEqualNoPositionDoubleBond(inchi1, inchi2)
            )

        # CIS/TRANS
        #TODO: if true then we apply it, if not, we dont
        if criteria["isomer_independence"]["cis_trans_independent_identity"]:
            results[InchiLayers.STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE] = (
                InChi.areEqualNoStereo(inchi1, inchi2)
            )

        # TAUTOMERS
        tautomer_cfg = criteria.get("tautomer_independence", {})

        if tautomer_cfg.get("tautomer_independent_identity", False):
            inchitrust_path = tautomer_cfg.get("inchitrust_path")
            results[InchiLayers.TAUTOMER_INDEPENDENCE] = (
                InChi.areEqualTautomers(inchi1, inchi2, inchitrust_path)
            )

        # SUBSTITUENT POSITION INDEPENDENCE
        if tautomer_cfg.get("substituent_position_independent_identity", False):
            results[InchiLayers.SUBSTITUENT_POSITION_INDEPENDENCE] = (
                InChi.areEqualSubstituentIndependent(inchi1, inchi2)
            )


        return results