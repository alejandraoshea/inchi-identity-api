from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdmolops
from collections import Counter
from rdkit.Chem import inchi, MolToSmiles
from rdkit.Chem.SaltRemover import SaltRemover
from backend.inchi.inchi_parser import InChIParser
from backend.inchi.inchi_layers_enum import InchiLayers
from backend.inchi.lipid_analysis import LipidAnalysis
import subprocess, os
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import Counter
from rdkit.Chem.MolStandardize import rdMolStandardize

class InChI:
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
        inchi1_isotopes = InChIParser.removeIsotopicLayers(inchi1)
        inchi2_isotopes = InChIParser.removeIsotopicLayers(inchi2)
        return inchi1_isotopes == inchi2_isotopes
    
    def areEqualDisolvedSalts(inchi1: str, inchi2: str) -> bool:
        mol1 = InChI.mol_from_inchi(inchi1)
        mol2 = InChI.mol_from_inchi(inchi2)

        if not mol1 or not mol2:
            return False    

        main1 = InChI.main_fragment(mol1)
        main2 = InChI.main_fragment(mol2)

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
        inchi1 = InChIParser.removeIsotopicLayers(inchi1)
        inchi2 = InChIParser.removeIsotopicLayers(inchi2)

        mol1 = InChI.mol_from_inchi(inchi1)
        mol2 = InChI.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        # STEP 2: remove salts
        mol1 = InChI.main_fragment(mol1)
        mol2 = InChI.main_fragment(mol2)

        # we convert back to InChI so we can inspect p/q layers
        inchi1 = Chem.MolToInchi(mol1)
        inchi2 = Chem.MolToInchi(mol2)

        inchi1 = inchi1.strip()
        inchi2 = inchi2.strip()

        # STEP 3: charge layers
        p1_plus, p1_minus, q1_plus, q1_minus = InChI.get_charge_info(inchi1)
        p2_plus, p2_minus, q2_plus, q2_minus = InChI.get_charge_info(inchi2)

        # CASE 1: p+N → remove p layer
        if (p1_plus or p2_plus) and not (p1_minus or p2_minus or q1_minus or q2_minus):
            inchi1 = InChI.remove_only_p_layer(inchi1)
            inchi2 = InChI.remove_only_p_layer(inchi2)

            return inchi1 == inchi2

        # CASE 2: q-N or p-N → neutralize molecules
        if p1_minus or p2_minus or q1_minus or q2_minus:

            mol1 = InChI.mol_from_inchi(inchi1)
            mol2 = InChI.mol_from_inchi(inchi2)

            if mol1 is None or mol2 is None:
                return False

            mol1 = InChI.neutralize_molecule(mol1)
            mol2 = InChI.neutralize_molecule(mol2)

            sig1 = Chem.MolToSmiles(mol1, canonical=True, isomericSmiles=False)
            sig2 = Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=False)

            return sig1 == sig2

        # CASE 3: q+N → we leave it
        return inchi1 == inchi2

 
    def areEqualNoStereo(inchi1: str, inchi2: str) -> bool:
        inchi1_no_stereo = InChIParser.removeStereoLayers(inchi1)
        inchi2_no_stereo = InChIParser.removeStereoLayers(inchi2)
        return inchi1_no_stereo == inchi2_no_stereo
 
    #stereochemical layer - sublayer
    def areEqualNoPositionDoubleBond(inchi1: str, inchi2: str) -> bool:
        if inchi1 == inchi2:
            return True

        # STEP 1: remove isotopes
        inchi1 = InChIParser.removeIsotopicLayers(inchi1)
        inchi2 = InChIParser.removeIsotopicLayers(inchi2)

        mol1 = InChI.mol_from_inchi(inchi1)
        mol2 = InChI.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        # STEP 2: remove salts
        mol1 = InChI.main_fragment(mol1)
        mol2 = InChI.main_fragment(mol2)

        # STEP 3: neutralize charges
        mol1 = InChI.neutralize_molecule(mol1)
        mol2 = InChI.neutralize_molecule(mol2)

        # STEP 4: lipid check
        if not (LipidAnalysis.is_lipid(inchi1, mol1) and
                LipidAnalysis.is_lipid(inchi2, mol2)):

            # fallback: simple structural comparison
            sig1 = Chem.MolToSmiles(mol1, canonical=True, isomericSmiles=False)
            sig2 = Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=False)
            return sig1 == sig2

        # LEVEL A: exact except cis/trans
        mol1 = LipidAnalysis.remove_cis_trans(mol1)
        mol2 = LipidAnalysis.remove_cis_trans(mol2)

        sig1 = MolToSmiles(mol1, canonical=True, isomericSmiles=False)
        sig2 = MolToSmiles(mol2, canonical=True, isomericSmiles=False)

        if sig1 == sig2:
            return True

        tails1 = LipidAnalysis.extract_tails(mol1)
        tails2 = LipidAnalysis.extract_tails(mol2)

        if not tails1 or not tails2:
            return False

        # LEVEL B: same chains, ignore position
        sig1 = sorted(LipidAnalysis.tail_sig_levelB(t) for t in tails1)
        sig2 = sorted(LipidAnalysis.tail_sig_levelB(t) for t in tails2)

        if sig1 == sig2:
            return True

        # LEVEL C: ignore DB/O positions
        sig1 = sorted(LipidAnalysis.tail_sig_levelC(t) for t in tails1)
        sig2 = sorted(LipidAnalysis.tail_sig_levelC(t) for t in tails2)

        if sig1 == sig2:
            return True

        # LEVEL D: global totals
        total1 = (
            sum(t["C"] for t in tails1),
            sum(t["DB"] for t in tails1),
        )

        total2 = (
            sum(t["C"] for t in tails2),
            sum(t["DB"] for t in tails2),
        )

        if total1 == total2:
            return True

        return False
    
   
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
        inchi1 = InChIParser.removeIsotopicLayers(inchi1)
        inchi2 = InChIParser.removeIsotopicLayers(inchi2)

        mol1 = InChI.mol_from_inchi(inchi1)
        mol2 = InChI.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        # STEP 2: remove salts
        mol1 = InChI.main_fragment(mol1)
        mol2 = InChI.main_fragment(mol2)

        # STEP 3: neutralize charges
        mol1 = InChI.neutralize_molecule(mol1)
        mol2 = InChI.neutralize_molecule(mol2)

        # STEP 4: generate canonical tautomer
        taut1 = InChI.run_inchitrust(mol1, inchitrust_path)
        taut2 = InChI.run_inchitrust(mol2, inchitrust_path)

        if taut1 is None or taut2 is None:
            return False

        return taut1 == taut2
    
    @staticmethod
    def get_substituent_signatures(mol):
        mol = InChI.main_fragment(mol)
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
        mol = InChI.mol_from_inchi(inchi)

        if mol is None:
            return None

        scaffold, subs = InChI.get_substituent_signatures(mol)
        return (scaffold, subs)
    
    
    @staticmethod
    def areEqualSubstituentIndependent(inchi1: str, inchi2: str) -> bool:
        # STEP 1: remove isotopes
        inchi1 = InChIParser.removeIsotopicLayers(inchi1)
        inchi2 = InChIParser.removeIsotopicLayers(inchi2)

        mol1 = InChI.mol_from_inchi(inchi1)
        mol2 = InChI.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        # STEP 2: remove salts
        mol1 = InChI.main_fragment(mol1)
        mol2 = InChI.main_fragment(mol2)

        # STEP 3: neutralize charges
        mol1 = InChI.neutralize_molecule(mol1)
        mol2 = InChI.neutralize_molecule(mol2)

        # STEP 4: remove cis/trans
        mol1 = LipidAnalysis.remove_cis_trans(mol1)
        mol2 = LipidAnalysis.remove_cis_trans(mol2)

        # STEP 5: canonicalize tautomers
        tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
        mol1 = tautomer_enumerator.Canonicalize(mol1)
        mol2 = tautomer_enumerator.Canonicalize(mol2)

        is_lipid1 = LipidAnalysis.is_lipid(inchi1, mol1, use_classyfire=False)
        is_lipid2 = LipidAnalysis.is_lipid(inchi2, mol2, use_classyfire=False)

        if is_lipid1 and is_lipid2:
            tails1 = LipidAnalysis.extract_tails(mol1)
            tails2 = LipidAnalysis.extract_tails(mol2)

            sigs1 = [LipidAnalysis.tail_sig_levelC(t) for t in tails1]
            sigs2 = [LipidAnalysis.tail_sig_levelC(t) for t in tails2]

            if Counter(sigs1) == Counter(sigs2):
                return True

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

            if total1 == total2:
                return True

            if LipidAnalysis.atom_count(mol1) == LipidAnalysis.atom_count(mol2):
                return True

            return False
        
        # CASE 2: not lipids - use Murcko scaffold substituent comparison
        sig1 = InChI.substituent_position_independent_signature(inchi1)
        sig2 = InChI.substituent_position_independent_signature(inchi2)

        return sig1 == sig2


    @staticmethod
    def get_ids(inchi1: str, inchi2: str, config: dict) -> dict:
        criteria = config.get("identity_criteria", {})

        results = {}

        # COMPLETE IDENTITY
        if criteria["complete_identity"]["enabled"]:
            results[InchiLayers.COMPLETE_IDENTITY] = (
                InChI.isCompleteIdentity(inchi1, inchi2)
            )

        # ISOTOPES
        if criteria["isotope_independence"]["isotope_independent_identity"]:
            results[InchiLayers.ISOTOPIC_INDEPENDENCE] = (
                InChI.areEqualNoIsotopes(inchi1, inchi2)
            )

        # SALTS
        if criteria["salt_independence"]["desalted_identity"]:
            results[InchiLayers.SALTS_INDEPENDENCE] = (
                InChI.areEqualDisolvedSalts(inchi1, inchi2)
            )

        # CHARGES
        if criteria["charge_independence"]["charge_independent_identity"]:
            results[InchiLayers.CHARGES_INDEPENDENCE] = (
                InChI.areEqualNoCharges(inchi1, inchi2)
            )

        #ISOMER
        # DOUBLE BOND POSITION
        if criteria["isomer_independence"]["double_bond_position_independent_identity"]:
            results[InchiLayers.DOUBLE_BONDS_INDEPENDENCE] = (
                InChI.areEqualNoPositionDoubleBond(inchi1, inchi2)
            )

        # CIS/TRANS
        if criteria["isomer_independence"]["cis_trans_independent_identity"]:
            results[InchiLayers.STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE] = (
                InChI.areEqualNoStereo(inchi1, inchi2)
            )

        # TAUTOMERS
        tautomer_cfg = criteria.get("tautomer_independence", {})

        if tautomer_cfg.get("tautomer_independent_identity", False):
            inchitrust_path = tautomer_cfg.get("inchitrust_path")
            results[InchiLayers.TAUTOMER_INDEPENDENCE] = (
                InChI.areEqualTautomers(inchi1, inchi2, inchitrust_path)
            )

        # SUBSTITUENT POSITION INDEPENDENCE
        if tautomer_cfg.get("substituent_position_independent_identity", False):
            results[InchiLayers.SUBSTITUENT_POSITION_INDEPENDENCE] = (
                InChI.areEqualSubstituentIndependent(inchi1, inchi2)
            )

        return results