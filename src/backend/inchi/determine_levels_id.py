from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdmolops
from collections import Counter
from rdkit.Chem import inchi, MolToSmiles
from rdkit.Chem.SaltRemover import SaltRemover
from backend.inchi.inchi_parser import InChIParser
from backend.inchi.inchi_layers_enum import InchiLayers
from backend.lipid.lipid_analysis import LipidAnalysis
from backend.lipid.lipid_analysis import LipidHeadValidator
from backend.lipid.lipid_tail_extraction import TailExtractor
import subprocess, os
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import Counter
from rdkit.Chem.MolStandardize import rdMolStandardize

class InChI:
    def isCompleteIdentity(inchi1: str, inchi2: str) -> bool:
        inchi1 = InChI.normalize_input(inchi1)
        inchi2 = InChI.normalize_input(inchi2)
        return inchi1 == inchi2

    def mol_from_inchi(inchi: str):
        try:
            mol = Chem.inchi.MolFromInchi(inchi.strip())
            if mol is not None:
                return mol
            
            mol = Chem.inchi.MolFromInchi(inchi.strip(), sanitize=False)
            if mol is None:
                raise ValueError("InChI conversion failed: Internal RDKit/InChI parser rejection")
            return mol
        except Exception as e:
            print(f"RDKit conversion failed for InChI: {inchi}\nError: {e}")
            return None
    
    @staticmethod
    def normalize_input(structure: str) -> str:
        if not structure:
            return structure
        
        structure = structure.strip()
        
        if structure.startswith("InChI="):
            return structure
        
        try:
            mol = Chem.MolFromSmiles(structure, sanitize=True)
            
            if mol is None:
                return structure
            
            Chem.SanitizeMol(mol)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            
            inchi_result = Chem.MolToInchi(mol)
            if inchi_result:
                return inchi_result.strip()
            else:
                print(f"InChI generation failed for valid SMILES: {structure}")
                return structure
                
        except Exception as e:
            print(f"\Error during SMILES→InChI conversion:")
            print(f"    Input: {structure}")
            print(f"    Error: {e}")
            return structure


    def areEqualNoIsotopes(inchi1: str, inchi2: str) -> bool:
        inchi1 = InChI.normalize_input(inchi1)
        inchi2 = InChI.normalize_input(inchi2)
    
        inchi1_isotopes = InChIParser.removeIsotopicLayers(inchi1)
        inchi2_isotopes = InChIParser.removeIsotopicLayers(inchi2)
        return inchi1_isotopes == inchi2_isotopes
    
    def main_fragment(mol):
        try:
            frags = Chem.GetMolFrags(mol, asMols=True)
            main = max(frags, key=lambda m: m.GetNumAtoms())
            Chem.SanitizeMol(main)
            Chem.AssignStereochemistry(main, cleanIt=True, force=True)
            return main
        except Exception:
            return mol

    def areEqualDisolvedSalts(inchi1: str, inchi2: str) -> bool:
        inchi1 = InChI.normalize_input(inchi1)
        inchi2 = InChI.normalize_input(inchi2)

        inchi1 = InChIParser.removeIsotopicLayers(inchi1)
        inchi2 = InChIParser.removeIsotopicLayers(inchi2)

        mol1 = InChI.mol_from_inchi(inchi1)
        mol2 = InChI.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        remover = SaltRemover()

        mol1_clean = remover.StripMol(mol1, dontRemoveEverything=True)
        mol2_clean = remover.StripMol(mol2, dontRemoveEverything=True)

        mol1_main = InChI.main_fragment(mol1_clean)
        mol2_main = InChI.main_fragment(mol2_clean)

        try:
            inchi1_final = inchi.MolToInchi(mol1_main)
            inchi2_final = inchi.MolToInchi(mol2_main)
        except Exception:
            return False

        return inchi1_final == inchi2_final
    
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
            if p: 
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
        inchi1 = InChI.normalize_input(inchi1)
        inchi2 = InChI.normalize_input(inchi2)

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

        # CASE 1: p+N - remove p layer
        if (p1_plus or p2_plus) and not (p1_minus or p2_minus or q1_minus or q2_minus):
            inchi1 = InChI.remove_only_p_layer(inchi1)
            inchi2 = InChI.remove_only_p_layer(inchi2)

            return inchi1 == inchi2

        # CASE 2: q-N or p-N - neutralize molecules
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

        # CASE 3: q+N - we leave it
        return inchi1 == inchi2
 
    def areEqualNoStereo(inchi1: str, inchi2: str) -> bool:
        if inchi1 == inchi2:
            return True

        inchi1 = InChI.normalize_input(inchi1)
        inchi2 = InChI.normalize_input(inchi2)

        inchi1 = InChIParser.removeIsotopicLayers(inchi1)
        inchi2 = InChIParser.removeIsotopicLayers(inchi2)

        mol1 = InChI.mol_from_inchi(inchi1)
        mol2 = InChI.mol_from_inchi(inchi2)

        if mol1 is None or mol2 is None:
            return False

        mol1 = InChI.main_fragment(mol1)
        mol2 = InChI.main_fragment(mol2)

        mol1 = InChI.neutralize_molecule(mol1)
        mol2 = InChI.neutralize_molecule(mol2)

        mol1 = InChI.remove_cis_trans(mol1)
        mol2 = InChI.remove_cis_trans(mol2)

        sig1 = Chem.MolToSmiles(mol1, canonical=True, isomericSmiles=False)
        sig2 = Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=False)
        return sig1 == sig2

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
    
    #stereochemical layer - sublayer
    def areEqualNoPositionDoubleBond(inchi1: str, inchi2: str) -> bool:
        if inchi1 == inchi2:
            return True
        
        inchi1 = InChI.normalize_input(inchi1)
        inchi2 = InChI.normalize_input(inchi2)

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
        mol1 = InChI.remove_cis_trans(mol1)
        mol2 = InChI.remove_cis_trans(mol2)

        sig1 = MolToSmiles(mol1, canonical=True, isomericSmiles=False)
        sig2 = MolToSmiles(mol2, canonical=True, isomericSmiles=False)

        if sig1 == sig2:
            return True

        tails1 = TailExtractor.extract_tails(mol1)
        tails2 = TailExtractor.extract_tails(mol2)

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
            import tempfile, os
            molblock = Chem.MolToMolBlock(mol)
            with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False, encoding='utf-8') as f:
                f.write(molblock)
                tmp_path = f.name
            
            result = subprocess.run(
                [inchitrust_path, tmp_path, "-", "/AuxNone", "/NoLabels", "/NoWarnings"],
                capture_output=True,
                text=True,
                encoding="utf-8"
            )
            os.unlink(tmp_path)
            
            output = "\n".join(
                line for line in result.stdout.splitlines()
                if line.startswith("InChI=")
            )
            return output if output else None
            
        except Exception as e:
            print(f"InChI Trust execution failed: {e}")
            return None
        
    def areEqualTautomers(inchi1: str, inchi2: str, inchitrust_path=None) -> bool:
        if inchitrust_path is None:
            inchitrust_path = os.environ.get("INCHITRUST_PATH")

        inchi1 = InChI.normalize_input(inchi1)
        inchi2 = InChI.normalize_input(inchi2)
        inchi1 = InChIParser.removeIsotopicLayers(inchi1)
        inchi2 = InChIParser.removeIsotopicLayers(inchi2)

        mol1 = InChI.mol_from_inchi(inchi1)
        mol2 = InChI.mol_from_inchi(inchi2)
        if mol1 is None or mol2 is None:
            return False

        mol1 = InChI.main_fragment(mol1)
        mol2 = InChI.main_fragment(mol2)
        mol1 = InChI.neutralize_molecule(mol1)
        mol2 = InChI.neutralize_molecule(mol2)

        if inchitrust_path:
            canon_inchi1 = InChI.run_inchitrust(mol1, inchitrust_path)
            canon_inchi2 = InChI.run_inchitrust(mol2, inchitrust_path)
            if canon_inchi1 is not None and canon_inchi2 is not None:
                return canon_inchi1 == canon_inchi2
        
        # fallback to RDKit
        enumerator = rdMolStandardize.TautomerEnumerator()
        canon1 = enumerator.Canonicalize(mol1)
        canon2 = enumerator.Canonicalize(mol2)
        sig1 = Chem.MolToSmiles(canon1, canonical=True, isomericSmiles=False)
        sig2 = Chem.MolToSmiles(canon2, canonical=True, isomericSmiles=False)
        return sig1 == sig2


    @staticmethod
    def areEqualSubstituentIndependent(inchi1: str, inchi2: str) -> bool:
        inchi1 = InChI.normalize_input(inchi1)
        inchi2 = InChI.normalize_input(inchi2)

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
        mol1 = InChI.remove_cis_trans(mol1)
        mol2 = InChI.remove_cis_trans(mol2)

        # STEP 5: canonicalize tautomers
        tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
        mol1 = tautomer_enumerator.Canonicalize(mol1)
        mol2 = tautomer_enumerator.Canonicalize(mol2)

        # STEP 6: lipid detection
        is_lipid1 = LipidAnalysis.is_lipid(inchi1, mol1, use_classyfire=False)
        is_lipid2 = LipidAnalysis.is_lipid(inchi2, mol2, use_classyfire=False)

        if is_lipid1 != is_lipid2:
            return False

        # CASE 1: both are lipids
        if is_lipid1 and is_lipid2:            
            tails1 = TailExtractor.extract_tails(mol1)
            tails2 = TailExtractor.extract_tails(mol2)

            if not tails1 or not tails2:
                return LipidAnalysis.atom_count(mol1) == LipidAnalysis.atom_count(mol2)

            sigs1 = sorted([LipidAnalysis.tail_sig_levelC(t) for t in tails1])
            sigs2 = sorted([LipidAnalysis.tail_sig_levelC(t) for t in tails2])

            if sigs1 == sigs2:
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

            # Atom count fallback
            if len(tails1) == len(tails2):
                if LipidAnalysis.atom_count(mol1) == LipidAnalysis.atom_count(mol2):
                    return True

            return False

        # CASE 2: both are not lipids
        sig1 = Chem.MolToSmiles(mol1, canonical=True, isomericSmiles=False)
        sig2 = Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=False)

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

        #STEREOCHEMICAL
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

        return results