from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.SaltRemover import SaltRemover
from inchi.inchi_parser import InChiParser
from inchi.inchi_layers_enum import InchiLayers

class InChi:
    def isCompleteIdentity(inchi1: str, inchi2: str) -> bool:
        return (
            InChi.main_layer(inchi1, inchi2) and
            InChi.charge_layer(inchi1, inchi2) and
            InChi.stereochemical_layer(inchi1, inchi2) and
            InChi.isotopic_layer(inchi1, inchi2) and
            InChi.fixed_H_layer(inchi1, inchi2) and
            InChi.reconnected_layer(inchi1, inchi2)
        )

    #main layer (sublayers: atom connections and hydrogen atoms)
    def main_layer(inchi1: str, inchi2: str) -> bool:
        return(
            InChiParser.getMainLayer(inchi1) == InChiParser.getMainLayer(inchi2) and
            InChiParser.getAtomConnectionsSublayer(inchi1) == InChiParser.getAtomConnectionsSublayer(inchi2) and
            InChiParser.getHydrogenAtomsSublayer(inchi1) == InChiParser.getHydrogenAtomsSublayer(inchi2)
        )
    
    #charge layer (sublayers: charge and proton)
    def charge_layer(inchi1: str, inchi2: str) -> bool:
        return (
            InChiParser.getChargeSublayer(inchi1) == InChiParser.getChargeSublayer(inchi2) and
            InChiParser.getProtonSublayer(inchi1) == InChiParser.getProtonSublayer(inchi2)
        )
    
    def compare_charge_effects(original_inchi, neutral_inchi):
        return {
            "orig_no_charge": InChiParser.removeChargeLayers(original_inchi),
            "neutral_no_charge": InChiParser.removeChargeLayers(neutral_inchi),
            "charge_independent_equal":
                InChiParser.removeChargeLayers(original_inchi)
                == InChiParser.removeChargeLayers(neutral_inchi)
        }
    
    #stereochemical layer (sublayers: double bonds, tetrahedrals, type)
    def stereochemical_layer(inchi1: str, inchi2: str) -> bool:
        return (
            InChiParser.getDoubleBondsSublayer(inchi1) == InChiParser.getDoubleBondsSublayer(inchi2) and
            InChiParser.getTetrahedralStereoSublayer(inchi1) == InChiParser.getTetrahedralStereoSublayer(inchi2) and
            InChiParser.getTypeStereoInfoSublayer(inchi1) == InChiParser.getTypeStereoInfoSublayer(inchi2)
        )
    
    def isotopic_layer(inchi1: str, inchi2: str) -> bool:
        return (
            InChiParser.getIsotopicLayer(inchi1) == InChiParser.getIsotopicLayer(inchi2) and
            InChiParser.getIsotopicHydrogenSublayer(inchi1) == InChiParser.getIsotopicHydrogenSublayer(inchi2) and
            InChiParser.getIsotopicStereoSublayer(inchi1) == InChiParser.getIsotopicStereoSublayer(inchi2)
        )

    def fixed_H_layer(inchi1: str, inchi2: str) -> bool:
        return InChiParser.getFixedHLayer(inchi1) == InChiParser.getFixedHLayer(inchi2)

    def reconnected_layer(inchi1: str, inchi2: str) -> bool:
        return InChiParser.getReconnectedLayer(inchi1) == InChiParser.getReconnectedLayer(inchi2)
    # never included in standard InChI

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
        remover = SaltRemover()
        try:
            mol_clean = remover.StripMol(mol, dontRemoveEverything=True)
            Chem.AssignStereochemistry(mol_clean, cleanIt=True, force=True)
            return mol_clean
        except Exception:
            return mol
    
    def areEqualNoIsotopes(inchi1: str, inchi2: str) -> bool:
        inchi1_isotopes = InChiParser.removeIsotopicLayers(inchi1)
        inchi2_isotopes = InChiParser.removeIsotopicLayers(inchi2)
        return inchi1_isotopes == inchi2_isotopes
    
    #TODO: ADD HIERARCHY
    def areEqualDisolvedSalts(inchi1: str, inchi2: str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)
        if not mol1 or not mol2:
            return False    
        main1 = InChi.main_fragment(mol1)
        main2 = InChi.main_fragment(mol2)      
        return Chem.MolToInchi(main1) == Chem.MolToInchi(main2)
    
    #helper method: detect negative charge to neutralize the mol
    def has_negative_charge(inchi: str) -> bool:
        for part in inchi.split("/"):
            if part.startswith("q-") or part.startswith("p-"):
                return True
        return False
    
    def areEqualNoCharges(inchi1: str, inchi2: str) -> bool:
        inchi1 = inchi1.strip()
        inchi2 = inchi2.strip()

        p1_plus, p1_minus, q1_plus, q1_minus = InChi.get_charge_info(inchi1)
        p2_plus, p2_minus, q2_plus, q2_minus = InChi.get_charge_info(inchi2)

        # CASE 1: p+N → remove p layer
        if (p1_plus or p2_plus) and not (p1_minus or p2_minus or q1_minus or q2_minus):
            inchi1 = InChi.remove_only_p_layer(inchi1)
            inchi2 = InChi.remove_only_p_layer(inchi2)
            return inchi1 == inchi2

        # CASE 2: q-N or p-N → neutralize
        if p1_minus or p2_minus or q1_minus or q2_minus:
            mol1 = InChi.mol_from_inchi(inchi1)
            mol2 = InChi.mol_from_inchi(inchi2)

            if mol1 is None or mol2 is None:
                return False

            mol1 = InChiParser.neutralize_molecule(mol1)
            mol2 = InChiParser.neutralize_molecule(mol2)

            sig1 = Chem.MolToSmiles(mol1, canonical=True, isomericSmiles=False)
            sig2 = Chem.MolToSmiles(mol2, canonical=True, isomericSmiles=False)

            return sig1 == sig2

        # CASE 3: q+N → leave as is
        return inchi1 == inchi2

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

    def areEqualNoStereo(inchi1: str, inchi2: str) -> bool:
        inchi1_no_stereo = InChiParser.removeStereoLayers(inchi1)
        inchi2_no_stereo = InChiParser.removeStereoLayers(inchi2)
        return inchi1_no_stereo == inchi2_no_stereo

    #stereochemical layer - sublayer
    def areEqualNoPositionDoubleBond(inchi1: str, inchi2:str) -> bool:
        inchi1_no_double_bonds = InChiParser.removeDoubleBondsSublayer(inchi1)
        inchi2_no_double_bonds = InChiParser.removeDoubleBondsSublayer(inchi2)
        return inchi1_no_double_bonds == inchi2_no_double_bonds
    
    def areEqualTautomers(inchi1: str, inchi2: str) -> bool:
        sig1 = InChiParser.getTautomerLayer(inchi1)
        sig2 = InChiParser.getTautomerLayer(inchi2)
        if sig1 is None or sig2 is None:
            return False
        return sig1 == sig2

    @staticmethod
    def get_ids(inchi1: str, inchi2: str) -> dict:
        #dict<InchiLayers, bool>: for every identity rule, returns whether true/false for each layer
            results = {}
            results[InchiLayers.COMPLETE_IDENTITY] = (
                InChi.isCompleteIdentity(inchi1, inchi2)
            )
            results[InchiLayers.INDEPENDENT_SALTS] = (
                InChi.areEqualDisolvedSalts(inchi1, inchi2) #check when removing salts: stereo
            )
            results[InchiLayers.INDEPENDENT_CHARGES] = (
                InChi.areEqualNoCharges(inchi1, inchi2)
            )
            results[InchiLayers.INDEPENDENT_DOUBLE_BONDS] = (
                InChi.areEqualNoPositionDoubleBond(inchi1, inchi2)
            )
            results[InchiLayers.TAUTOMERIC] = (
                InChi.areEqualTautomers(inchi1, inchi2)
            )
            results[InchiLayers.STEREOCHEMICAL] = (
                InChi.areEqualNoStereo(inchi1, inchi2)
            )
            results[InchiLayers.ISOTOPIC] = (
                InChi.areEqualNoIsotopes(inchi1, inchi2)
            )
            return results