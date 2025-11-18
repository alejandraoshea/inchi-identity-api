from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.SaltRemover import SaltRemover

class InChi:
    def mol_from_inchi(inchi: str):
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                raise ValueError(f"Invalid InChI: {inchi}")
            return mol
        except Exception:
            return None

    def _split(inchi: str):
        #split InChI into '/' components
        return inchi.split('/')

    def isCompleteIdentity(inchi1: str, inchi2: str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)
        if not mol1 or not mol2:
            return False
        return Chem.MolToInchi(mol1) == Chem.MolToInchi(mol2)

    def main_layer(inchi1: str, inchi2: str) -> bool:
        def extract(parts):
            return [p for p in parts if (not p[0].isalpha() and not p.startswith("InChI=")) or p.startswith(('c', 'h'))]
        p1 = extract(InChi._split(inchi1))
        p2 = extract(InChi._split(inchi2))
        return p1 == p2

    #Atom connections (prefix: "c"). The atoms in the chemical formula (except for hydrogens) are numbered in sequence; this sublayer describes which atoms are connected by bonds to which other ones. The type of those bonds is later specified in the stereochemical layer prefixed by "b".
    def charge_layer(inchi1: str, inchi2: str) -> bool:
        def extract(parts):
            return [p for p in parts if p.startswith(('q', 'p'))]
        p1 = extract(InChi._split(inchi1))
        p2 = extract(InChi._split(inchi2))
        return p1 == p2
    #charge sublayer (prefix: "q")
    #proton sublayer (prefix: "p" for "protons")

    
    #Hydrogen atoms (prefix: "h")
    def stereochemical_layer(inchi1: str, inchi2: str) -> bool:
        def extract(parts):
            return [p for p in parts if p.startswith(('b', 't', 'm', 's'))]
        p1 = extract(InChi._split(inchi1))
        p2 = extract(InChi._split(inchi2))
        return p1 == p2
    #double bonds and cumulenes (prefix: "b")
    #tetrahedral stereochemistry of atoms and allenes (prefixes: "t", "m")
    #type of stereochemistry information (prefix: "s")

    def isotopic_layer(inchi1: str, inchi2: str) -> bool:
        def extract(parts):
            return [p for p in parts if p.startswith('i')]
        p1 = extract(InChi._split(inchi1))
        p2 = extract(InChi._split(inchi2))
        return p1 == p2
    #Isotopic layer (prefix: "i"), may include sublayers:[13]
    #sublayer "h" for isotopic hydrogen
    #sublayers "b", "t", "m", "s" for isotopic stereochemistry


    def fixed_H_layer(inchi1: str, inchi2: str) -> bool:
        def extract(parts):
            return [p for p in parts if p.startswith('f')]
        p1 = extract(InChi._split(inchi1))
        p2 = extract(InChi._split(inchi2))
        return p1 == p2
    #Fixed-H layer (prefix: "f") for tautomeric hydrogens; contains some or all of the above types of layers except atom connections; 
    # may end with "o" sublayer; never included in standard InChI[13]

    def reconnected_layer(inchi1: str, inchi2: str) -> bool:
        def extract(parts):
            return [p for p in parts if p.startswith('r')]
        p1 = extract(InChi._split(inchi1))
        p2 = extract(InChi._split(inchi2))
        return p1 == p2
    #Reconnected layer (prefix: "r"); contains the whole InChI of a structure with reconnected metal atoms; 
    # never included in standard InChI

    def main_fragment(mol):
        remover = SaltRemover()
        try:
            mol_clean = remover.StripMol(mol, dontRemoveEverything=True) #strip the salt
            return mol_clean #contains the desalted molecule
        except Exception:
            return mol  

    def areEqualDisolvedSalts(inchi1: str, inchi2: str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)
        
        if not mol1 or not mol2:
            return False
        
        main1 = InChi.main_fragment(mol1)
        main2 = InChi.main_fragment(mol2)
        
        return Chem.MolToInchi(main1) == Chem.MolToInchi(main2)


    def areEqualNoPositionDoubleBond(inchi1: str, inchi2:str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)
        if not mol1 or not mol2:
            return False

        mcs = rdFMCS.FindMCS(
            [mol1, mol2],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder
        )
        return mcs.numAtoms == mol1.GetNumAtoms() == mol2.GetNumAtoms()


    def areEqualNoIsotopes(inchi1: str, inchi2: str) -> bool:
        mol1 = InChi.mol_from_inchi(inchi1)
        mol2 = InChi.mol_from_inchi(inchi2)
        
        if not mol1 or not mol2:
            return False
        
        mcs = rdFMCS.FindMCS(
            [mol1, mol2],
            atomCompare=rdFMCS.AtomCompare.CompareElements,  
            bondCompare=rdFMCS.BondCompare.CompareOrderExact
        )
        
        return mcs.numAtoms == mol1.GetNumAtoms() == mol2.GetNumAtoms()



#TODO:
#MAP<IDS,BOOLEAN>GETIDS(INCHI1,INCHI2)
#areEqualNoIsotopes(): RDKit puede sustituir
#areEqualNoSustituyentes(): Smiles? ??
#areEqualNoTautomeric(): ver Inchi Layers
#areEqualNoPositionDoubleBond(inchi1,inchi2): formula igual que el resto????