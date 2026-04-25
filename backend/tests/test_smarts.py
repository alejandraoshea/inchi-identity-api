from rdkit import Chem
from inchi.lipid_analysis import LipidAnalysis
from inchi.lipid_structure_detector import LipidHeadValidator

if __name__ == "__main__":
    print("LIPID STRUCTURE VALIDATOR - Usage Examples")
    
    # Example 1: 
    print("\n1. glycosylglycerol:")
    
    core = "O[C@H]1[C@H](OC[C@]([H])(O)CO)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O"
    good = "O[C@H]1[C@H](OC[C@]([H])(O)COC(CCCCCC)=O)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O"
    wrong = "O[C@H]1[C@H](OC[C@]([H])(O)CO)O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)[C@H](OC(CCCCCC)=O)[C@H]2O)[C@H](O)[C@@H]1O"
    
    validator = LipidHeadValidator()
    
    mol_core = Chem.MolFromSmiles(core)
    mol_good = Chem.MolFromSmiles(good)
    mol_wrong = Chem.MolFromSmiles(wrong)
    
    print(f"CORE (no FA):  Valid = {validator.matches_any_valid_head(mol_core)}")
    print(f"GOOD (correct): Valid = {validator.matches_any_valid_head(mol_good)}")
    print(f"WRONG (incorrect): Valid = {validator.matches_any_valid_head(mol_wrong)}")
    
    print("\nExpected: CORE=False, GOOD=True, WRONG=False")
    
    # Example 2: Detailed validation
    print("\n2. Detailed validation of GOOD molecule:")
    result = validator.validate_structure(mol_good, verbose=True)
    
    # Example 3: Multiple patterns
    print("\n3. Testing different lipid classes:")
    
    test_cases = {
        "1,2-DG": "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC",
        "1,3-DG": "CCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCC",
        "TG": "CCCCCCCC(=O)OCC(COC(=O)CCCCCCCC)OC(=O)CCCCCCCC",
    }
    
    for name, smiles in test_cases.items():
        mol = Chem.MolFromSmiles(smiles)
        lipid_class = validator.identify_lipid_class(mol)
        print(f"{name:10s} → {lipid_class or 'Not recognized'}")
    
