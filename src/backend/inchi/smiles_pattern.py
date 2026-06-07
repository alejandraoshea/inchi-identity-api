from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import re


class SmilesCorrector:    
    @staticmethod
    def attempt_sanitize_strategies(smiles: str) -> tuple:
        # Strategy 1: Standard parsing (no special handling)
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=True)
            if mol is not None:
                return mol, "standard_parsing", None
        except Exception as e:
            pass
        
        # Strategy 2: Parse without sanitization, then sanitize
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                    return mol, "sanitize_after_parse", None
                except Exception:
                    pass
        except Exception as e:
            pass
        
        # Strategy 3: Try kekulizing (aromatic issues)
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                try:
                    Chem.Kekulize(mol, clearAromaticFlags=False)
                    Chem.SanitizeMol(mol)
                    return mol, "kekulize", None
                except Exception:
                    pass
        except Exception as e:
            pass
        
        # Strategy 4: Remove charges and retry (for over-charged molecules)
        try:
            smiles_no_charge = re.sub(r'\[([A-Za-z]+)([+-]\d*)\]', r'[\1]', smiles)
            if smiles_no_charge != smiles:
                mol = Chem.MolFromSmiles(smiles_no_charge, sanitize=True)
                if mol is not None:
                    return mol, "remove_charges", None
        except Exception as e:
            pass
        
        # Strategy 5: Fix hypervalent atoms then InChI roundtrip.
        # Handles sloppy SMILES like betaine CN(C)(C)CC(=O)O where N has 4 bonds
        # but no explicit charge — we infer the charge from excess bond order so
        # the full sanitization and InChI generation capture the right charge layers.
        try:
            _NEUTRAL_VALENCE = {6: 4, 7: 3, 8: 2, 16: 2, 15: 3}
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                rw = Chem.RWMol(mol)
                for atom in rw.GetAtoms():
                    an = atom.GetAtomicNum()
                    nv = _NEUTRAL_VALENCE.get(an)
                    if nv is not None and atom.GetFormalCharge() == 0:
                        bond_order = int(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()))
                        bond_order += atom.GetNumExplicitHs()
                        excess = bond_order - nv
                        if excess > 0:
                            atom.SetFormalCharge(excess)
                mol = rw.GetMol()
                Chem.SanitizeMol(mol)
                inchi = Chem.MolToInchi(mol)
                if inchi:
                    mol_from_inchi = Chem.MolFromInchi(inchi)
                    if mol_from_inchi is not None:
                        return mol_from_inchi, "hypervalent_fix_inchi_roundtrip", None
        except Exception as e:
            pass
        
        # Strategy 6: Cleanup with RDKit's standardizer
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize
            cleaner = rdMolStandardize.Cleanup(Chem.MolFromSmiles(smiles, sanitize=False))
            if cleaner is not None:
                return cleaner, "rdkit_cleanup", None
        except Exception as e:
            pass
        
        return None, None, "All sanitization strategies failed"
    
    @staticmethod
    def auto_correct(smiles: str, verbose: bool = False) -> dict:
        result = {
            "original": smiles,
            "corrected": smiles,
            "was_corrected": False,
            "parse_result": "unknown",
            "message": None,
            "strategy": None,
            "mol": None
        }
        
        if not smiles or not isinstance(smiles, str):
            result["parse_result"] = "error"
            result["message"] = "Empty or invalid input"
            return result
        
        smiles = smiles.strip()
        
        if smiles.startswith("InChI="):
            result["parse_result"] = "success"
            result["message"] = "Already InChI format (no correction needed)"
            return result
        
        mol, strategy, error = SmilesCorrector.attempt_sanitize_strategies(smiles)
        
        if mol is not None:
            try:
                corrected_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                result["mol"] = mol
                result["corrected"] = corrected_smiles
                result["strategy"] = strategy
                result["parse_result"] = "success"
                result["was_corrected"] = (corrected_smiles != smiles)
                
                if result["was_corrected"]:
                    result["message"] = f"Auto-corrected using {strategy}: '{smiles}' → '{corrected_smiles}'"
                else:
                    result["message"] = f"Valid SMILES (parsed with {strategy})"
                
                if verbose:
                    print(f"Success {result['message']}")
                
                return result
            
            except Exception as e:
                result["parse_result"] = "partial_success"
                result["message"] = f"Parsed but could not convert back to SMILES: {e}"
                result["mol"] = mol
                return result
        
        result["parse_result"] = "error"
        result["message"] = f"Could not parse SMILES with any strategy. Last error: {error}"
        
        if verbose:
            print(f"Error {result['message']}")
        
        return result
    
    @staticmethod
    def validate_chemistry(mol) -> dict:
        if mol is None:
            return {"valid": False, "issues": ["Molecule is None"]}
        
        issues = []
        
        try:
            bad_valence = Chem.DetectChemistryProblems(mol)
            if bad_valence:
                issues.append(f"Valence issues: {bad_valence}")
        except:
            pass
        
        if mol.GetNumAtoms() == 0:
            issues.append("Molecule has no atoms")
        
        if mol.GetNumAtoms() > 500:
            issues.append(f"Suspiciously large molecule ({mol.GetNumAtoms()} atoms)")
        
        return {
            "valid": len(issues) == 0,
            "issues": issues,
            "num_atoms": mol.GetNumAtoms() if mol else 0
        }

