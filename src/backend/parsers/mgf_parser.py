import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, asdict

from backend.inchi.determine_levels_id import InChI
from backend.inchi.inchi_layers_enum import InchiLayers
from backend.inchi.inchi_parser import InChIParser
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
        

class MgfParser:

    def parse_mgf(file_path):
        entries = []
        current = {}

        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()

                if line == "BEGIN IONS":
                    current = {}

                elif line == "END IONS":
                    if current:
                        entries.append(current)

                elif "=" in line:
                    key, value = line.split("=", 1)
                    current[key.upper()] = value

        return entries

    def extract_inchis(entries):
        inchis = []

        for entry in entries:
            if "INCHI" in entry:
                inchis.append(entry["INCHI"])
            elif "SMILES" in entry:
                inchis.append(entry["SMILES"]) 

        return inchis
    
    def extract_structures(entries):
        result = []

        for entry in entries:
            structure = None

            if "INCHI" in entry:
                structure = entry["INCHI"]
            elif "SMILES" in entry:
                structure = entry["SMILES"]

            if structure:
                result.append({
                    "structure": structure,
                    "entry": entry
                })

        return result
    

@dataclass
class UnificationChange:
    original_structure: str
    intermediate_inchi: Optional[str]
    canonical_inchi: str
    structure_type: str


class SimpleMgfDeduplicator:
    
    def __init__(self, level: str = "COMPLETE_IDENTITY", config: dict = None):
        self.level = level
        self.config = config
        self.changes_log: List[UnificationChange] = []
        
        self.level_map = {
            "COMPLETE_IDENTITY": InchiLayers.COMPLETE_IDENTITY,
            "ISOTOPIC_INDEPENDENCE": InchiLayers.ISOTOPIC_INDEPENDENCE,
            "SALTS_INDEPENDENCE": InchiLayers.SALTS_INDEPENDENCE,
            "CHARGES_INDEPENDENCE": InchiLayers.CHARGES_INDEPENDENCE,
            "DOUBLE_BONDS_INDEPENDENCE": InchiLayers.DOUBLE_BONDS_INDEPENDENCE,
            "STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE": InchiLayers.STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE,
            "TAUTOMER_INDEPENDENCE": InchiLayers.TAUTOMER_INDEPENDENCE,
            "SUBSTITUENT_POSITION_INDEPENDENCE": InchiLayers.SUBSTITUENT_POSITION_INDEPENDENCE
        }
    
    def parse_mgf(self, file_path: str) -> List[Dict]:
        entries = []
        current = {}
        current_peaks = []

        with open(file_path, "r") as f:
            for line in f:
                line_stripped = line.strip()

                if line_stripped == "BEGIN IONS":
                    current = {}
                    current_peaks = []

                elif line_stripped == "END IONS":
                    if current:
                        if current_peaks:
                            current["_PEAKS"] = current_peaks
                        entries.append(current)

                elif "=" in line_stripped:
                    key, value = line_stripped.split("=", 1)
                    current[key.upper()] = value
                
                else:
                    if line_stripped and current:
                        current_peaks.append(line_stripped)

        return entries
    
    def extract_structure(self, entry: Dict) -> Tuple[Optional[str], Optional[str]]:
        if "INCHI" in entry:
            return entry["INCHI"], "INCHI"
        elif "SMILES" in entry:
            return entry["SMILES"], "SMILES"
        else:
            return None, None
    
    def extract_inchi(self, entry: Dict) -> Optional[str]:
        structure, _ = self.extract_structure(entry)
        return structure
    
    def get_canonical_form(self, inchi: str) -> str:
        if not inchi:
            return inchi
        
        inchi = inchi.strip()
        
        if not inchi.startswith("InChI="):
            try:
                mol = Chem.MolFromSmiles(inchi)
                if mol:
                    inchi = Chem.MolToInchi(mol)
                else:
                    return inchi
            except:
                return inchi
        
        if self.level == "COMPLETE_IDENTITY":
            return inchi
        
        if self.level == "ISOTOPIC_INDEPENDENCE":
            try:
                return InChIParser.removeIsotopicLayers(inchi)
            except Exception as e:
                return inchi
        
        try:    
            if self.level == "SALTS_INDEPENDENCE":
                inchi = InChIParser.removeIsotopicLayers(inchi)
                mol = InChI.mol_from_inchi(inchi)
                if mol is None:
                    return inchi
                
                remover = SaltRemover()
                mol_clean = remover.StripMol(mol, dontRemoveEverything=True)
                mol_main = InChI.main_fragment(mol_clean)
                
                return Chem.MolToInchi(mol_main)
            
            if self.level == "CHARGES_INDEPENDENCE":
                inchi = InChIParser.removeIsotopicLayers(inchi)
                mol = InChI.mol_from_inchi(inchi)
                if mol is None:
                    return inchi
                
                mol = InChI.main_fragment(mol)
                inchi_temp = Chem.MolToInchi(mol)
                
                p_plus, p_minus, q_plus, q_minus = InChI.get_charge_info(inchi_temp)
                
                if (p_plus) and not (p_minus or q_minus):
                    return InChI.remove_only_p_layer(inchi_temp)
                
                elif p_minus or q_minus:
                    mol = InChI.mol_from_inchi(inchi_temp)
                    if mol is None:
                        return inchi_temp
                    mol = InChI.neutralize_molecule(mol)
                    normalized_inchi = Chem.MolToInchi(mol)
                    return normalized_inchi
                
                else:
                    return inchi_temp

            
            if self.level == "STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE":
                inchi = InChIParser.removeIsotopicLayers(inchi)
                mol = InChI.mol_from_inchi(inchi)
                if mol is None:
                    return inchi
                
                mol = InChI.main_fragment(mol)
                mol = InChI.neutralize_molecule(mol)
                mol = InChI.remove_cis_trans(mol)
                
                return Chem.MolToInchi(mol)
            
            if self.level == "DOUBLE_BONDS_INDEPENDENCE":
                inchi = InChIParser.removeIsotopicLayers(inchi)
                mol = InChI.mol_from_inchi(inchi)
                if mol is None:
                    return inchi
                
                mol = InChI.main_fragment(mol)
                mol = InChI.neutralize_molecule(mol)
                mol = InChI.remove_cis_trans(mol)
                
                return Chem.MolToInchi(mol)
            
            if self.level == "TAUTOMER_INDEPENDENCE":
                inchi = InChIParser.removeIsotopicLayers(inchi)
                mol = InChI.mol_from_inchi(inchi)
                if mol is None:
                    return inchi
                
                mol = InChI.main_fragment(mol)
                mol = InChI.neutralize_molecule(mol)
                mol = InChI.remove_cis_trans(mol)
                tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
                mol = tautomer_enumerator.Canonicalize(mol)
                
                return Chem.MolToInchi(mol)
            
            if self.level == "SUBSTITUENT_POSITION_INDEPENDENCE":
                inchi = InChIParser.removeIsotopicLayers(inchi)
                mol = InChI.mol_from_inchi(inchi)
                if mol is None:
                    return inchi
                
                mol = InChI.main_fragment(mol)
                mol = InChI.neutralize_molecule(mol)
                mol = InChI.remove_cis_trans(mol)
                
                tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
                mol = tautomer_enumerator.Canonicalize(mol)
                
                return Chem.MolToInchi(mol)
            
        except Exception as e:
            print(f"Warning: Could not get canonical form at level {self.level}: {e}")
            return inchi
        
        return inchi
    
    def structures_match(self, struct1: str, struct2: str) -> bool:
        if not struct1 or not struct2:
            return False
        
        canonical1 = self.get_canonical_form(struct1)
        canonical2 = self.get_canonical_form(struct2)
        
        if self.level == "COMPLETE_IDENTITY":
            return canonical1 == canonical2
        
        if not self.config:
            return canonical1 == canonical2
        
        try:
            comparison = InChI.get_ids(canonical1, canonical2, self.config)
            level_enum = self.level_map.get(self.level, InchiLayers.COMPLETE_IDENTITY)
            return comparison.get(level_enum, False)
        except Exception as e:
            print(f"Error comparing structures at level {self.level}: {e}")
            return canonical1 == canonical2
    
    def unify_inchis_in_file(
        self, 
        entries: List[Dict], 
        source_file: str
    ) -> List[Dict]:
        canonical_map = {}
        
        for entry in entries:
            structure, field_name = self.extract_structure(entry)
            if not structure:
                continue
            
            found_canonical = None
            found_field = None
            
            for canonical_struct, canonical_field in canonical_map.values():
                if self.structures_match(structure, canonical_struct):
                    found_canonical = canonical_struct
                    found_field = canonical_field
                    break
            
            if found_canonical:
                if structure != found_canonical:
                    intermediate = None
                    if field_name == "SMILES":
                        mol = Chem.MolFromSmiles(structure)
                        if mol:
                            intermediate = Chem.MolToInchi(mol)
                            if intermediate != found_canonical:
                                self.changes_log.append(UnificationChange(
                                    original_structure=structure,
                                    intermediate_inchi=intermediate,
                                    canonical_inchi=found_canonical,
                                    structure_type=field_name
                                ))
                    else:
                        self.changes_log.append(UnificationChange(
                            original_structure=structure,
                            intermediate_inchi=None,
                            canonical_inchi=found_canonical,
                            structure_type=field_name
                        ))
                
                canonical_map[structure] = (found_canonical, found_field)
            else:
                canonical_struct = self.get_canonical_form(structure)
                
                if structure != canonical_struct:
                    intermediate = None
                    if field_name == "SMILES":
                        mol = Chem.MolFromSmiles(structure)
                        if mol:
                            intermediate = Chem.MolToInchi(mol)
                            if intermediate != canonical_struct:
                                self.changes_log.append(UnificationChange(
                                    original_structure=structure,
                                    intermediate_inchi=intermediate,
                                    canonical_inchi=canonical_struct,
                                    structure_type=field_name
                                ))
                    else:
                        self.changes_log.append(UnificationChange(
                            original_structure=structure,
                            intermediate_inchi=None,
                            canonical_inchi=canonical_struct,
                            structure_type=field_name
                        ))
                
                canonical_map[structure] = (canonical_struct, field_name)
        
        modified_entries = []
        for entry in entries:
            entry_copy = entry.copy()
            structure, field_name = self.extract_structure(entry_copy)
            
            if structure and structure in canonical_map:
                canonical_struct, canonical_field = canonical_map[structure]
                
                if "SMILES" in entry_copy:
                    del entry_copy["SMILES"]
                if "INCHI" in entry_copy:
                    del entry_copy["INCHI"]
                
                entry_copy["INCHI"] = canonical_struct
            
            modified_entries.append(entry_copy)

        return modified_entries
      
    def cross_unify(
        self,
        entries_a: List[Dict],
        entries_b: List[Dict],
        source_a: str = "File A",
        source_b: str = "File B"
    ) -> List[Dict]:
        canonical_from_a = {}
        for entry in entries_a:
            structure, field_name = self.extract_structure(entry)
            if structure:
                canonical_from_a[structure] = (structure, field_name)
        
        modified_b = []
        for entry in entries_b:
            entry_copy = entry.copy()
            structure_b, field_b = self.extract_structure(entry_copy)
            
            if structure_b:
                matched_canonical = None
                matched_field = None
                
                for struct_a, (canonical_a, field_a) in canonical_from_a.items():
                    if self.structures_match(structure_b, struct_a):
                        matched_canonical = struct_a
                        matched_field = field_a
                        break
                
                if matched_canonical:
                    if structure_b != matched_canonical:
                        intermediate = None
                        if field_b == "SMILES":
                            mol = Chem.MolFromSmiles(structure_b)
                            if mol:
                                intermediate = Chem.MolToInchi(mol)
                                if intermediate != matched_canonical:
                                    self.changes_log.append(UnificationChange(
                                        original_structure=structure_b,
                                        intermediate_inchi=intermediate,
                                        canonical_inchi=matched_canonical,
                                        structure_type=field_b
                                    ))
                        else:
                            self.changes_log.append(UnificationChange(
                                original_structure=structure_b,
                                intermediate_inchi=None,
                                canonical_inchi=matched_canonical,
                                structure_type=field_b
                            ))
                    
                    if "SMILES" in entry_copy:
                        del entry_copy["SMILES"]
                    if "INCHI" in entry_copy:
                        del entry_copy["INCHI"]
                    entry_copy["INCHI"] = matched_canonical  # Siempre InChI
            
            modified_b.append(entry_copy)
        
        return entries_a + modified_b
    
    def write_mgf(self, entries: List[Dict], output_path: str):
        with open(output_path, "w") as f:
            for entry in entries:
                f.write("BEGIN IONS\n")
                
                for key, value in entry.items():
                    if key != "_PEAKS":
                        f.write(f"{key}={value}\n")
                
                if "_PEAKS" in entry:
                    for peak_line in entry["_PEAKS"]:
                        f.write(f"{peak_line}\n")
                
                f.write("END IONS\n\n")
    
    def write_log(self, output_path: str):
        log_data = {
            "level": self.level,
            "total_changes": len(self.changes_log),
            "changes": [asdict(change) for change in self.changes_log]
        }
        
        with open(output_path, "w") as f:
            json.dump(log_data, f, indent=2)
    
    def process_files(
        self, file_path_a: str, file_path_b: str,
        output_mgf: Optional[str] = None, output_log: Optional[str] = None) -> Dict:
        entries_a = self.parse_mgf(file_path_a)
        entries_b = self.parse_mgf(file_path_b)
        
        source_a = Path(file_path_a).name
        source_b = Path(file_path_b).name
        
        changes_before = len(self.changes_log)
        entries_a_unified = self.unify_inchis_in_file(entries_a, source_a)
        changes_a = len(self.changes_log) - changes_before
        
        changes_before = len(self.changes_log)
        entries_b_unified = self.unify_inchis_in_file(entries_b, source_b)
        changes_b = len(self.changes_log) - changes_before
        
        changes_before = len(self.changes_log)
        all_entries = self.cross_unify(
            entries_a_unified, entries_b_unified, source_a, source_b
        )
        changes_cross = len(self.changes_log) - changes_before
        
        if output_mgf:
            self.write_mgf(all_entries, output_mgf)

        if output_log:
            self.write_log(output_log)
            
        return {
            "input_counts": {
                source_a: len(entries_a),
                source_b: len(entries_b)
            },
            "output_count": len(all_entries),
            "changes_count": len(self.changes_log),
            "changes_breakdown": {
                f"{source_a}_internal": changes_a,
                f"{source_b}_internal": changes_b,
                f"{source_b}_to_{source_a}": changes_cross
            },
            "changes_log": [asdict(change) for change in self.changes_log],
            "level": self.level
        }