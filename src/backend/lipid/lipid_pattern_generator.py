import openpyxl
from pathlib import Path
from typing import Dict, List
from dataclasses import dataclass

@dataclass
class HeadgroupPattern:
    """
    Lipid headgroup pattern definition.
    
    This is defined here (duplicated from lipid_analysis.py) 
    to avoid circular import issues.
    """
    name: str
    smarts: str
    lipid_class: str
    fa_positions: List[str]
    description: str = ""


class PatternGenerator:
    """
    Generate lipid SMARTS patterns from Excel templates and sugars.
    """
    
    @staticmethod
    def load_templates_from_excel(excel_path: str) -> Dict[str, dict]:
        """
        Load templates from "Hoja1" column D.
        
        Returns:
            Dict of {template_id: {smarts, lipid_class, fa_positions, name}}
        """
        templates = {}
        
        try:
            wb = openpyxl.load_workbook(excel_path)
            ws = wb["Hoja1"]
            
            for row_idx, row in enumerate(ws.iter_rows(min_row=2, values_only=True), start=2):
                if not row or len(row) < 5:
                    continue
                
                example_num = row[0]  
                template_smarts = row[3] 
                lipid_type = row[4]  
                lipid_subtype = row[5] 
                
                if not template_smarts or not isinstance(template_smarts, str):
                    continue
                
                if "[G]" not in template_smarts:
                    continue
                
                template_id = f"template_{example_num}".replace(" ", "_").replace("(", "").replace(")", "")
                
                templates[template_id] = {
                    "smarts": template_smarts,
                    "lipid_class": lipid_type if lipid_type else "Unknown",
                    "lipid_subtype": lipid_subtype if lipid_subtype else "",
                    "fa_positions": ["N-acyl"] if "SP" in str(lipid_type) else ["sn-1", "sn-2"],  
                    "example_num": example_num
                }
            
            print(f"✓ Loaded {len(templates)} templates from Excel")
            
        except Exception as e:
            print(f"Error loading templates: {e}")
        
        return templates
    
    @staticmethod
    def load_sugars_from_excel(excel_path: str, sheet_name: str = None) -> Dict[str, str]:
        """
        Load sugar library from "Azúcares (G)" sheet.
        
        Returns:
            Dictionary: {sugar_name: sugar_smarts}
        """
        sugars = {}
        
        try:
            wb = openpyxl.load_workbook(excel_path)
            
            if sheet_name is None:
                for name in ["Azúcares (G)", "Azúcares", "Sugars", "azucares"]:
                    if name in wb.sheetnames:
                        sheet_name = name
                        break
            
            if sheet_name not in wb.sheetnames:
                print(f"Warning: Sugar sheet not found. Available: {wb.sheetnames}")
                return {}
            
            ws = wb[sheet_name]
            
            for row in ws.iter_rows(min_row=2, values_only=True):
                if not row or not row[0]:
                    continue
                
                sugar_name = str(row[0]).strip()
                sugar_smarts = str(row[1]).strip() if len(row) > 1 and row[1] else None
                
                if sugar_name and sugar_smarts:
                    sugars[sugar_name] = sugar_smarts
            
            if sugars:
                sugar_list = ', '.join(list(sugars.keys())[:5])
                print(f"✓ Loaded {len(sugars)} sugars: {sugar_list}...")
            
        except Exception as e:
            print(f"Error reading sugars: {e}")
        
        return sugars
    
    @staticmethod
    def generate_patterns(templates: Dict, sugars: Dict[str, str]) -> Dict[str, HeadgroupPattern]:
        """
        Generate all pattern combinations: templates × sugars.
        
        Key insight: Replace FIRST [G] in template with sugar.
        Result may still have [G] from the sugar (those are sugar-sugar links).
        
        Args:
            templates: From load_templates_from_excel()
            sugars: From load_sugars_from_excel()
        
        Returns:
            {pattern_id: HeadgroupPattern}
        """
        patterns = {}
        
        for template_id, template_info in templates.items():
            template_smarts = template_info["smarts"]
            
            if "[G]" not in template_smarts:
                continue
            
            for sugar_name, sugar_smarts in sugars.items():
                complete_smarts = template_smarts.replace("[G]", sugar_smarts, 1)
                
                safe_sugar_name = sugar_name.replace(" ", "_").replace("(", "").replace(")", "").replace("‑", "_")
                pattern_id = f"{template_id}_{safe_sugar_name}"
                
                lipid_subtype = template_info.get("lipid_subtype", "")
                pattern_name = f"{lipid_subtype} ({sugar_name})" if lipid_subtype else f"Example {template_info['example_num']} ({sugar_name})"
                
                patterns[pattern_id] = HeadgroupPattern(
                    name=pattern_name,
                    smarts=complete_smarts,
                    lipid_class=template_info["lipid_class"],
                    fa_positions=template_info["fa_positions"],
                    description=f"Auto-generated from template {template_info['example_num']} + {sugar_name}"
                )
        
        return patterns


def build_combined_patterns(manual_patterns: Dict[str, HeadgroupPattern],
                           excel_path: str = None) -> Dict[str, HeadgroupPattern]:
    """
    Combine manual patterns + auto-generated patterns from Excel.
    
    This is what LipidHeadValidator will call.
    
    Args:
        manual_patterns: Existing manual patterns
        excel_path: Path to Excel file (optional)
    
    Returns:
        Complete pattern dictionary
    """
    all_patterns = dict(manual_patterns)
    
    if excel_path and Path(excel_path).exists():
        try:
            templates = PatternGenerator.load_templates_from_excel(excel_path)
            sugars = PatternGenerator.load_sugars_from_excel(excel_path)
            
            if templates and sugars:
                generated = PatternGenerator.generate_patterns(templates, sugars)
                
                added = 0
                for pattern_id, pattern in generated.items():
                    if pattern_id not in all_patterns:
                        all_patterns[pattern_id] = pattern
                        added += 1
                
                print(f"[OK] Total patterns: {len(all_patterns)} "
                      f"({len(manual_patterns)} manual + {added} generated)")
        except Exception as e:
            print(f"Pattern generation failed: {e}")
            import traceback
            traceback.print_exc()
            print(f"Using manual patterns only ({len(manual_patterns)} total)")
    else:
        if excel_path:
            print(f"Excel not found: {excel_path}")
        print(f"Using manual patterns only ({len(all_patterns)} total)")
    
    return all_patterns