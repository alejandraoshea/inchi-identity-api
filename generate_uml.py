import zlib
import urllib.request
import os

def encode6bit(b):
    if b < 10: return chr(48 + b)
    b -= 10
    if b < 26: return chr(65 + b)
    b -= 26
    if b < 26: return chr(97 + b)
    b -= 26
    if b == 0: return '-'
    if b == 1: return '_'
    return '?'

def append3bytes(b1, b2, b3):
    c1 = b1 >> 2
    c2 = ((b1 & 0x3) << 4) | (b2 >> 4)
    c3 = ((b2 & 0xF) << 2) | (b3 >> 6)
    c4 = b3 & 0x3F
    return encode6bit(c1) + encode6bit(c2) + encode6bit(c3) + encode6bit(c4)

def encode_plantuml(data: bytes) -> str:
    result = ""
    i = 0
    while i < len(data):
        b1 = data[i]
        b2 = data[i+1] if i+1 < len(data) else 0
        b3 = data[i+2] if i+2 < len(data) else 0
        result += append3bytes(b1, b2, b3)
        i += 3
    return result

def plantuml_to_png(uml_text: str, output_path: str):
    data = uml_text.encode('utf-8')
    obj = zlib.compressobj(9, zlib.DEFLATED, -15)
    compressed = obj.compress(data) + obj.flush()
    encoded = encode_plantuml(compressed)
    url = f"https://www.plantuml.com/plantuml/png/{encoded}"
    print(f"Fetching ({len(compressed)} bytes compressed)...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        img_data = resp.read()
    with open(output_path, "wb") as f:
        f.write(img_data)
    print(f"  -> {output_path}  ({len(img_data)} bytes)")

ROOT = r"c:\Users\34611\Desktop\INGENIERIA_BIOMEDICA\CUARTO_CARRERA\TFG\identity-levels-inchi"

UML = """
@startuml uml_diagram
skinparam dpi 150
skinparam defaultFontSize 13
skinparam classFontSize 14
skinparam classFontStyle bold
skinparam classAttributeFontSize 12
skinparam classStereotypeFontSize 11
skinparam packageFontSize 15
skinparam packageFontStyle bold
skinparam ArrowFontSize 11
skinparam classBackgroundColor #FEFEFE
skinparam classBorderColor #44607A
skinparam classArrowColor #2D6A9F
skinparam classBorderThickness 1.5
skinparam nodesep 60
skinparam ranksep 70
hide empty members

left to right direction

' ══════════════════════════════════════════════
'  INCHI PACKAGE
' ══════════════════════════════════════════════
package "inchi" #E3F2FD {

  enum InchiLayers {
    COMPLETE_IDENTITY
    ISOTOPIC_INDEPENDENCE
    SALTS_INDEPENDENCE
    CHARGES_INDEPENDENCE
    DOUBLE_BONDS_INDEPENDENCE
    STEREO_CIS_TRANS_INDEPENDENCE
    TAUTOMER_INDEPENDENCE
    SUBSTITUENT_POSITION_INDEPENDENCE
  }

  class InChI {
    {static} isCompleteIdentity(inchi1, inchi2) : bool
    {static} normalize_input(structure) : str
    {static} mol_from_inchi(inchi) : Mol
    {static} areEqualNoIsotopes(inchi1, inchi2) : bool
    {static} main_fragment(mol) : Mol
    {static} areEqualDisolvedSalts(inchi1, inchi2) : bool
    {static} has_negative_charge(inchi) : bool
    {static} get_charge_info(inchi) : dict
    {static} remove_only_p_layer(inchi) : str
    {static} neutralize_molecule(mol) : Mol
    {static} areEqualNoCharges(inchi1, inchi2) : bool
    {static} areEqualNoStereo(inchi1, inchi2) : bool
    {static} remove_cis_trans(mol) : Mol
    {static} areEqualNoPositionDoubleBond(inchi1, inchi2) : bool
    {static} run_inchitrust(mol, inchitrust_path) : str
    {static} areEqualTautomers(inchi1, inchi2, inchitrust_path) : bool
    {static} get_substituent_signatures(mol) : list
    {static} substituent_position_independent_signature(inchi) : str
    {static} areEqualSubstituentIndependent(inchi1, inchi2) : bool
    {static} get_ids(inchi1, inchi2, config) : dict
  }

  class InChIParser {
    {static} getMainLayer(inchi) : str
    {static} getAtomConnectionsSublayer(inchi)
    {static} getHydrogenAtomsSublayer(inchi)
    {static} getChargeSublayer(inchi)
    {static} getProtonSublayer(inchi)
    {static} removeChargeLayers(inchi) : str
    {static} getDoubleBondsSublayer(inchi)
    {static} getTetrahedralStereoSublayer(inchi)
    {static} getTypeStereoInfoSublayer(inchi)
    {static} getIsotopicLayer(inchi)
    {static} getIsotopicHydrogenSublayer(inchi)
    {static} getIsotopicStereoSublayer(inchi)
    {static} getFixedHLayer(inchi)
    {static} getReconnectedLayer(inchi)
    {static} removeDoubleBondsSublayer(inchi) : str
    {static} getTautomerLayer(inchi)
    {static} removeStereoLayers(inchi) : str
    {static} removeIsotopicLayers(inchi) : str
  }

  InChI ..> InChIParser : uses
  InChI ..> InchiLayers : uses
}

' ══════════════════════════════════════════════
'  LIPID PACKAGE
' ══════════════════════════════════════════════
package "lipid" #E8F5E9 {

  class HeadgroupPattern <<dataclass>> {
    + name : str
    + smarts : str
    + lipid_class : str
    + fa_positions : List[str]
    + description : str
  }

  class PatternGenerator {
    {static} load_templates_from_excel(excel_path) : Dict
    {static} load_sugars_from_excel(excel_path, sheet_name) : Dict
    {static} generate_patterns(templates, sugars) : Dict[str, HeadgroupPattern]
  }

  class LipidHeadValidator {
    + MANUAL_PATTERNS : Dict[str, HeadgroupPattern]
    + HEADGROUP_PATTERNS : Dict[str, HeadgroupPattern]
    + compiled_patterns : Dict
    --
    + __init__()
    + compile_patterns()
    + matches_pattern(mol, pattern_id) : bool
    + matches_any_valid_head(mol) : bool
    + get_matching_patterns(mol) : List[HeadgroupPattern]
    + identify_lipid_class(mol) : List[str]
    + validate_structure(mol, verbose) : Dict
    {static} is_valid_lipid_structure(mol) : bool
  }

  class LipidAnalysis {
    {static} CLASSYFIRE_URL : str
    {static} MIN_TAIL_CARBONS : int
    + HEAD_ANCHORS : Dict
    --
    {static} is_lipid_rdkit(mol) : bool
    {static} is_lipid_classyfire(inchi) : bool
    {static} is_lipid(inchi, mol, use_classyfire) : bool
    {static} has_long_carbon_chain(mol, min_len) : bool
    {static} count_lipid_linkages(mol)
    {static} tail_sig_levelB(t)
    {static} tail_sig_levelC(t)
    {static} atom_count(mol)
  }

  class TailExtractor {
    {static} detect_head_atoms(mol)
    {static} extract_tails(mol) : list
    {static} walk_chain(start_atom, coming_from) : list
  }

  class LipidComparator {
    {static} lipid_signature(mol)
  }

  PatternGenerator ..> HeadgroupPattern      : creates
  LipidHeadValidator "1" *-- "many" HeadgroupPattern : contains
  LipidHeadValidator ..> LipidAnalysis       : uses
  LipidHeadValidator ..> PatternGenerator    : uses
  LipidComparator    ..> LipidAnalysis       : uses
  LipidComparator    ..> TailExtractor       : uses
}

' ══════════════════════════════════════════════
'  CROSS-PACKAGE RELATIONSHIPS
' ══════════════════════════════════════════════
InChI ..> LipidAnalysis  : uses
InChI ..> TailExtractor  : uses

@enduml
"""

plantuml_to_png(UML.strip(), os.path.join(ROOT, "uml_diagram.png"))
print("Done.")
