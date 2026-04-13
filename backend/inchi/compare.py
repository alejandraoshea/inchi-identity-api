import json
from pathlib import Path
from backend.inchi.config_loader import load_config
from backend.inchi.determine_levels_id import InChI
from backend.inchi.mgf_parser import MgfParser

def compare_pair(inchi1, inchi2, config):
    comparison = InChI.get_ids(inchi1, inchi2, config)

    return {
        "inchi_1": inchi1,
        "inchi_2": inchi2,
        "results": {k.name: v for k, v in comparison.items()}
    }

def compare_text_files(list1, list2, config):
    results = []

    for i1 in list1:
        for i2 in list2:
            comparison = InChI.get_ids(i1, i2, config)

            results.append({
                "inchi_1": i1,
                "inchi_2": i2,
                "results": {k.name: v for k, v in comparison.items()}
            })

    return {"comparisons": results}

def compare_mgf_files(file1, file2, config):
    entries1 = MgfParser.parse_mgf(file1)
    entries2 = MgfParser.parse_mgf(file2)

    list1 = MgfParser.extract_inchis(entries1)
    list2 = MgfParser.extract_inchis(entries2)

    results = []

    for i1 in list1:
        for i2 in list2:
            comparison = InChI.get_ids(i1, i2, config)

            results.append({
                "inchi_1": i1,
                "inchi_2": i2,
                "results": {k.name: v for k, v in comparison.items()}
            })

    return {"comparisons": results}