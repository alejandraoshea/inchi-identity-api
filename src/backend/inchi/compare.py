import json
from pathlib import Path
from backend.inchi.determine_levels_id import InChI
from backend.inchi.smiles_pattern import SmilesCorrector  # ← NEW
from backend.parsers.mgf_parser import MgfParser, SimpleMgfDeduplicator
from typing import Dict, Optional


def compare_pair(inchi1, inchi2, config):
    if not inchi1.startswith("InChI="):
        correction = SmilesCorrector.auto_correct(inchi1, verbose=True)
        if correction["parse_result"] == "error":
            return {
                "inchi_1": inchi1,
                "inchi_2": inchi2,
                "error": f"Invalid SMILES: {correction['message']}",
                "results": {}
            }
        inchi1 = correction["corrected"]
    
    if not inchi2.startswith("InChI="):
        correction = SmilesCorrector.auto_correct(inchi2, verbose=True)
        if correction["parse_result"] == "error":
            return {
                "inchi_1": inchi1,
                "inchi_2": inchi2,
                "error": f"Invalid SMILES: {correction['message']}",
                "results": {}
            }
        inchi2 = correction["corrected"]
    
    inchi1 = InChI.normalize_input(inchi1)
    inchi2 = InChI.normalize_input(inchi2)
    
    comparison = InChI.get_ids(inchi1, inchi2, config)

    return {
        "inchi_1": inchi1,
        "inchi_2": inchi2,
        "results": {k.name: bool(v) for k, v in comparison.items()}
    }


def read_file_lines(file_path):
    return [
        line.strip()
        for line in Path(file_path).read_text().splitlines()
        if line.strip()
    ]


def extract_matches(comparison):
    return {
        k.name: True
        for k, v in comparison.items()
        if bool(v)
    }


def compare_text_files(list1, list2, config, mode="pairwise", only_equal=False):
    results = []

    def process(i1, i2):
        if not i1.startswith("InChI="):
            correction = SmilesCorrector.auto_correct(i1)
            if correction["parse_result"] != "error":
                i1 = correction["corrected"]
        
        if not i2.startswith("InChI="):
            correction = SmilesCorrector.auto_correct(i2)
            if correction["parse_result"] != "error":
                i2 = correction["corrected"]
        
        i1 = InChI.normalize_input(i1)
        i2 = InChI.normalize_input(i2)
        
        comparison = InChI.get_ids(i1, i2, config)

        if only_equal:
            matches = extract_matches(comparison)
            if not matches:
                return None
            return {
                "inchi_1": i1,
                "inchi_2": i2,
                "matches": matches
            }
        else:
            return {
                "inchi_1": i1,
                "inchi_2": i2,
                "results": {k.name: bool(v) for k, v in comparison.items()}
            }

    if mode == "pairwise":
        for i1, i2 in zip(list1, list2):
            res = process(i1, i2)
            if res is not None:
                results.append(res)

    elif mode == "cross":
        for i1 in list1:
            for i2 in list2:
                res = process(i1, i2)
                if res is not None:
                    results.append(res)

    return {"comparisons": results}


def compare_mgf_files(
    file1: str,
    file2: str,
    config: dict,
    level: str = "COMPLETE_IDENTITY",
    output_mgf: Optional[str] = None,
    output_log: Optional[str] = None
) -> Dict:
    deduplicator = SimpleMgfDeduplicator(level=level, config=config)
    
    return deduplicator.process_files(
        file1,
        file2,
        output_mgf=output_mgf,
        output_log=output_log
    )