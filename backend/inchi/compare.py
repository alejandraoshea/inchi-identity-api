import json
from pathlib import Path
from backend.inchi.determine_levels_id import InChI
from backend.inchi.mgf_parser import MgfParser


def compare_pair(inchi1, inchi2, config):
    comparison = InChI.get_ids(inchi1, inchi2, config)

    return {
        "inchi_1": inchi1,
        "inchi_2": inchi2,
        "results": {k.name: bool(v) for k, v in comparison.items()}
    }


def _read_file_lines(file_path):
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
        comparison = InChI.get_ids(i1, i2, config)
        matches = extract_matches(comparison)

        if not matches:
            return None

        return {
            "inchi_1": i1,
            "inchi_2": i2,
            "matches": matches
        }

    if mode == "pairwise":
        for i1, i2 in zip(list1, list2):
            res = process(i1, i2)
            if res:
                results.append(res)

    elif mode == "cross":
        for i1 in list1:
            for i2 in list2:
                res = process(i1, i2)
                if res:
                    results.append(res)

    return {"comparisons": results}


def compare_mgf_files(file1, file2, config, level="COMPLETE_IDENTITY"):
    entries1 = MgfParser.parse_mgf(file1)
    entries2 = MgfParser.parse_mgf(file2)

    structs1 = MgfParser.extract_structures(entries1)
    structs2 = MgfParser.extract_structures(entries2)

    all_structs = structs1 + structs2

    groups = []

    for item in all_structs:
        placed = False

        for group in groups:
            rep = group["representative"]

            comparison = InChI.get_ids(
                item["structure"],
                rep,
                config
            )

            if comparison.get(level):
                group["entries"].append(item["entry"])
                placed = True
                break

        if not placed:
            groups.append({
                "representative": item["structure"],
                "entries": [item["entry"]]
            })

    return {"groups": groups}
