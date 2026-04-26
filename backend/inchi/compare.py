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
    return [line.strip() for line in Path(file_path).read_text().splitlines() if line.strip()]


def _extract_differences(comparison):
    return {
        k.name: False
        for k, v in comparison.items()
        if not bool(v)
    }


def compare_text_files(file1, file2, config, mode="pairwise", only_differences=False):
    list1 = _read_file_lines(file1)
    list2 = _read_file_lines(file2)

    results = []

    def process(i1, i2):
        comparison = InChI.get_ids(i1, i2, config)

        if only_differences:
            diffs = _extract_differences(comparison)
            if not diffs:
                return None

            return {
                "inchi_1": i1,
                "inchi_2": i2,
                "differences": diffs
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
            if res:
                results.append(res)

    elif mode == "cross":
        for i1 in list1:
            for i2 in list2:
                res = process(i1, i2)
                if res:
                    results.append(res)

    return {"comparisons": results}


def compare_mgf_files(file1, file2, config, only_differences=False):
    entries1 = MgfParser.parse_mgf(file1)
    entries2 = MgfParser.parse_mgf(file2)

    list1 = MgfParser.extract_inchis(entries1)
    list2 = MgfParser.extract_inchis(entries2)

    results = []

    for i1 in list1:
        for i2 in list2:
            comparison = InChI.get_ids(i1, i2, config)

            if only_differences:
                diffs = _extract_differences(comparison)
                if not diffs:
                    continue

                results.append({
                    "inchi_1": i1,
                    "inchi_2": i2,
                    "differences": diffs
                })
            else:
                results.append({
                    "inchi_1": i1,
                    "inchi_2": i2,
                    "results": {k.name: bool(v) for k, v in comparison.items()}
                })

    return {"comparisons": results}