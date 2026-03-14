import json
from pathlib import Path
from inchi.config_loader import load_config
from .determine_levels_id import InChi

def compare_files(file1, file2, config_path=None, output_file=None):
    config = load_config(config_path)

    with open(file1) as f:
        inchi_list1 = [line.strip() for line in f if line.strip()]

    with open(file2) as f:
        inchi_list2 = [line.strip() for line in f if line.strip()]

    results = []

    for i1 in inchi_list1:
        for i2 in inchi_list2:
            comparison = InChi.get_ids(i1, i2, config)

            results.append({
                "inchi_1": i1,
                "inchi_2": i2,
                "results": {k.name: v for k, v in comparison.items()}
            })

    output = {
        "comparisons": results
    }

    if output_file:
        with open(output_file, "w") as f:
            json.dump(output, f, indent=4)

    return output