from pathlib import Path
import json
import os

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parents[2] / "backend/configs/default_config.json"

def load_config(config_path=None):
    if config_path is None:
        config_path = DEFAULT_CONFIG_PATH

    with open(config_path) as f:
        return json.load(f)


def validate_inchitrust_path(path):
    if path is None:
        return None

    if path == "/path/to/inchitrust_executable":
        raise ValueError(
            "Please set the correct InChI Trust path in your config file."
        )

    if not os.path.isfile(path):
        raise ValueError(
            f"InChI Trust executable not found at: {path}"
        )

    return path

def build_config_from_levels(selected_levels, base_config):
    import copy
    config = copy.deepcopy(base_config)

    for section in config["identity_criteria"].values():
        for key, value in section.items():
            if isinstance(value, bool):
                section[key] = False

    mapping = {
        "complete_identity": ("complete_identity", "enabled"),
        "isotope": ("isotope_independence", "isotope_independent_identity"),
        "salt": ("salt_independence", "desalted_identity"),
        "charge": ("charge_independence", "charge_independent_identity"),
        "double_bond": ("isomer_independence", "double_bond_position_independent_identity"),
        "stereo_cis_trans": ("isomer_independence", "cis_trans_independent_identity"),
        "tautomer": ("tautomer_independence", "tautomer_independent_identity"),
        "substituent": ("tautomer_independence", "substituent_position_independent_identity"),
    }

    for lvl in selected_levels:
        if lvl in mapping:
            section, key = mapping[lvl]
            config["identity_criteria"][section][key] = True

    return config


def apply_inchitrust(config, path):
    path = validate_inchitrust_path(path)

    if path:
        config["identity_criteria"]["tautomer_independence"]["inchitrust_path"] = path

    return config