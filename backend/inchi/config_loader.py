from pathlib import Path
import json
import os

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parents[2] / "configs/default_config.json"

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