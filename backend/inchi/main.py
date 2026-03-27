from input.input_inchi import InputInChI
from inchi.determine_levels_id import InChi
from inchi.config_loader import load_config
from inchi import InChi
from input.input_inchi import InputInChI

config = load_config()

InputInChI.input_inchi(InChi.get_ids, config)
InputInChI.input_inchi(InChi.get_ids)
