from input.input_inchi import InputInChI
from inchi.determine_levels_id import InChI
from inchi.config_loader import load_config
from inchi import InChI
from input.input_inchi import InputInChI

config = load_config()

InputInChI.input_inchi(InChI.get_ids, config)
InputInChI.input_inchi(InChI.get_ids)
