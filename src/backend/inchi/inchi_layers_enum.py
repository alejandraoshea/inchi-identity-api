from enum import Enum

class InchiLayers(Enum):
    COMPLETE_IDENTITY = "complete_identity"
    ISOTOPIC_INDEPENDENCE = "isotopic_independent"
    SALTS_INDEPENDENCE = "salt_independent"
    CHARGES_INDEPENDENCE = "charge_independent"
    STEREOCHEMICAL_INDEPENDENCE = "setereochemical_independent"
    DOUBLE_BONDS_INDEPENDENCE = "double_bond_position_independent"
    STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE = "cis_trans_independent" 
    TAUTOMER_INDEPENDENCE = "tautomer_independent"
    SUBSTITUENT_POSITION_INDEPENDENCE = "substituent_position_independent"
