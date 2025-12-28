from enum import Enum

class InchiLayers(Enum):
    COMPLETE_IDENTITY = "complete_identity"
    INDEPENDENT_SALTS = "independent_salts"
    INDEPENDENT_CHARGES = "independent_charges"
    TAUTOMERIC = "tautomeric"
    TAUTOMERIC_1 = "tautomeric_1"
    TAUTOMERIC_2 = "tautomeric_2"
    STEREOCHEMICAL = "stereochemical" 
    STEREOCHEMICAL_CIS_TRANS = "stereochemical_cis_trans" 
    INDEPENDENT_DOUBLE_BONDS = "independent_double_bonds"
    ISOTOPIC = "isotopic" #presencia heavy metals

#TODO: tautomeric inchi.exe    
#Definimos tautomeric 1 y 2 (definidos en inchi.exe: ver cómo funciona)
#para carga: si es p+N quitamos, si es q+N lo dejamos! , si es p-N o q-N tenemos que neutralizar el inchi usando RDKit
