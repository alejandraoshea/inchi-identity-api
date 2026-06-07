from flask import Blueprint, jsonify, request
from rdkit import Chem
from rdkit.Chem import AllChem
from inchi_identity.inchi.smiles_pattern import SmilesCorrector
from dotenv import load_dotenv

load_dotenv()

generation_3d_routes = Blueprint("generation_3d_routes", __name__)


def normalize_to_inchi(input_str: str) -> str:
    if not input_str:
        raise ValueError("Empty input")

    input_str = input_str.strip()

    if input_str.startswith("InChI="):
        return input_str

    try:
        correction = SmilesCorrector.auto_correct(input_str, verbose=False)
        corrected_smiles = correction["corrected"]

        mol = Chem.MolFromSmiles(corrected_smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {input_str}")

        inchi = Chem.MolToInchi(mol)
        if not inchi:
            raise ValueError(f"Could not convert SMILES to InChI: {input_str}")

        return inchi

    except Exception as e:
        raise ValueError(f"Invalid input (not SMILES or InChI): {input_str} - {e}")


def inchi_to_mol3d(inchi):
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        return None, None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        return None, None
    AllChem.UFFOptimizeMolecule(mol)
    return mol, Chem.MolToMolBlock(mol)


@generation_3d_routes.route("/api/generate_3d", methods=["POST"])
def generate_3d():
    """
    Accepts an InChI or SMILES string, generates a 3D molecular conformation
    using RDKit, and returns the structure serialized in SDF format.
    All rendering is handled client-side by the frontend.
    """
    try:
        input_str = request.json.get("inchi")

        inchi = normalize_to_inchi(input_str)

        mol, sdf = inchi_to_mol3d(inchi)
        if mol is None:
            return jsonify({"error": "Could not generate 3D conformation for the provided structure"}), 400

        return jsonify({"sdf": sdf, "inchi": inchi})

    except Exception as e:
        return jsonify({"error": str(e)}), 500