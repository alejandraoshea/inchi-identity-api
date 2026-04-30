from flask import Blueprint, jsonify, request
from rdkit import Chem
from rdkit.Chem import AllChem
from dotenv import load_dotenv

load_dotenv()

generation_3d_routes = Blueprint("generation_3d_routes", __name__)

@generation_3d_routes.route("/api/generate_3d", methods=["POST"])
def generate_3d():
    try:
        inchi = request.json.get("inchi")

        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            return jsonify({"error": "Couldn't Generate the 3D Molecule"}), 400

        mol = Chem.AddHs(mol)

        #generate the 3d mol
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)

        sdf = Chem.MolToMolBlock(mol)

        return jsonify({"sdf": sdf})

    except Exception as e:
        return jsonify({"error": str(e)}), 500