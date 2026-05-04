from flask import Blueprint, jsonify, request
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from dotenv import load_dotenv
import base64

load_dotenv()

generation_3d_routes = Blueprint("generation_3d_routes", __name__)


def _inchi_to_mol3d(inchi):
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
    try:
        inchi = request.json.get("inchi")
        mol, sdf = _inchi_to_mol3d(inchi)
        if mol is None:
            return jsonify({"error": "Could not generate 3D molecule"}), 400
        return jsonify({"sdf": sdf})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@generation_3d_routes.route("/api/render_3d_image", methods=["POST"])
def render_3d_image():
    try:
        inchi  = request.json.get("inchi")
        width  = int(request.json.get("width",  400))
        height = int(request.json.get("height", 260))

        mol, sdf = _inchi_to_mol3d(inchi)
        if mol is None:
            return jsonify({"error": "Could not generate 3D molecule"}), 400

        mol_draw = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(mol_draw)

        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.DrawMolecule(mol_draw)
        drawer.FinishDrawing()

        b64 = base64.b64encode(drawer.GetDrawingText()).decode("utf-8")
        return jsonify({"image": b64, "sdf": sdf})

    except Exception as e:
        return jsonify({"error": str(e)}), 500