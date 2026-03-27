from flask import Blueprint, jsonify, request
from inchi.determine_levels_id import InChi
from inchi.compare import compare_files
from inchi.config_loader import load_config
import tempfile

inchi_comparison_routes = Blueprint("inchi_comparison_routes", __name__)

@inchi_comparison_routes.route("/api/compare_inchis", methods=["POST"])
def compare_inchis():
    data = request.json

    inchi1 = data.get("inchi1")
    inchi2 = data.get("inchi2")

    if not inchi1 or not inchi2:
        return jsonify({"message": "Both InChIs are required"}), 400

    if not inchi1.startswith("InChI=") or not inchi2.startswith("InChI="):
        return jsonify({"message": "Invalid InChI format"}), 400
    
    config = load_config()

    results = InChi.get_ids(inchi1, inchi2, config)

    formatted = {}
    for k, v in results.items():
        formatted[k.name] = v if isinstance(v, dict) else bool(v)

    return jsonify({
        "inchi_1": inchi1,
        "inchi_2": inchi2,
        "results": formatted
    })

@inchi_comparison_routes.route("/api/compare_inchis_custom", methods=["POST"])
def compare_inchis_custom():
    data = request.json

    inchi1 = data.get("inchi1")
    inchi2 = data.get("inchi2")
    levels = data.get("levels", [])

    if not inchi1 or not inchi2:
        return jsonify({"message": "Both InChIs are required"}), 400

    if not inchi1.startswith("InChI=") or not inchi2.startswith("InChI="):
        return jsonify({"message": "Invalid InChI format"}), 400
    
    config = load_config()

    # disable everything first
    for section in config["identity_criteria"].values():
        for key in section:
            if isinstance(section[key], bool):
                section[key] = False

    # mapping (cleaner than many ifs)
    mapping = {
        "complete_identity": ("complete_identity", "enabled"),
        "isotope": ("isotope_independence", "isotope_independent_identity"),
        "salt": ("salt_independence", "desalted_identity"),
        "charge": ("charge_independence", "charge_independent_identity"),
        "double_bond": ("isomer_independence", "double_bond_position_independent_identity"),
        "cis_trans": ("isomer_independence", "cis_trans_independent_identity"),
        "tautomer": ("tautomer_independence", "tautomer_independent_identity"),
        "substituent": ("tautomer_independence", "substituent_position_independent_identity"),
    }

    if not levels:
        return jsonify({"message": "No levels selected"}), 400

    for level in levels:
        if level in mapping:
            section, key = mapping[level]
            config["identity_criteria"][section][key] = True

    results = InChi.get_ids(inchi1, inchi2, config)

    formatted = {}
    for k, v in results.items():
        formatted[k.name] = v if isinstance(v, dict) else bool(v)

    return jsonify({"results": formatted})

@inchi_comparison_routes.route("/api/compare_inchi_files", methods=["POST"])
def compare_inchi_files():
    file1 = request.files.get("file1")
    file2 = request.files.get("file2")

    if not file1 or not file2:
        return jsonify({"message": "Both files are required"}), 400

    try:
        with tempfile.NamedTemporaryFile(delete=False) as f1:
            file1.save(f1.name)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(delete=False) as f2:
            file2.save(f2.name)
            path2 = f2.name

        result = compare_files(path1, path2)

        return jsonify(result)

    except Exception as e:
        return jsonify({"message": str(e)}), 500
    
@inchi_comparison_routes.route("/api/inchi_levels", methods=["GET"])
def get_inchi_levels():
    return jsonify([
        {"key": "complete_identity", "label": "Complete Identity"},
        {"key": "isotope", "label": "Isotope Independence"},
        {"key": "salt", "label": "Salt Independence"},
        {"key": "charge", "label": "Charge Independence"},
        {"key": "double_bond", "label": "Double Bond Independence"},
        {"key": "cis_trans", "label": "Cis/Trans Independence"},
        {"key": "tautomer", "label": "Tautomer Independence"},
        {"key": "substituent", "label": "Substituent Independence"}
    ])