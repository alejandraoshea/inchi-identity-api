from flask import Blueprint, jsonify, request
from backend.inchi.determine_levels_id import InChI
from backend.inchi.compare import compare_pair, compare_text_files, compare_mgf_files
from backend.inchi.config_loader import load_config, build_config_from_levels
import tempfile, traceback, os
from dotenv import load_dotenv

load_dotenv()

inchi_comparison_routes = Blueprint("inchi_comparison_routes", __name__)


@inchi_comparison_routes.route("/api/compare_inchis", methods=["POST"])
def compare_inchis():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"message": "No JSON received"}), 400
        inchi1 = data.get("inchi1")
        inchi2 = data.get("inchi2")
        if not inchi1 or not inchi2:
            return jsonify({"message": "Missing InChIs"}), 400
        config = load_config()
        result = compare_pair(inchi1, inchi2, config)
        return jsonify({"results": result["results"]})
    except Exception as e:
        return jsonify({"message": str(e)}), 500


@inchi_comparison_routes.route("/api/compare_inchis_custom", methods=["POST"])
def compare_inchis_custom():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"message": "No JSON received"}), 400
        inchi1 = data.get("inchi1")
        inchi2 = data.get("inchi2")
        selected_levels = data.get("levels", [])
        if not inchi1 or not inchi2:
            return jsonify({"message": "Missing InChIs"}), 400
        base_config = load_config()
        config = build_config_from_levels(selected_levels, base_config)
        result = compare_pair(inchi1, inchi2, config)
        return jsonify({"results": result["results"]})
    except Exception as e:
        return jsonify({"message": str(e)}), 500


@inchi_comparison_routes.route("/api/inchi_levels", methods=["GET"])
def get_inchi_levels():
    return jsonify([
        {"key": "complete_identity",  "label": "Complete Identity"},
        {"key": "isotope",            "label": "Isotope Independence"},
        {"key": "salt",               "label": "Salt Independence"},
        {"key": "charge",             "label": "Charge Independence"},
        {"key": "double_bond",        "label": "Double Bond Independence"},
        {"key": "cis_trans",          "label": "Cis/Trans Independence"},
        {"key": "tautomer",           "label": "Tautomer Independence"},
        {"key": "substituent",        "label": "Substituent Independence"}
    ])


@inchi_comparison_routes.route("/api/compare_files", methods=["POST"])
def compare_files_api():
    try:
        data = request.get_json()
        list1 = data.get("list1", [])
        list2 = data.get("list2", [])
        mode  = data.get("mode", "pairwise")
        if not list1 or not list2:
            return jsonify({"message": "Both lists required"}), 400
        config = load_config()
        result = compare_text_files(list1, list2, config, mode=mode, only_equal=True)
        return jsonify(result)
    except Exception as e:
        return jsonify({"message": str(e)}), 500


@inchi_comparison_routes.route("/api/compare_mgf_files", methods=["POST"])
def compare_mgf_files_upload():
    try:
        if "file1" not in request.files or "file2" not in request.files:
            return jsonify({"message": "Both MGF files required"}), 400

        file1 = request.files["file1"]
        file2 = request.files["file2"]
        level = request.form.get("level") or "COMPLETE_IDENTITY" 

        with tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".mgf") as tmp1:
            file1.save(tmp1.name)
            tmp1_path = tmp1.name

        with tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".mgf") as tmp2:
            file2.save(tmp2.name)
            tmp2_path = tmp2.name

        config = load_config()
        result = compare_mgf_files(tmp1_path, tmp2_path, config, level=level)

        os.unlink(tmp1_path)
        os.unlink(tmp2_path)

        return jsonify(result)

    except Exception as e:
        traceback.print_exc()
        return jsonify({"message": str(e)}), 500