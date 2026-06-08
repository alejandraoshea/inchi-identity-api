from flask import Blueprint, jsonify, request
from inchi_identity.inchi.determine_levels_id import InChI
from inchi_identity.inchi.compare import compare_pair, compare_text_files, compare_mgf_files
from inchi_identity.inchi.config_loader import load_config, build_config_from_layers
from importlib.resources import files as _pkg_files
import tempfile, traceback, os, base64
from dotenv import load_dotenv

load_dotenv()

_DEFAULT_CONFIG_PATH = _pkg_files('inchi_identity') / 'configs' / 'default_config.json'

inchi_comparison_routes = Blueprint("inchi_comparison_routes", __name__)

@inchi_comparison_routes.route("/api/compare_inchis", methods=["POST"])
def compare_inchis():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"message": "No JSON received"}), 400
        
        input1 = data.get("inchi1", "").strip()
        input2 = data.get("inchi2", "").strip()
        
        if not input1 or not input2:
            return jsonify({"message": "Missing InChIs"}), 400
        
        inchi1 = InChI.normalize_input(input1)
        inchi2 = InChI.normalize_input(input2)
        
        config = load_config(_DEFAULT_CONFIG_PATH)
        result = compare_pair(inchi1, inchi2, config)

        return jsonify({
            "inchi_1": inchi1,
            "inchi_2": inchi2,
            "results": result["results"]
        })

    except Exception as e:
        return jsonify({"message": str(e)}), 500


@inchi_comparison_routes.route("/api/compare_inchis_custom", methods=["POST"])
def compare_inchis_custom():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"message": "No JSON received"}), 400
        
        input1 = data.get("inchi1", "").strip()
        input2 = data.get("inchi2", "").strip()
        selected_layers = data.get("layers", [])
        
        if not input1 or not input2:
            return jsonify({"message": "Missing InChIs"}), 400
        
        inchi1 = InChI.normalize_input(input1)
        inchi2 = InChI.normalize_input(input2)
        
        base_config = load_config(_DEFAULT_CONFIG_PATH)
        config = build_config_from_layers(selected_layers, base_config)
        result = compare_pair(inchi1, inchi2, config)
        
        return jsonify({
            "inchi_1": inchi1,
            "inchi_2": inchi2,
            "results": result["results"]
        })
    
    except Exception as e:
        return jsonify({"message": str(e)}), 500


@inchi_comparison_routes.route("/api/inchi_layers", methods=["GET"])
def get_inchi_layers():
    return jsonify([
        {"key": "complete_identity",  "label": "Complete Identity"},
        {"key": "isotope",            "label": "Isotope Independence"},
        {"key": "salt",               "label": "Salt Independence"},
        {"key": "charge",             "label": "Charge Independence"},
        {"key": "double_bond",         "label": "Double Bond Independence"},
        {"key": "cis_trans",           "label": "Cis/Trans Independence"},
        {"key": "sn_position",         "label": "sn-Position Independence"},
        {"key": "chain_position",      "label": "Chain Position Independence"},
        {"key": "sum_composition",     "label": "Sum Composition Independence"},
        {"key": "tautomer",           "label": "Tautomer Independence"},
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
        config = load_config(_DEFAULT_CONFIG_PATH)
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
        layer = request.form.get("layer") or "COMPLETE_IDENTITY" 

        original_name1 = file1.filename or "file_a.mgf"
        original_name2 = file2.filename or "file_b.mgf"

        with tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".mgf") as tmp1:
            file1.save(tmp1.name)
            tmp1_path = tmp1.name

        with tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".mgf") as tmp2:
            file2.save(tmp2.name)
            tmp2_path = tmp2.name

        output_tmp = tempfile.mktemp(suffix=".mgf")

        config = load_config(_DEFAULT_CONFIG_PATH)
        result = compare_mgf_files(tmp1_path, tmp2_path, config, layer=layer, output_mgf=output_tmp)

        os.unlink(tmp1_path)
        os.unlink(tmp2_path)

        old_keys = list(result.get("input_counts", {}).keys())
        for old_key in old_keys:
            if tmp1_path.replace("\\", "/").split("/")[-1] in old_key.replace("\\", "/") or os.path.basename(tmp1_path) in old_key:
                result["input_counts"][original_name1] = result["input_counts"].pop(old_key)
                if "changes_breakdown" in result:
                    for bk in list(result["changes_breakdown"].keys()):
                        if os.path.basename(tmp1_path) in bk:
                            new_bk = bk.replace(os.path.basename(tmp1_path), original_name1)
                            result["changes_breakdown"][new_bk] = result["changes_breakdown"].pop(bk)
            elif tmp2_path.replace("\\", "/").split("/")[-1] in old_key.replace("\\", "/") or os.path.basename(tmp2_path) in old_key:
                result["input_counts"][original_name2] = result["input_counts"].pop(old_key)
                if "changes_breakdown" in result:
                    for bk in list(result["changes_breakdown"].keys()):
                        if os.path.basename(tmp2_path) in bk:
                            new_bk = bk.replace(os.path.basename(tmp2_path), original_name2)
                            result["changes_breakdown"][new_bk] = result["changes_breakdown"].pop(bk)

        if os.path.exists(output_tmp):
            with open(output_tmp, "rb") as f_out:
                result["output_mgf_b64"] = base64.b64encode(f_out.read()).decode("utf-8")
            os.unlink(output_tmp)

        return jsonify(result)

    except Exception as e:
        traceback.print_exc()
        return jsonify({"message": str(e)}), 500