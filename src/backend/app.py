from flask import Flask, request, jsonify
from inchi.compare import compare_text_files, compare_pair
from inchi.determine_levels_id import InChI
from backend.inchi.config_loader import load_config
from backend.routes.inchi_comparison_routes import inchi_comparison_routes
from backend.routes.generation_3d_routes import generation_3d_routes
import tempfile
import json
from flask_cors import CORS

app = Flask(__name__)
CORS(app)
CORS(app, resources={r"/*": {"origins": "http://localhost:3000"}})

app.register_blueprint(inchi_comparison_routes)
app.register_blueprint(generation_3d_routes)

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080)