from flask import Flask, request, jsonify
from inchi.compare import compare_files
from inchi.determine_levels_id import InChi
from inchi.config_loader import load_config
from routes.inchi_comparison_routes import inchi_comparison_routes
import tempfile
import json

app = Flask(__name__)

app.register_blueprint(inchi_comparison_routes)