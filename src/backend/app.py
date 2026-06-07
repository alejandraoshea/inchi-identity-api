from flask import Flask
from flask_cors import CORS
from backend.routes.inchi_comparison_routes import inchi_comparison_routes
from backend.routes.generation_3d_routes import generation_3d_routes

app = Flask(__name__)
CORS(app)

app.register_blueprint(inchi_comparison_routes)
app.register_blueprint(generation_3d_routes)

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080)