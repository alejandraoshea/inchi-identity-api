import os

from flask import Flask, jsonify
from flask_cors import CORS
from werkzeug.middleware.proxy_fix import ProxyFix

from backend.routes.inchi_comparison_routes import inchi_comparison_routes
from backend.routes.generation_3d_routes import generation_3d_routes

app = Flask(__name__)

app.config["MAX_CONTENT_LENGTH"] = int(os.getenv("MAX_CONTENT_LENGTH", 25 * 1024 * 1024))
app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1, x_port=1)

cors_origins = os.getenv("CORS_ORIGINS")
CORS(app, resources={r"/api/*": {"origins": cors_origins.split(",") if cors_origins else "*"}})

app.register_blueprint(inchi_comparison_routes)
app.register_blueprint(generation_3d_routes)


@app.get("/api/health")
def health():
    return jsonify({"status": "ok"})


if __name__ == "__main__":
    app.run(
        host=os.getenv("FLASK_RUN_HOST", "127.0.0.1"),
        port=int(os.getenv("PORT", "8080")),
    )
