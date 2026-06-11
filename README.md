# InChI Identity API

Flask-based REST API exposing the [InChI Identity](https://github.com/alejandraoshea/identity-levels-inchi) 
hierarchical molecular comparison framework for interactive and programmatic use.

Any system capable of making HTTP requests can submit molecular identifiers, retrieve 
equivalence profiles, and trigger MGF spectral library unification — independently of 
the web frontend.

---

## Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/inchi_layers` | List all available identity layers |
| POST | `/api/compare_inchis` | Full hierarchical comparison across all layers |
| POST | `/api/compare_inchis_custom` | Comparison restricted to selected layers |
| POST | `/api/compare_files` | Pairwise or cross-comparison of identifier files |
| POST | `/api/compare_mgf_files` | MGF spectral library unification |
| POST | `/api/generate_3d` | Generate 3D SDF structure from InChI or SMILES |
| GET | `/api/health` | Health check for reverse proxies and service monitors |

---

## Installation

### Option A — pip (recommended)

```bash
git clone https://github.com/alejandraoshea/inchi-identity-api.git
cd inchi-identity-api
pip install -r requirements.txt
```

### Option B — conda (alternative)

```bash
conda env create -f conda_env.yml
conda activate inchi-identity-api
pip install -e .
```

### InChI Trust (for Layer 6)

```bash
export INCHITRUST_PATH=/path/to/inchi-1
```

If not available, Layer 6 falls back to RDKit's TautomerEnumerator automatically.

---

## Running

### Local development

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -e .
PYTHONPATH=src python3 -m backend.app
```

The API will be available at `http://localhost:8080` by default.

Runtime settings can be configured with environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `PORT` | `8080` | Port used by the Flask development entry point |
| `FLASK_RUN_HOST` | `127.0.0.1` | Host used by the Flask development entry point |
| `MAX_CONTENT_LENGTH` | `26214400` | Maximum request body size in bytes |
| `CORS_ORIGINS` | `*` | Comma-separated list of allowed frontend origins |
| `INCHITRUST_PATH` | unset | Optional InChI Trust installation path |

### Production-style local run

Use Gunicorn instead of Flask's development server:

```bash
. .venv/bin/activate
gunicorn --workers 2 --timeout 120 --bind 127.0.0.1:8080 backend.app:app
```

Run this command from the repository root so `Naming_Example.xlsx` is available to
the `inchi-identity` dependency.

### Tests

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

---

## Deployment with Nginx and systemd

The API is intended to run behind Nginx on the same domain as the frontend:

- `/api/` is proxied to Gunicorn on `127.0.0.1:8080`.
- `/` serves the frontend `index.html`.
- Other frontend routes are handled by the frontend static/SPA fallback.
- HTTP requests are redirected to HTTPS.

Example templates are included in:

- `deploy/systemd/inchi-identity-api.service`
- `deploy/nginx/inchi-identity-api.conf`
- `deploy/inchi-identity-api.env.example`

1. Install the API on the server:

```bash
sudo mkdir -p /var/www/inchi-identity-api
sudo chown -R "$USER":"$USER" /var/www/inchi-identity-api
git clone https://github.com/alejandraoshea/inchi-identity-api.git /var/www/inchi-identity-api
cd /var/www/inchi-identity-api
python3 -m venv .venv
. .venv/bin/activate
pip install -e .
```

2. Configure the service environment:

```bash
sudo cp deploy/inchi-identity-api.env.example /etc/inchi-identity-api.env
sudo editor /etc/inchi-identity-api.env
```

Set `CORS_ORIGINS` to the deployed frontend origin, for example:

```bash
CORS_ORIGINS=https://metaboidentity.eps.uspceu.es
```

3. Install and start the systemd service:

```bash
sudo cp deploy/systemd/inchi-identity-api.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now inchi-identity-api
sudo systemctl status inchi-identity-api
```

4. Configure Nginx for the HTTPS domain:

```bash
sudo cp deploy/nginx/inchi-identity-api.conf /etc/nginx/sites-available/inchi-identity-api
sudo ln -s /etc/nginx/sites-available/inchi-identity-api /etc/nginx/sites-enabled/inchi-identity-api
sudo editor /etc/nginx/sites-available/inchi-identity-api
```

The included template uses `metaboidentity.eps.uspceu.es`, serves the frontend
from `/var/www/inchi-identity-app`, and expects Let's Encrypt certificates at
`/etc/letsencrypt/live/metaboidentity.eps.uspceu.es/`:

```nginx
listen 443 ssl http2;
server_name metaboidentity.eps.uspceu.es;
root /var/www/inchi-identity-app;
index index.html;

location /api/ {
    proxy_pass http://127.0.0.1:8080;
}

location / {
    try_files $uri $uri/ /index.html;
}
```

5. Reload Nginx:

```bash
sudo nginx -t
sudo systemctl reload nginx
```

6. Verify the deployment:

```bash
curl https://metaboidentity.eps.uspceu.es/api/health
curl https://metaboidentity.eps.uspceu.es/api/inchi_layers
```

The frontend remains available on the same domain outside `/api/`.

---

## Example requests

### Full comparison

```bash
curl -X POST http://localhost:8080/api/compare_inchis \
  -H "Content-Type: application/json" \
  -d '{"inchi1": "InChI=1S/C5H11NO2/...", "inchi2": "InChI=1S/C5H11NO2/..."}'
```

### Selected layers only

```bash
curl -X POST http://localhost:8080/api/compare_inchis_custom \
  -H "Content-Type: application/json" \
  -d '{"inchi1": "...", "inchi2": "...", "levels": ["isotope", "salt", "charge"]}'
```

### MGF unification

```bash
curl -X POST http://localhost:8080/api/compare_mgf_files \
  -F "file1=@file1.mgf" \
  -F "file2=@file2.mgf" \
  -F "layer=CHARGES_INDEPENDENCE"
```

---

## CORS

CORS is enabled by default to allow the frontend to connect from any origin during
development. For production, restrict allowed origins with `CORS_ORIGINS`.

---

## Related repositories

- [inchi-identity](https://github.com/alejandraoshea/identity-levels-inchi) — Python comparison engine and CLI
- [inchi-identity-app](https://github.com/alejandraoshea/inchi-identity-app) — Web frontend
