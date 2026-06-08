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

---

## Installation

RDKit must be installed via conda:

```bash
git clone https://github.com/alejandraoshea/inchi-identity-api.git
cd inchi-identity-api
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

```bash
flask --app src/api/app.py run --port 8080
```

The API will be available at `http://localhost:8080`.

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
development. For production, restrict allowed origins in `app.py`.

---

## Related repositories

- [inchi-identity](https://github.com/alejandraoshea/identity-levels-inchi) — Python comparison engine and CLI
- [inchi-identity-app](https://github.com/alejandraoshea/inchi-identity-app) — Web frontend
