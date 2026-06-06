# InChI Identity

A hierarchical molecular identity comparison framework for untargeted metabolomics.

Instead of binary exact matching, this tool evaluates molecular equivalence across six 
progressive normalization layers, returning a structured equivalence profile that reflects 
the structural resolution actually supported by the experimental evidence. For lipids, 
a dedicated four-level hierarchy (Levels A–D) addresses cis/trans geometry, sn-position, 
intra-chain double bond position, and global sum composition.
---

## Comparison Layers

| Layer | Name | Description |
|-------|------|-------------|
| 1 | Complete Identity | Exact InChI string equality |
| 2 | Isotopic Independence | Equality after /i layer removal |
| 3 | Salt Independence | Equality after counterion removal |
| 4 | Charge Independence | Equality after charge/protonation normalization |
| 5 | Stereochemical/Isomeric Independence | Levels A–D (all molecules / lipids only) |
| 6 | Tautomeric Independence | Equality after canonical tautomer generation |

---

## Installation

RDKit must be installed via conda:

```bash
git clone https://github.com/alejandraoshea/identity-levels-inchi.git
cd identity-levels-inchi
conda env create -f conda_env.yml
conda activate inchi-identity
pip install -e .
```

Verify:

```bash
inchi --help
```

### InChI Trust (optional, for Layer 6)

Download from https://www.inchi-trust.org/downloads/ then set the path:

```bash
export INCHITRUST_PATH=/path/to/inchi-1
```

If not available, Layer 6 falls back to RDKit's TautomerEnumerator automatically.

---

## CLI Usage

### Compare two structures across all layers

```bash
inchi compare-pair "<inchi_1>" "<inchi_2>"
```

### Compare with selected layers only

```bash
inchi compare-pair-layers "<inchi_1>" "<inchi_2>" --layers isotope salt charge
```

### File-based comparison
# Pairwise
```bash
inchi compare file1.txt file2.txt
```
# Cross-comparison (all vs all)
```bash
inchi compare file1.txt file2.txt --mode cross
```

# Only show equivalent pairs
```bash
inchi compare file1.txt file2.txt --only-equal
inchi compare file1.txt file2.txt --mode cross --only-equal
```

### MGF spectral library unification
```bash
inchi compare-mgf file1.mgf file2.mgf \
    --layer CHARGES_INDEPENDENCE \
    --output-mgf unified.mgf \
    --output-log unified_log.json
```

Available layers: `COMPLETE_IDENTITY`, `ISOTOPIC_INDEPENDENCE`, `SALTS_INDEPENDENCE`,
`CHARGES_INDEPENDENCE`, `DOUBLE_BONDS_INDEPENDENCE`,
`STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE`, `TAUTOMER_INDEPENDENCE`

---

## Example output

```json
{
    "inchi_1": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1",
    "inchi_2": "InChI=1S/C18H34O2/...",
    "results": {
        "COMPLETE_IDENTITY": false,
        "ISOTOPIC_INDEPENDENCE": false,
        "SALTS_INDEPENDENCE": false,
        "CHARGES_INDEPENDENCE": false,
        "DOUBLE_BONDS_INDEPENDENCE": false,
        "STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE": false,
        "TAUTOMER_INDEPENDENCE": false
    }
}
```

---

## Related repositories

- [inchi-identity-api](https://github.com/alejandraoshea/inchi-identity-api) — Flask REST backend
- [inchi-identity-app](https://github.com/alejandraoshea/inchi-identity-app) — Web frontend

---
