# 🧪 InChI Multi-Level Comparison Tool

A cheminformatics framework for **advanced comparison of chemical compounds** using hierarchical normalization of InChI representations.

This tool goes beyond exact string matching by evaluating molecular identity across multiple chemically meaningful levels, making it especially useful for **metabolomics, lipidomics, and chemical database integration**.

---

## Features

- Multi-level molecular comparison
- Intelligent normalization pipeline:
  - isotope removal
  - salt stripping
  - charge neutralization
  - stereochemistry abstraction
  - tautomer canonicalization
  - lipid-aware comparison logic (tail-based analysis)
- Command-Line Interface (CLI) for batch processing
- Web interface for interactive exploration (Flask)
- Supports `.txt` and `.mgf` metabolomics files 

---

## Comparison Levels

- Complete Identity
- Isotope Independence
- Salt Independence
- Charge Independence
- Stereochemistry Independence
- Double Bond Position Independence
- Tautomer Independence (via InChI Trust)
- Substituent Position Independence

---

## Installation
This project uses a Conda environment (recommended) with additional Python dependencies installed via pip.


### 1. Clone the repository

```bash
git clone https://github.com/alejandraoshea/identity-levels-inchi.git
cd identity-levels-inchi
```

### 2. Create and activate the conda environment

```bash
conda env create -f conda_env.yml
conda activate id-levels
```

### 3. Install Chromium

```bash
playwright install chromium
```

### 4. Install project (CLI)
```bash
pip install .
```

### 5. Verify CLI
```bash
inchi --help
```

### Requirements
1. Download InChI Trust (Required for Tautomer Level)
- Download from: https://www.inchi-trust.org/downloads/
Then either:
- Option A: Add to PATH:
```bash
    export PATH=$PATH:/path/to/inchi
```

- Option B: Pass manually in CLI
```bash
--inchitrust /path/to/inchi-1
```

### CLI Usage
1. Compare two InChIs directly
```bash
inchi compare-pair "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1" "InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h8-9H,2-7,10-17H2,1H3,(H,19,20)/b9-8+" > src/backend/output/result1.json
```

2. Compare two InChIs with selected levels
```bash
inchi compare-pair-levels "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1" "InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h8-9H,2-7,10-17H2,1H3,(H,19,20)/b9-8+" --levels isotope salt charge > src/backend/output/result2.json
```

3. Compare two text files (pairwise)
```bash
inchi compare file1.txt file2.txt > src/backend/output/result_pairwise.json
```

4. Cross comparison (all from file 1 with all from file 2)
```bash
inchi compare file1.txt file2.txt --mode cross > src/backend/output/result_cross.json
```

5. Show only true values when comparing .txt files
```bash
inchi compare file1.txt file2.txt --only-equal > src/backend/output/result_pairwise_true.json
inchi compare file1.txt file2.txt --mode cross --only-equal > src/backend/output/result_cross_true.json
```

6. Compare two mgf files
```bash

```

### Example Input
```bash
InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1
InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h8-9H,2-7,10-17H2,1H3,(H,19,20)/b9-8+
```

### Example Output
```bash
    {
        "inchi_1": "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1",
        "inchi_2": "InChI=1S/C18H34O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18(19)20/h8-9H,2-7,10-17H2,1H3,(H,19,20)/b9-8+",
        "results": {
            "COMPLETE_IDENTITY": false,
            "ISOTOPIC_INDEPENDENCE": false,
            "SALTS_INDEPENDENCE": false,
            "CHARGES_INDEPENDENCE": false,
            "DOUBLE_BONDS_INDEPENDENCE": false,
            "STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE": false,
            "TAUTOMER_INDEPENDENCE": true,
            "SUBSTITUENT_POSITION_INDEPENDENCE": false
        }
    }

```

### Author
Alejandra O'Shea Fernández