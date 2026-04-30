class MgfParser:

    def parse_mgf(file_path):
        entries = []
        current = {}

        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()

                if line == "BEGIN IONS":
                    current = {}

                elif line == "END IONS":
                    if current:
                        entries.append(current)

                elif "=" in line:
                    key, value = line.split("=", 1)
                    current[key.upper()] = value

        return entries

    def extract_inchis(entries):
        inchis = []

        for entry in entries:
            if "INCHI" in entry:
                inchis.append(entry["INCHI"])
            elif "SMILES" in entry:
                inchis.append(entry["SMILES"]) 

        return inchis
    
    def extract_structures(entries):
        result = []

        for entry in entries:
            structure = None

            if "INCHI" in entry:
                structure = entry["INCHI"]
            elif "SMILES" in entry:
                structure = entry["SMILES"]

            if structure:
                result.append({
                    "structure": structure,
                    "entry": entry
                })

        return result