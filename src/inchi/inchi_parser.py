from rdkit import Chem

class InChiParser:
    def getMainLayer(inchi: str): #obtains chemical formula
        parts = inchi.split("/")
        for layer in parts:
            if layer.startswith("InChI="):
                continue
            if layer and layer[0].isalpha():
                return layer
        return None

    def getAtomConnectionsSublayer(inchi: str):
        parts = inchi.split("/")
        for sublayer in parts:
            if sublayer.startswith("c"):
                return sublayer
        return None

    def getHydrogenAtomsSublayer(inchi: str):
        parts = inchi.split("/")
        for sublayer in parts:
            if sublayer.startswith("h"):
                return sublayer
        return None

    def getChargeSublayer(inchi: str):
        parts = inchi.split("/")
        for p in parts:
            if p.startswith("q"):
                return p
        return None

    def getProtonSublayer(inchi: str):
        parts = inchi.split("/")
        for p in parts:
            if p.startswith("p"):
                return p
        return None
    
    from rdkit import Chem

    def neutralize_molecule(mol):
        #neutralize common charged functional groups and returns a new molecule
        # (reactant_smarts, product_smiles)
        neutralization_reactions = [
            ('[n+;H]', 'n'),           # protonated pyridine
            ('[N+;!H0]', 'N'),         # ammonium
            ('[$([O-]);!$([O-][#7])]', 'O'),  # carboxylate, phosphate oxygens
            ('[S-]', 'S'),
            ('[N-]', 'N'),
        ]

        mol = Chem.RWMol(mol)

        for smarts, replacement in neutralization_reactions:
            patt = Chem.MolFromSmarts(smarts)
            while mol.HasSubstructMatch(patt):
                matches = mol.GetSubstructMatches(patt)
                for match in matches:
                    atom = mol.GetAtomWithIdx(match[0])
                    atom.SetFormalCharge(0)
                    atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)

        return mol.GetMol()

    def neutralize_inchi_and_compare(inchi: str):
        mol = Chem.MolFromInchi(inchi, sanitize=True)
        if mol is None:
            return None  # invalid InChI

        neutral_mol = InChiParser.neutralize_molecule(mol)

        neutral_inchi = Chem.MolToInchi(neutral_mol)

        return {
            "original": inchi,
            "neutral": neutral_inchi,
            "same": inchi == neutral_inchi
        }

    #stereochemical layer: double bonds and cumulenes (sublayer)
    def getDoubleBondsSublayer(inchi: str):
        parts = inchi.split("/")
        for p in parts:
            if p.startswith("b"):
                return p
        return None

    def getTetrahedralStereoSublayer(inchi: str):
        parts = inchi.split("/")
        for p in parts:
            if p.startswith("t") or p.startswith("m"):
                return p
        return None
    
    def getTypeStereoInfoSublayer(inchi: str):
        parts = inchi.split("/")
        for p in parts:
            if p.startswith("s"):
                return p
        return None
    
    def getIsotopicLayer(inchi: str):
        parts = inchi.split("/")
        for p in parts:
            if p.startswith("i"):
                return p
        return None

    def getIsotopicHydrogenSublayer(inchi: str):
        parts = inchi.split("/")

        for i, p in enumerate(parts):
            if p.startswith("i"):
                # CASE 1: hydrogen is in the next block: /i/h1H
                if i + 1 < len(parts) and parts[i + 1].startswith("h"):
                    return parts[i + 1]

                # CASE 2: hydrogen is inside the isotopic block: i2+2/h1H
                body = p[1:] 
                if "h" in body:
                    return "h" + body.split("h", 1)[1]

                return None
        return None

    def getIsotopicStereoSublayer(inchi: str):
        parts = inchi.split("/")

        idx = None #find isotopic index
        for i, p in enumerate(parts):
            if p.startswith("i"):
                idx = i
                break

        if idx is None:
            return None

        # we collect all sublayers starting with b, t, m, s after isotopic layer
        stereo_sublayers = []
        for p in parts[idx+1:]:
            if p and p[0] in ("b", "t", "m", "s"):
                stereo_sublayers.append(p)
            else:
                break

        if stereo_sublayers:
            return "/".join(stereo_sublayers)
        else:
            return None

    def getFixedHLayer(inchi: str):
        parts = inchi.split("/")
        for p in parts:
            if p.startswith("f"):
                return p
        return None

    def getReconnectedLayer(inchi: str):
        parts = inchi.split("/")
        for p in parts:
            if p.startswith("r"):
                return p
        return None

    def removeChargeLayersUsingParser(inchi: str) -> str:
        if inchi is None:
            return None

        parts = inchi.split("/")
        filtered_parts = []

        for p in parts:
            # remove p+N
            if p.startswith("p+"):
                continue

            # keep everything else for now: q+N 
            # and for p-N or q-N TODO
            filtered_parts.append(p)

        return "/".join(filtered_parts)
    
    def removeDoubleBondsSublayer(inchi: str) -> str:
        parts = inchi.split("/")
        filtered = [p for p in parts if not (p.startswith("b") or p.startswith("h"))]
        return "/".join(filtered)
    
    def getTautomerLayer(inchi: str):
        if inchi is None:
            return None

        main = InChiParser.getMainLayer(inchi)
        conn = InChiParser.getAtomConnectionsSublayer(inchi)

        if main is None or conn is None:
            return None

        # return a canonical tuple
        return (main, conn)
    
    def removeStereoLayersUsingParser(inchi: str) -> str:
        if inchi is None:
            return None

        db = InChiParser.getDoubleBondsSublayer(inchi)          # /b
        tet = InChiParser.getTetrahedralStereoSublayer(inchi)    # /t or /m
        typ = InChiParser.getTypeStereoInfoSublayer(inchi)       # /s

        parts = inchi.split("/")
        filtered_parts = []

        for p in parts:
            if p == db:
                continue
            if p == tet:
                continue
            if p == typ:
                continue
            filtered_parts.append(p)

        return "/".join(filtered_parts)

    def removeIsotopicLayersUsingParser(inchi: str) -> str:
        if inchi is None:
            return None

        isotopic_layer = InChiParser.getIsotopicLayer(inchi)
        isotopic_stereo = InChiParser.getIsotopicStereoSublayer(inchi)

        parts = inchi.split("/")
        filtered_parts = []

        for p in parts:
            # we skip isotopic layer and isotopic stereo sublayers
            if p == isotopic_layer:
                continue
            if isotopic_stereo:
                # isotopic_stereo might be joined with "/", split to individual parts to check
                stereo_parts = isotopic_stereo.split("/")
                if p in stereo_parts:
                    continue
            filtered_parts.append(p)

        return "/".join(filtered_parts)
    