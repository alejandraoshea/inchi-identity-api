class InChiParser:
    def getMainLayer(inchi: str): #obtains chemical formula
        parts = inchi.split("/")
        for layer in parts:
            if layer.startswith("InChI="):
                continue
            if layer and layer[0].isalpha():
                return layer
        return None

    def getAtomConnectionsLayerSublayer(inchi: str):
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
