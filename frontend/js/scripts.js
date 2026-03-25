let RDKit;

initRDKitModule().then(instance => {
    RDKit = instance;
});

async function draw() {
    const inchi1 = document.getElementById("inchi1").value.trim();
    const inchi2 = document.getElementById("inchi2").value.trim();

    visualizeFromInchi("mol1", inchi1);
    visualizeFromInchi("mol2", inchi2);
}

async function visualizeFromInchi(containerId, inchi) {
    const element = document.getElementById(containerId);
    element.innerHTML = "";

    try {
        const res = await fetch("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/SDF?record_type=3d", {
            method: "POST",
            headers: {
                "Content-Type": "application/x-www-form-urlencoded"
            },
            body: `inchi=${encodeURIComponent(inchi)}`
        });

        if (!res.ok) throw new Error();

        const data = await res.text();

        const viewer = $3Dmol.createViewer(element, { backgroundColor: 'white' });

        viewer.addModel(data, "sdf");
        viewer.setStyle({}, {
            stick: { radius: 0.15 },
            sphere: { scale: 0.25 }
        });

        viewer.zoomTo();
        viewer.render();

    } catch {
        try {
            const mol = RDKit.get_mol_from_inchi(inchi);
            element.innerHTML = mol.get_svg();
        } catch {
            element.innerHTML = "Invalid InChI";
        }
    }
}

function compare() {
    const inchi1 = document.getElementById("inchi1").value.trim();
    const inchi2 = document.getElementById("inchi2").value.trim();

    if (!inchi1 || !inchi2) {
        alert("Please enter both InChI strings.");
        return;
    }

    draw();

    const results = {
        complete: false,
        isotope: true,
        salt: true,
        charge: true,
        stereo: true,
        double: false,
        tautomer: true,
        substituent: false
    };

    updateLayers(results);
}

function updateLayers(results) {
    document.querySelectorAll(".layer").forEach(layer => {
        const key = layer.dataset.key;
        const match = results[key];

        layer.classList.remove("match", "nomatch");
        layer.innerHTML = key;

        const badge = document.createElement("span");
        badge.classList.add("badge");

        if (match) {
            layer.classList.add("match");
            badge.classList.add("green");
            badge.innerText = "INDEPENDENT";
        } else {
            layer.classList.add("nomatch");
            badge.classList.add("red");
            badge.innerText = "NOT INDEPENDENT";
        }

        layer.appendChild(badge);
    });
}

function switchMode(mode, el) {
    document.querySelectorAll(".tab").forEach(t => t.classList.remove("active"));
    el.classList.add("active");

    document.querySelectorAll(".tab-content").forEach(c => c.classList.remove("active"));
    document.getElementById("mode-" + mode).classList.add("active");
}