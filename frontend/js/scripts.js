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

async function compare(isAdvanced = false) {
    const inchi1 = document.getElementById("inchi1").value.trim();
    const inchi2 = document.getElementById("inchi2").value.trim();

    if (!inchi1 || !inchi2) {
        alert("Please enter both InChI strings.");
        return;
    }

    draw();

    let url = "/api/compare_inchis";
    let body = { inchi1, inchi2 };

    // ADVANCED MODE
    if (isAdvanced) {
        url = "/api/compare_inchis_custom";

        const selectedLevels = [];
        document.querySelectorAll(".layer input:checked").forEach(cb => {
            selectedLevels.push(cb.value);
        });

        body.levels = selectedLevels;
    }

    try {
        const res = await fetch(url, {
            method: "POST",
            headers: {
                "Content-Type": "application/json"
            },
            body: JSON.stringify(body)
        });

        const data = await res.json();

        if (!res.ok) {
            alert(data.message || "Error");
            return;
        }

        updateLayers(data.results);

    } catch (err) {
        console.error(err);
        alert("Server error");
    }
}


function updateLayers(results) {
    document.querySelectorAll(".layer").forEach(layer => {
        const checkbox = layer.querySelector("input");
        if (!checkbox) return;

        const key = checkbox.value; 
        const match = results[key];

        layer.classList.remove("match", "nomatch");

        let badge = layer.querySelector(".badge");
        if (!badge) {
            badge = document.createElement("span");
            badge.classList.add("badge");
            layer.appendChild(badge);
        }

        if (match) {
            layer.classList.add("match");
            badge.className = "badge green";
            badge.innerText = "MATCH";
        } else {
            layer.classList.add("nomatch");
            badge.className = "badge red";
            badge.innerText = "DIFFERENT";
        }
    });
}

document.querySelectorAll(".tab").forEach(tab => {
    tab.addEventListener("click", () => {
        const mode = tab.dataset.mode;

        document.querySelectorAll(".tab").forEach(t => t.classList.remove("active"));
        tab.classList.add("active");

        document.querySelectorAll(".tab-content").forEach(c => c.classList.remove("active"));
        document.getElementById("mode-" + mode).classList.add("active");
    });
});

function updateLayerMode(mode) {
    const checkboxes = document.querySelectorAll(".level-checkbox");

    checkboxes.forEach(cb => {
        if (mode === "advanced") {
            cb.style.display = "inline-block";
        } else {
            cb.style.display = "none";
        }
    });
}