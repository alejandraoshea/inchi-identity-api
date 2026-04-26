let RDKit;

initRDKitModule().then(instance => {
    RDKit = instance;
});

function updateLayerVisibility(mode) {
    const container = document.getElementById("layers-container");
    if (mode === "advanced") {
        container.querySelectorAll("input.level-checkbox").forEach(cb => cb.disabled = false);
    } else {
        container.querySelectorAll("input.level-checkbox").forEach(cb => cb.disabled = true);
    }
}

function draw(inchi1, inchi2, id1 = "mol1", id2 = "mol2") {
    const el1 = document.getElementById(id1);
    const el2 = document.getElementById(id2);

    if (el1) el1.innerHTML = "<div class='mol-loading'>Loading 3D...</div>";
    if (el2) el2.innerHTML = "<div class='mol-loading'>Loading 3D...</div>";

    visualizeFromInchi(id1, inchi1);
    visualizeFromInchi(id2, inchi2);
}

async function visualizeFromInchi(containerId, inchi) {
    const element = document.getElementById(containerId);

    if (!element) {
        console.error("Missing container:", containerId);
        return;
    }

    try {
        const res = await fetch("http://127.0.0.1:5000/api/generate_3d", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ inchi })
        });

        const data = await res.json();
        if (!res.ok) throw new Error(data.error);

        element.innerHTML = "";

        const viewer = $3Dmol.createViewer(element, { backgroundColor: 'white' });
        viewer.addModel(data.sdf, "sdf");
        viewer.setStyle({}, { stick: {}, sphere: { scale: 0.3 } });
        viewer.zoomTo();
        viewer.render();

    } catch (e) {
        console.error("3D failed → fallback 2D:", e);

        try {
            const mol = RDKit.get_mol(inchi);
            element.innerHTML = mol ? mol.get_svg() : "Invalid InChI";
        } catch {
            element.innerHTML = "Invalid InChI";
        }
    }
}

async function compare(isAdvanced = false) {
    const inchi1 = isAdvanced 
        ? document.getElementById("inchi1_adv").value.trim() 
        : document.getElementById("inchi1").value.trim();

    const inchi2 = isAdvanced 
        ? document.getElementById("inchi2_adv").value.trim() 
        : document.getElementById("inchi2").value.trim();

    if (!inchi1 || !inchi2) {
        showToast("Please enter both InChIs", "error");
        return;
    }

    updateLayers({}, isAdvanced);
    setLoadingState(true);

    let url = "http://127.0.0.1:5000/api/compare_inchis";
    let body = { inchi1, inchi2 };

    if (isAdvanced) {
        url = "http://127.0.0.1:5000/api/compare_inchis_custom";
        body.levels = Array.from(document.querySelectorAll(".level-checkbox:checked"))
            .map(cb => cb.value);
    }

    try {
        const res = await fetch(url, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(body)
        });

        const text = await res.text();
        let data;

        try {
            data = JSON.parse(text);
        } catch {
            showToast("Server returned invalid JSON", "error");
            return;
        }

        if (!res.ok) {
            showToast(data.message || "Error", "error");
            return;
        }

        draw(inchi1, inchi2);
        updateLayers(mapBackendResults(data.results), isAdvanced);

    } catch (err) {
        console.error(err);
        showToast("Server error", "error");
    } finally {
        setLoadingState(false);
    }
}

function setLoadingState(isLoading) {
    const layers = document.querySelectorAll(".layer");
    const button = document.querySelector("button[onclick='compare()']");

    layers.forEach(layer => {
        let badge = layer.querySelector(".badge");

        if (!badge) {
            badge = document.createElement("span");
            badge.classList.add("badge");
            layer.appendChild(badge);
        }

        if (isLoading) {
            layer.classList.remove("match", "nomatch");

            badge.className = "badge loading";
            badge.innerText = "LOADING...";
        } else {
            badge.classList.remove("loading");
        }
    });

    if (button) {
        button.disabled = isLoading;
        button.innerText = isLoading ? "Comparing..." : "Compare";
    }
}

function mapBackendResults(results) {
    return {
        complete_identity: results.COMPLETE_IDENTITY ?? null,
        isotope: results.ISOTOPIC_INDEPENDENCE ?? null,
        salt: results.SALTS_INDEPENDENCE ?? null,
        charge: results.CHARGES_INDEPENDENCE ?? null,
        stereo_cis_trans: results.STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE ?? null,
        double_bond: results.DOUBLE_BONDS_INDEPENDENCE ?? null,
        tautomer: results.TAUTOMER_INDEPENDENCE ?? null,
        substituent: results.SUBSTITUENT_POSITION_INDEPENDENCE ?? null
    };
}

const layerLabels = {
    complete_identity: "Complete Identity",
    isotope: "Isotope Independence",
    salt: "Salt Independence",
    charge: "Charge Independence",
    double_bond: "Double Bond Independence",
    stereo_cis_trans: "Cis/Trans Independence",
    tautomer: "Tautomer Independence",
    substituent: "Substituent Independence"
};


function updateLayers(results, isAdvanced = false) {
    document.querySelectorAll(".layer").forEach(layer => {

        let key = null;
        let checkbox = layer.querySelector("input");

        if (checkbox) {
            key = checkbox.value;

            if (isAdvanced && !checkbox.checked) {
                layer.style.display = "none";
                return;
            } else {
                layer.style.display = "block";
            }

        } else {
            key = layer.dataset.key;
        }

        if (!key) return;

        const match = results[key];

        let badge = layer.querySelector(".badge");
        if (!badge) {
            badge = document.createElement("span");
            badge.classList.add("badge");
            layer.appendChild(badge);
        }

        if (match === true) {
            layer.classList.add("match");
            layer.classList.remove("nomatch");
            badge.className = "badge green";
            badge.innerText = "EQUAL";
        } else if (match === false) {
            layer.classList.add("nomatch");
            layer.classList.remove("match");
            badge.className = "badge red";
            badge.innerText = "DIFF";
        } else {
            layer.classList.remove("match", "nomatch");
            badge.className = "badge";
            badge.innerText = "N/A";
        }
    });
}

function showToast(message, type = "info") {
    const container = document.getElementById("toast-container");
    if (!container) return;

    const toast = document.createElement("div");
    toast.classList.add("toast", `toast-${type}`);
    toast.textContent = message;

    container.appendChild(toast);

    setTimeout(() => {
        toast.style.animation = "fadeOut 0.3s ease forwards";

        setTimeout(() => {
            toast.remove();
        }, 300);
    }, 1500);
}

let file1Data = [];
let file2Data = [];
let comparisonsData = [];

const file1Input = document.getElementById("file1");
if (file1Input) {
    file1Input.addEventListener("change", async (e) => {
        const file = e.target.files[0];
        file1Data = await readFile(file);
        e.target.parentElement.textContent = file.name;
    });
}

const file2Input = document.getElementById("file2");
if (file2Input) {
    file2Input.addEventListener("change", async (e) => {
        const file = e.target.files[0];
        file2Data = await readFile(file);
        e.target.parentElement.textContent = file.name;
    });
}

document.getElementById("compare-files-btn").addEventListener("click", async () => {

    const compareMode = document.getElementById("mode-select").value; 

    const res = await fetch("http://127.0.0.1:5000/api/compare_files", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
            list1: file1Data,
            list2: file2Data,
            mode: compareMode
        })
    });

    const data = await res.json();
    comparisonsData = data.comparisons;

    if (compareMode === "pairwise") {
        renderPairwise(comparisonsData);
    } else {
        populateDropdown(comparisonsData);
    }
});

async function readFile(file) {
    const text = await file.text();
    return text.split("\n").map(l => l.trim()).filter(Boolean);
}

function populateDropdown(comparisons) {
    const container = document.getElementById("file-results-container");
    container.innerHTML = "";

    const grouped = groupByInchi1(comparisons);

    Object.keys(grouped).forEach((inchi1, groupIndex) => {
        const item = document.createElement("div");
        item.className = "file-item";

        const header = document.createElement("div");
        header.className = "file-header";

        header.innerHTML = `
            <span class="arrow">▶</span>
            <span class="inchi-text">${inchi1}</span>
        `;

        const content = document.createElement("div");
        content.style.display = "none";
        content.style.padding = "10px";

        header.onclick = () => {
            const isOpen = item.classList.contains("open");

            item.classList.toggle("open");
            content.style.display = isOpen ? "none" : "block";

            const arrow = header.querySelector(".arrow");
            arrow.style.transform = isOpen ? "rotate(0deg)" : "rotate(90deg)";
        };

        grouped[inchi1].forEach((comp, index) => {

            const card = document.createElement("div");
            card.className = "comparison-card";

            card.innerHTML = `
                <div class="molecule-row">
                    <div class="molecule-card">
                        <div class="molecule-header inchi-text">${comp.inchi_1}</div>
                        <div id="mol1-${groupIndex}-${index}" class="molecule"></div>
                    </div>

                    <div class="molecule-card">
                        <div class="molecule-header inchi-text">${comp.inchi_2}</div>
                        <div id="mol2-${groupIndex}-${index}" class="molecule"></div>
                    </div>
                </div>

                <div id="layers-${groupIndex}-${index}" class="layers-grid"></div>
            `;

            content.appendChild(card);

            setTimeout(() => {
                draw(
                    comp.inchi_1,
                    comp.inchi_2,
                    `mol1-${groupIndex}-${index}`,
                    `mol2-${groupIndex}-${index}`
                );
            }, 0);

            const layersDiv = card.querySelector(`#layers-${groupIndex}-${index}`);
            const mapped = mapBackendResults(comp.results || comp.matches || {});

            Object.keys(layerLabels).forEach(key => {
                const val = mapped[key];
                if (val === null) return;

                const div = document.createElement("div");
                div.className = "layer";

                div.innerHTML = `
                    <span class="layer-label">${layerLabels[key]}</span>
                    <span class="badge green">EQUAL</span>
                `;

                layersDiv.appendChild(div);
            });

            if (layersDiv.children.length === 0) {
                layersDiv.innerHTML = `<div class="layer">No matching identity levels</div>`;
            }
        });

        item.appendChild(header);
        item.appendChild(content);
        container.appendChild(item);
    });
}

function renderPairwise(comparisons) {
    const container = document.getElementById("file-results-container");
    container.innerHTML = "";

    comparisons.forEach((comp, index) => {

        const card = document.createElement("div");
        card.className = "comparison-card";

        card.innerHTML = `
            <div class="molecule-row">
                <div class="molecule-card">
                    <div class="molecule-header inchi-text">${comp.inchi_1}</div>
                    <div id="mol1-${index}" class="molecule"></div>
                </div>

                <div class="molecule-card">
                    <div class="molecule-header inchi-text">${comp.inchi_2}</div>
                    <div id="mol2-${index}" class="molecule"></div>
                </div>
            </div>

            <div id="layers-${index}" class="layers-grid"></div>
        `;

        container.appendChild(card);

        draw(comp.inchi_1, comp.inchi_2, `mol1-${index}`, `mol2-${index}`);

        const layersDiv = card.querySelector(`#layers-${index}`);
        const mapped = mapBackendResults(comp.results || comp.matches || {});

        Object.keys(layerLabels).forEach(key => {
            const val = mapped[key];

            if (val !== true) return; 

            const div = document.createElement("div");
            div.className = "layer match";

            div.innerHTML = `
                <span class="layer-label">${layerLabels[key]}</span>
                <span class="badge green">EQUAL</span>
            `;

            layersDiv.appendChild(div);
        });
    });
}



function renderComparison(comp) {
    const { inchi_1, inchi_2, results } = comp;

    draw(inchi_1, inchi_2);
    updateLayers(mapBackendResults(results), false);
    console.log("Selected:", comp);
}

function autoResizeTextarea(el) {
    el.style.height = "auto";
    el.style.height = el.scrollHeight + "px";
}

document.querySelectorAll("textarea").forEach(textarea => {
    textarea.addEventListener("input", () => autoResizeTextarea(textarea));
    autoResizeTextarea(textarea);
});

function clearAdvancedSelection() {
    document.querySelectorAll(".level-checkbox").forEach(cb => {
        cb.checked = false;
    });

    document.querySelectorAll("#layers-advanced .layer").forEach(layer => {
        layer.style.display = "block";
        layer.classList.remove("match", "nomatch");

        const badge = layer.querySelector(".badge");
        if (badge) {
            badge.innerText = "";
            badge.className = "badge";
        }
    });

    showToast("Selection cleared", "info");
}

function groupByInchi1(comparisons) {
    const grouped = {};

    comparisons.forEach(comp => {
        if (!grouped[comp.inchi_1]) {
            grouped[comp.inchi_1] = [];
        }
        grouped[comp.inchi_1].push(comp);
    });

    return grouped;
}

