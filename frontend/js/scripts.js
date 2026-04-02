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
    visualizeFromInchi(id1, inchi1);
    visualizeFromInchi(id2, inchi2);
}

async function visualizeFromInchi(containerId, inchi) {
    const element = document.getElementById(containerId);

    if (!element) {
        console.error("Missing container:", containerId);
        return;
    }

    element.innerHTML = "Loading...";

    try {
        const res = await fetch("http://127.0.0.1:5000/api/generate_3d", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ inchi })
        });

        const data = await res.json();

        if (!res.ok) throw new Error(data.error);

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
    const inchi1 = isAdvanced ? document.getElementById("inchi1_adv").value.trim() : document.getElementById("inchi1").value.trim();
    const inchi2 = isAdvanced ? document.getElementById("inchi2_adv").value.trim() : document.getElementById("inchi2").value.trim();

    if (!inchi1 || !inchi2) {
        showToast("Please enter both InChIs", "error");
        return;
    }

    if (isAdvanced) {
        const selectedLevels = document.querySelectorAll(".level-checkbox:checked");

        if (selectedLevels.length === 0) {
            showToast("Please select at least one identity level", "error");

            const panel = document.getElementById("layers-advanced");

            panel.classList.add("error-highlight");

            setTimeout(() => {
                panel.classList.remove("error-highlight");
            }, 800);

            return;
        }
    }

    draw(inchi1, inchi2);

    let url = "http://127.0.0.1:5000/api/compare_inchis";
    let body = { inchi1, inchi2 };

    if (isAdvanced) {
        url = "http://127.0.0.1:5000/api/compare_inchis_custom";
        body.levels = Array.from(document.querySelectorAll(".level-checkbox:checked")).map(cb => cb.value);
        console.log("SELECTED LEVELS:", body.levels);
    }

    try {
        const res = await fetch(url, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(body)
        });

        const text = await res.text();

        console.log("RAW RESPONSE:", text);

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
        if (!res.ok) { showToast(data.message || "Error", "error"); return; }

        updateLayers(mapBackendResults(data.results), isAdvanced);

    } catch (err) {
        console.error(err);
        showToast("Server error", "error");
    }
}

function mapBackendResults(results) {
    return {
        complete_identity: results.COMPLETE_IDENTITY,
        isotope: results.ISOTOPIC_INDEPENDENCE,
        salt: results.SALTS_INDEPENDENCE,
        charge: results.CHARGES_INDEPENDENCE,
        stereo_cis_trans: results.STEREOCHEMICAL_CIS_TRANS_INDEPENDENCE,
        double_bond: results.DOUBLE_BONDS_INDEPENDENCE,
        tautomer: results.TAUTOMER_INDEPENDENCE,
        substituent: results.SUBSTITUENT_POSITION_INDEPENDENCE
    };
}

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
            badge.innerText = "INDEPENDENT";
        } else if (match === false) {
            layer.classList.add("nomatch");
            layer.classList.remove("match");
            badge.className = "badge red";
            badge.innerText = "NOT INDEPENDENT";
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

const compareBtn = document.getElementById("compare-files-btn");
if (compareBtn) {
    compareBtn.addEventListener("click", async () => {
        if (!file1Data.length || !file2Data.length) {
            showToast("Upload both files", "error");
            return;
        }

        const res = await fetch("http://127.0.0.1:5000/api/compare_files", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
                list1: file1Data,
                list2: file2Data
            })
        });

        const data = await res.json();
        comparisonsData = data.comparisons;
        populateDropdown(comparisonsData);
    });
}

async function readFile(file) {
    const text = await file.text();
    return text.split("\n").map(l => l.trim()).filter(Boolean);
}

function populateDropdown(comparisons) {
    const container = document.getElementById("file-results-container");
    container.innerHTML = "";

    comparisons.forEach((comp, index) => {
        const wrapper = document.createElement("div");
        wrapper.className = "file-item";

        const header = document.createElement("div");
        header.className = "file-header";
        header.textContent = `${index + 1}. InChI 1 vs InChI 2`;

        const content = document.createElement("div");
        content.className = "file-content";
        content.style.display = "none";

        content.innerHTML = `
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

        header.onclick = () => {
            content.style.display = content.style.display === "none" ? "block" : "none";
            if (content.dataset.loaded) return;

            requestAnimationFrame(() => {
                draw(comp.inchi_1, comp.inchi_2, `mol1-${index}`, `mol2-${index}`);
            });

            const layersDiv = content.querySelector(`#layers-${index}`);
            layersDiv.innerHTML = "";
            layersDiv.className = "layers-grid";

            const mapped = mapBackendResults(comp.results);
            Object.entries(mapped).forEach(([key, val]) => {
                const div = document.createElement("div");
                div.className = "layer";

                div.innerHTML = `
                    <span>${key.replace(/_/g, " ")}</span>
                    <span class="badge ${val ? "green" : "red"}">
                        ${val ? "INDEPENDENT" : "NOT"}
                    </span>
                `;
                layersDiv.appendChild(div);
            });

            content.dataset.loaded = true;
        };

        wrapper.appendChild(header);
        wrapper.appendChild(content);
        container.appendChild(wrapper);
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