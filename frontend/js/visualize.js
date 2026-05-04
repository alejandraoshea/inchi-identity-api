var RDKit = null;
var rdkitReady = initRDKitModule().then(function(inst) { RDKit = inst; }).catch(function() {});
var sdfCache = {};
var fetching = {};

var modalBackdrop   = null;
var modalInchiLabel = null;
var modalCanvas     = null;
var modalViewer     = null;
var modalBuilt      = false;

function visualizeFromInchi(containerId, inchi) {
    var el = document.getElementById(containerId);
    if (!el || !inchi) 
        return;
    renderCard(el, inchi);
}

function draw(inchi1, inchi2, id1, id2) {
    visualizeFromInchi(id1 || "mol1", inchi1);
    visualizeFromInchi(id2 || "mol2", inchi2);
}

function drawPair(leftEl, rightEl, inchi1, inchi2) {
    if (leftEl  && inchi1) 
        renderCard(leftEl,  inchi1);
    if (rightEl && inchi2) 
        renderCard(rightEl, inchi2);
}

function renderCard(el, inchi) {
    el.innerHTML = "<div class='mol-loading'></div>";

    fetch("http://127.0.0.1:5000/api/render_3d_image", {
        method: "POST",
        headers: { 
            "Content-Type": "application/json" 
        },
        body: JSON.stringify({ 
            inchi: inchi 
        })
    })
    .then(function(r) { 
        return r.json(); 
    })
    .then(function(data) {
        if (data.error) 
            throw new Error(data.error);
        if (data.sdf) 
            sdfCache[inchi] = data.sdf;

        el.innerHTML = "";
        var img = document.createElement("img");
        img.src   = "data:image/png;base64," + data.image;
        img.style.cssText = "width:100%;height:100%;object-fit:contain;display:block;cursor:pointer;";
        img.title = "Click for interactive 3D";

        var hint = document.createElement("div");
        hint.className = "mol-hint";
        hint.textContent = "click for 3D \u2197";

        el.appendChild(img);
        el.appendChild(hint);

        img.addEventListener("click", function(e) {
            e.stopPropagation();
            openModal(inchi);
        });
    })
    .catch(function(err) {
        console.info("PNG failed, 2D fallback:", err.message);
        render2D(el, inchi);
    });
}

function render2D(el, inchi) {
    fetchSDF(inchi).then(function(sdf) {
        if (!sdf) { el.innerHTML = "<div class='mol-error'>Could not render</div>"; return; }
        rdkitReady.then(function() {
            try {
                var mol = RDKit.get_mol(sdf);
                if (!mol || !mol.is_valid()) throw new Error();
                var svg = mol.get_svg();
                mol.delete();
                el.innerHTML = svg;
                var s = el.querySelector("svg");
                if (s) {
                    s.removeAttribute("width");
                    s.removeAttribute("height");
                    s.style.cssText = "width:100%;height:100%;cursor:pointer;display:block;";
                }
                el.style.cursor = "pointer";
                el.title = "Click for interactive 3D";
                el.addEventListener("click", function(e) { e.stopPropagation(); openModal(inchi); });
            } catch(e) {
                el.innerHTML = "<div class='mol-error'>Could not render</div>";
            }
        });
    });
}

function buildModal() {
    if (modalBuilt) return;
    modalBuilt = true;

    var backdrop = document.createElement("div");
    backdrop.className = "mol-modal-backdrop";
    backdrop.style.display = "none";

    var box = document.createElement("div");
    box.className = "mol-modal";
    box.style.padding = "16px";

    var headerRow = document.createElement("div");
    headerRow.className = "mol-modal-header";

    var inchiLabel = document.createElement("div");
    inchiLabel.className = "mol-modal-header-text";

    var closeBtn = document.createElement("button");
    closeBtn.className = "mol-modal-close";
    closeBtn.innerHTML = "&#x2715;";

    headerRow.appendChild(inchiLabel);
    headerRow.appendChild(closeBtn);

    var canvasWrap = document.createElement("div");
    canvasWrap.className = "mol-modal-body";

    box.appendChild(headerRow);
    box.appendChild(canvasWrap);
    backdrop.appendChild(box);
    document.body.appendChild(backdrop);

    closeBtn.addEventListener("click", function() { 
        backdrop.style.display = "none"; 
    });

    backdrop.addEventListener("click", function(ev) { 
        if (ev.target === backdrop) 
            backdrop.style.display = "none"; 
        });
    document.addEventListener("keydown", function(ev) {
        if (ev.key === "Escape" && backdrop.style.display !== "none") backdrop.style.display = "none";
    });

    modalBackdrop   = backdrop;
    modalInchiLabel = inchiLabel;
    modalCanvas     = canvasWrap;

    modalViewer = $3Dmol.createViewer(canvasWrap, { backgroundColor: "white" });
}

function openModal(inchi) {
    buildModal();

    modalInchiLabel.textContent = inchi;
    modalBackdrop.style.display = "flex";

    var loadInto = function(sdf) {
        try { 
            modalViewer.clear(); 
        } catch(e) {}
        try { 
            modalViewer.removeAllModels(); 
        } catch(e) {}
        modalViewer.addModel(sdf, "sdf");
        modalViewer.setStyle({}, { 
            stick: {}, 
            sphere: { 
                scale: 0.3 
            } 
        });
        modalViewer.zoomTo();
        modalViewer.render();
        setTimeout(function() {
            try { 
                modalViewer.resize(); 
                modalViewer.render(); 
            } catch(e) {}
        }, 80);
    };

    if (sdfCache[inchi]) {
        loadInto(sdfCache[inchi]);
        return;
    }

    fetchSDF(inchi).then(function(sdf) {
        if (!sdf) { 
            return; 
        }
        loadInto(sdf);
    });
}

function fetchSDF(inchi) {
    if (!inchi) 
        return Promise.resolve(null);
    if (sdfCache[inchi]) 
        return Promise.resolve(sdfCache[inchi]);
    if (fetching[inchi]) 
        return fetching[inchi];

    var p = fetch("http://127.0.0.1:5000/api/generate_3d", {
        method: "POST",
        headers: { 
            "Content-Type": "application/json" 
        },
        body: JSON.stringify({ 
            inchi: inchi 
        })
    })
    .then(function(r) { 
        return r.json(); 
    })
    .then(function(data) {
        delete fetching[inchi];
        if (data.error) 
            throw new Error(data.error);
        sdfCache[inchi] = data.sdf;
        return data.sdf;
    })
    .catch(function() { 
        delete fetching[inchi]; 
        return null; 
    });

    fetching[inchi] = p;
    return p;
}