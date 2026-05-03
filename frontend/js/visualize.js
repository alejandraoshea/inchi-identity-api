var RDKit = null;
var _rdkitReady = initRDKitModule().then(function (inst) {
    RDKit = inst;
    return inst;
}).catch(function (err) {
    console.warn("RDKit failed to load:", err);
});

var _viewers = [];
var _WEBGL_LIMIT = 14; 

var VISUALIZE_2D_ONLY = false;

function visualizeFromInchi(containerId, inchi) {
    var element = document.getElementById(containerId);
    if (!element) { console.warn("visualize: missing container", containerId); return; }
    if (!inchi)   { element.innerHTML = _errorHTML("No InChI"); return; }

    element.innerHTML = "<div class='mol-loading'>Loading...</div>";

    if (VISUALIZE_2D_ONLY) { _render2D(element, inchi); return; }

    if (element.offsetWidth === 0 || element.offsetHeight === 0) {
        _waitUntilVisible(element, function () { _render3DSingle(element, inchi); });
    } else {
        _render3DSingle(element, inchi);
    }
}

function drawGrid(gridContainerId, inchi1, inchi2) {
    var container = document.getElementById(gridContainerId);
    if (!container) { console.warn("drawGrid: missing container", gridContainerId); return function(){}; }

    container.innerHTML = "<div class='mol-loading'>Loading...</div>";

    if (VISUALIZE_2D_ONLY) {
        _render2DPair(container, inchi1, inchi2);
        return function(){};
    }

    if (container.offsetWidth === 0 || container.offsetHeight === 0) {
        var destroy = function(){};
        _waitUntilVisible(container, function () {
            destroy = _renderGridViewer(container, inchi1, inchi2);
        });
        return function() { destroy(); };
    }

    return _renderGridViewer(container, inchi1, inchi2);
}

function draw(inchi1, inchi2, id1, id2) {
    id1 = id1 || "mol1";
    id2 = id2 || "mol2";
    visualizeFromInchi(id1, inchi1);
    visualizeFromInchi(id2, inchi2);
}

function _renderGridViewer(container, inchi1, inchi2) {
    var p1 = _fetchSDF(inchi1);
    var p2 = _fetchSDF(inchi2);

    var destroyed = false;
    var viewer    = null;

    Promise.all([p1, p2]).then(function (sdfs) {
        if (destroyed) return;
        if (!sdfs[0] && !sdfs[1]) throw new Error("Both SDFs failed");

        if (container.offsetWidth === 0 || container.offsetHeight === 0) {
            throw new Error("Container has zero size after fetch");
        }

        while (_viewers.length >= _WEBGL_LIMIT) {
            _demoteOldest();
        }

        container.innerHTML = "";
        container.style.display = "flex";

        var config = {
            rows: 1, cols: 2,
            viewer_config: { backgroundColor: "white" }
        };

        var viewers = $3Dmol.createViewerGrid(container, config);
        var vLeft  = viewers[0][0];
        var vRight = viewers[0][1];

        if (sdfs[0]) {
            vLeft.addModel(sdfs[0], "sdf");
            vLeft.setStyle({}, { stick: {}, sphere: { scale: 0.3 } });
            vLeft.zoomTo();
            vLeft.render();
        } else {
            _render2D(_leftDiv(container), inchi1);
        }

        if (sdfs[1]) {
            vRight.addModel(sdfs[1], "sdf");
            vRight.setStyle({}, { stick: {}, sphere: { scale: 0.3 } });
            vRight.zoomTo();
            vRight.render();
        } else {
            _render2D(_rightDiv(container), inchi2);
        }

        viewer = { vLeft: vLeft, vRight: vRight };
        _viewers.push({ container: container, inchi1: inchi1, inchi2: inchi2, viewer: viewer });

    }).catch(function (err) {
        console.info("Grid 3D failed, falling back to 2D:", err.message);
        if (!destroyed) _render2DPair(container, inchi1, inchi2);
    });

    return function destroy() {
        destroyed = true;
        _viewers = _viewers.filter(function (v) { return v.container !== container; });
        _destroyContainer(container);
    };
}

function _render3DSingle(element, inchi) {
    if (element.offsetWidth === 0 || element.offsetHeight === 0) {
        return _render2D(element, inchi);
    }
    _fetchSDF(inchi).then(function (sdf) {
        if (!sdf) throw new Error("No SDF");
        if (element.offsetWidth === 0 || element.offsetHeight === 0) throw new Error("zero-size");

        while (_viewers.length >= _WEBGL_LIMIT) { _demoteOldest(); }

        element.innerHTML = "";
        var v = $3Dmol.createViewer(element, { backgroundColor: "white" });
        v.addModel(sdf, "sdf");
        v.setStyle({}, { stick: {}, sphere: { scale: 0.3 } });
        v.zoomTo();
        v.render();
        _viewers.push({ container: element, inchi1: inchi, inchi2: null, viewer: v });
    }).catch(function (err) {
        console.info("Single 3D failed:", err.message);
        _render2D(element, inchi);
    });
}

function _demoteOldest() {
    if (_viewers.length === 0) return;
    var oldest = _viewers.shift();
    _destroyContainer(oldest.container);
    _render2DPair(oldest.container, oldest.inchi1, oldest.inchi2);
}

function _destroyContainer(container) {
    try {
        var canvas = container.querySelector("canvas");
        if (canvas && canvas.parentNode) {
            canvas.parentNode.removeChild(canvas);
        }
        container.innerHTML = "";
    } catch(e) { /* ignore */ }
}

function _render2DPair(container, inchi1, inchi2) {
    container.innerHTML =
        "<div id='_2d_left_"  + container.id + "' style='flex:1;overflow:hidden;'></div>" +
        "<div id='_2d_right_" + container.id + "' style='flex:1;overflow:hidden;'></div>";
    container.style.display = "flex";
    container.style.gap     = "8px";

    if (inchi1) _render2D(document.getElementById("_2d_left_"  + container.id), inchi1);
    if (inchi2) _render2D(document.getElementById("_2d_right_" + container.id), inchi2);
}

function _render2D(element, inchi) {
    if (!element || !inchi) return;
    _fetchSDF(inchi).then(function (sdf) {
        if (!sdf) throw new Error("No SDF");
        return _rdkitReady.then(function () {
            if (!RDKit) throw new Error("RDKit not ready");
            var mol = RDKit.get_mol(sdf);
            if (!mol || !mol.is_valid()) throw new Error("Invalid molecule");
            var svg = mol.get_svg();
            mol.delete();
            element.innerHTML = svg;
            var svgEl = element.querySelector("svg");
            if (svgEl) {
                svgEl.removeAttribute("width");
                svgEl.removeAttribute("height");
                svgEl.style.width  = "100%";
                svgEl.style.height = "100%";
            }
            _attachClickModal(element, inchi);
        });
    }).catch(function (err) {
        console.error("_render2D error:", err.message);
        element.innerHTML = _errorHTML("Could not render");
    });
}

var _sdfCache = {};

function _fetchSDF(inchi) {
    if (_sdfCache[inchi]) return Promise.resolve(_sdfCache[inchi]);
    return fetch("http://127.0.0.1:5000/api/generate_3d", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ inchi: inchi })
    })
    .then(function (res) { return res.json(); })
    .then(function (data) {
        if (data.error) throw new Error(data.error);
        _sdfCache[inchi] = data.sdf; 
        return data.sdf;
    })
    .catch(function () { return null; });
}

function _leftDiv(container) {
    return container.querySelector("div:first-child") || container;
}
function _rightDiv(container) {
    return container.querySelector("div:last-child") || container;
}

function _waitUntilVisible(element, callback) {
    var observer = new IntersectionObserver(function (entries, obs) {
        entries.forEach(function (entry) {
            if (entry.isIntersecting) { obs.disconnect(); callback(); }
        });
    }, { threshold: 0.01 });
    observer.observe(element);
}

function _errorHTML(msg) {
    return "<div style='font-size:11px;color:#aaa;text-align:center;padding:8px'>" + msg + "</div>";
}

function _attachClickModal(element, inchi) {
    var svgEl = element.querySelector("svg");
    if (!svgEl) return;
    svgEl.style.cursor = "pointer";

    svgEl.addEventListener("click", function (e) {
        e.stopPropagation();

        var backdrop = document.createElement("div");
        backdrop.className = "mol-modal-backdrop";

        var modal = document.createElement("div");
        modal.className = "mol-modal";

        var closeBtn = document.createElement("button");
        closeBtn.className = "mol-modal-close";
        closeBtn.innerHTML = "&#x2715;";

        var headerEl = document.createElement("div");
        headerEl.className = "mol-modal-header";
        headerEl.textContent = inchi;

        var body = document.createElement("div");
        body.className = "mol-modal-body";
        var svgClone = svgEl.cloneNode(true);
        svgClone.removeAttribute("width");
        svgClone.removeAttribute("height");
        svgClone.style.width  = "100%";
        svgClone.style.height = "100%";
        svgClone.style.cursor = "default";
        body.appendChild(svgClone);

        modal.appendChild(closeBtn);
        modal.appendChild(headerEl);
        modal.appendChild(body);
        backdrop.appendChild(modal);
        document.body.appendChild(backdrop);

        function close() { backdrop.remove(); }
        closeBtn.addEventListener("click", close);
        backdrop.addEventListener("click", function (ev) {
            if (ev.target === backdrop) close();
        });
        document.addEventListener("keydown", function esc(ev) {
            if (ev.key === "Escape") { close(); document.removeEventListener("keydown", esc); }
        });
    });
}