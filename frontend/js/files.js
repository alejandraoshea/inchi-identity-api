var file1Data = [];
var file2Data = [];

document.addEventListener("DOMContentLoaded", function () {

    var file1Input = document.getElementById("file1");
    if (file1Input) {
        file1Input.addEventListener("change", function (e) {
            var file = e.target.files[0];
            if (!file) return;
            readFile(file, function (lines) {
                file1Data = lines;
                var label = file1Input.closest("label");
                if (label) {
                    label.classList.add("has-file");
                    var nameEl = label.querySelector(".file-btn-name");
                    if (nameEl) nameEl.textContent = file.name;
                }
                showToast("Loaded " + lines.length + " entries from File A", "success");
            });
        });
    }

    var file2Input = document.getElementById("file2");
    if (file2Input) {
        file2Input.addEventListener("change", function (e) {
            var file = e.target.files[0];
            if (!file) return;
            readFile(file, function (lines) {
                file2Data = lines;
                var label = file2Input.closest("label");
                if (label) {
                    label.classList.add("has-file");
                    var nameEl = label.querySelector(".file-btn-name");
                    if (nameEl) nameEl.textContent = file.name;
                }
                showToast("Loaded " + lines.length + " entries from File B", "success");
            });
        });
    }

    var compareBtn = document.getElementById("compare-files-btn");
    if (compareBtn) {
        compareBtn.addEventListener("click", function () {
            if (!file1Data.length || !file2Data.length) {
                showToast("Please upload both files first", "error");
                return;
            }
            var mode = document.getElementById("mode-select")
                ? document.getElementById("mode-select").value
                : "pairwise";

            compareBtn.disabled    = true;
            compareBtn.textContent = "Comparing...";

            fetch("http://127.0.0.1:5000/api/compare_files", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({ list1: file1Data, list2: file2Data, mode: mode })
            })
            .then(function (res) {
                if (!res.ok) return res.json().then(function (e) { throw new Error(e.message || "Server error"); });
                return res.json();
            })
            .then(function (data) {
                if (mode === "pairwise") renderPairwise(data.comparisons);
                else                     renderCross(data.comparisons);
            })
            .catch(function (err) {
                console.error(err);
                showToast(err.message || "Network error", "error");
            })
            .finally(function () {
                compareBtn.disabled    = false;
                compareBtn.textContent = "Compare Files";
            });
        });
    }

});

function renderPairwise(comparisons) {
    var container = document.getElementById("file-results-container");
    container.innerHTML = "";

    comparisons.forEach(function (comp, index) {
        var mapped = mapResults(comp.results || comp.matches || {});
        var card   = _buildCard(comp, index, "pairwise", mapped);
        container.appendChild(card);

        requestAnimationFrame(function () {
            requestAnimationFrame(function () {
                drawGrid("grid-" + index, comp.inchi_1, comp.inchi_2);
            });
        });
    });
}

function renderCross(comparisons) {
    var container = document.getElementById("file-results-container");
    container.innerHTML = "";

    var grouped = _groupByInchi1(comparisons);

    Object.keys(grouped).forEach(function (inchi1, groupIndex) {
        var item = document.createElement("div");
        item.className = "file-item";

        var matches = grouped[inchi1];

        var header = document.createElement("div");
        header.className = "file-header";
        header.innerHTML =
            "<span class='arrow'>&#9658;</span>" +
            "<span class='inchi-text'>" + inchi1 + "</span>" +
            "<span class='badge' style='margin-left:auto;flex-shrink:0'>" +
                matches.length + " match" + (matches.length !== 1 ? "es" : "") +
            "</span>";

        var content = document.createElement("div");
        content.style.display    = "none";
        content.style.padding    = "10px";
        content.style.flexDirection = "column";
        content.style.gap        = "12px";

        var destroyFns = [];

        header.addEventListener("click", function () {
            var open = item.classList.toggle("open");
            header.querySelector(".arrow").style.transform = open ? "rotate(90deg)" : "rotate(0deg)";

            if (open) {
                content.style.display = "flex";
                requestAnimationFrame(function () {
                    requestAnimationFrame(function () {
                        matches.forEach(function (comp, index) {
                            var gridId = "grid-" + groupIndex + "-" + index;
                            var gridEl = document.getElementById(gridId);
                            if (gridEl) {
                                gridEl.innerHTML = "<div class='mol-loading'>Loading...</div>";
                                var destroy = drawGrid(gridId, comp.inchi_1, comp.inchi_2);
                                destroyFns[index] = destroy;
                            }
                        });
                    });
                });
            } else {
                content.style.display = "none";
                setTimeout(function () {
                    destroyFns.forEach(function (fn) { if (fn) fn(); });
                    destroyFns = [];
                }, 50);
            }
        });

        matches.forEach(function (comp, index) {
            var mapped = mapResults(comp.results || comp.matches || {});
            var card   = _buildCard(comp, groupIndex + "-" + index, "cross", mapped);
            content.appendChild(card);
        });

        item.appendChild(header);
        item.appendChild(content);
        container.appendChild(item);
    });
}

function _buildCard(comp, index, mode, mapped) {
    var card = document.createElement("div");
    card.className = "comparison-card";

    card.innerHTML =
        "<div class='mol-grid-headers'>" +
            "<div class='molecule-header inchi-text'>" + comp.inchi_1 + "</div>" +
            "<div class='molecule-header inchi-text'>" + comp.inchi_2 + "</div>" +
        "</div>" +
        "<div id='grid-" + index + "' class='mol-grid-container'></div>" +
        "<div class='layers-grid' id='layers-" + index + "'></div>";

    var layersDiv = card.querySelector("#layers-" + index);

    Object.keys(layerLabels).forEach(function (key) {
        var val = mapped[key];
        if (mode === "pairwise" && val === null) return;
        if (mode === "cross"    && val !== true) return;

        var div = document.createElement("div");
        var stateClass = val === true ? "match" : val === false ? "nomatch" : "";
        var badgeClass = val === true ? "green"  : val === false ? "red"    : "";
        var badgeText  = val === true ? "EQUAL"  : val === false ? "DIFF"   : "N/A";

        div.className = "layer " + stateClass;
        div.innerHTML =
            "<span class='layer-label'>" + layerLabels[key] + "</span>" +
            "<span class='badge " + badgeClass + "'>" + badgeText + "</span>";
        layersDiv.appendChild(div);
    });

    if (layersDiv.children.length === 0) {
        layersDiv.innerHTML =
            "<div class='layer' style='grid-column:1/-1;color:var(--muted);font-size:12px'>No matching identity levels</div>";
    }

    return card;
}

function _groupByInchi1(comparisons) {
    return comparisons.reduce(function (acc, comp) {
        if (!acc[comp.inchi_1]) acc[comp.inchi_1] = [];
        acc[comp.inchi_1].push(comp);
        return acc;
    }, {});
}

function readFile(file, callback) {
    var reader = new FileReader();
    reader.onload = function (e) {
        var lines = e.target.result
            .split("\n")
            .map(function (l) { return l.trim(); })
            .filter(Boolean);
        callback(lines);
    };
    reader.readAsText(file);
}