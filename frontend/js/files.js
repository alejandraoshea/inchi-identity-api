var file1Data = [];
var file2Data = [];

document.addEventListener("DOMContentLoaded", function() {

    var f1 = document.getElementById("file1");
    if (f1) f1.addEventListener("change", function(e) {
        var file = e.target.files[0]; if (!file) return;
        readFile(file, function(lines) {
            file1Data = lines;
            var lbl = f1.closest("label");
            if (lbl) { 
                lbl.classList.add("has-file"); 
                var n = lbl.querySelector(".file-btn-name"); 
                if (n) n.textContent = file.name; 
            }
            showToast("Loaded " + lines.length + " entries from File A", "success");
        });
    });

    var f2 = document.getElementById("file2");
    if (f2) f2.addEventListener("change", function(e) {
        var file = e.target.files[0]; if (!file) return;
        readFile(file, function(lines) {
            file2Data = lines;
            var lbl = f2.closest("label");
            if (lbl) { 
                lbl.classList.add("has-file"); 
                var n = lbl.querySelector(".file-btn-name"); 
                if (n) n.textContent = file.name; 
            }
            showToast("Loaded " + lines.length + " entries from File B", "success");
        });
    });

    var btn = document.getElementById("compare-files-btn");
    if (btn) btn.addEventListener("click", function() {
        if (!file1Data.length || !file2Data.length) { 
            showToast("Please upload both files first", "error"); 
            return; 
        }
        var mode = (document.getElementById("mode-select") || {}).value || "pairwise";
        btn.disabled = true; btn.textContent = "Comparing...";
        fetch("http://127.0.0.1:5000/api/compare_files", {
            method: "POST",
            headers: { 
                "Content-Type": "application/json" 
            },
            body: JSON.stringify({ 
                list1: file1Data, list2: file2Data, mode: mode 
            })
        })
        .then(function(r) { 
            if (!r.ok) return r.json().then(function(e) { 
                throw new Error(e.message || "Error"); }); 
                return r.json(); 
            })
        .then(function(data) { 
            mode === "pairwise" ? renderPairwise(data.comparisons) : renderCross(data.comparisons); 
        })
        .catch(function(err) { 
            showToast(err.message || "Network error", "error"); 
        })
        .finally(function() { 
            btn.disabled = false; btn.textContent = "Compare Files"; 
        });
    });
});

function renderPairwise(comparisons) {
    var container = document.getElementById("file-results-container");
    container.innerHTML = "";
    comparisons.forEach(function(comp, i) {
        var card = buildCard(comp, i, "pairwise", mapResults(comp.results || comp.matches || {}));
        container.appendChild(card);
        requestAnimationFrame(function() {
            drawPair(
                document.getElementById("left-"  + i),
                document.getElementById("right-" + i),
                comp.inchi_1, comp.inchi_2
            );
        });
    });
}

function renderCross(comparisons) {
    var container = document.getElementById("file-results-container");
    container.innerHTML = "";
    var grouped = group(comparisons);
    Object.keys(grouped).forEach(function(inchi1, gi) {
        var matches = grouped[inchi1];
        var item    = document.createElement("div"); item.className = "file-item";
        var header  = document.createElement("div"); header.className = "file-header";
        header.innerHTML = "<span class='arrow'>&#9658;</span><span class='inchi-text'>" + inchi1 + "</span>" +
            "<span class='badge' style='margin-left:auto;flex-shrink:0'>" + matches.length + " match" + (matches.length !== 1 ? "es" : "") + "</span>";
        var content = document.createElement("div");
        content.style.cssText = "display:none;flex-direction:column;gap:12px;padding:10px;";
        var rendered = false;
        header.addEventListener("click", function() {
            var open = item.classList.toggle("open");
            content.style.display = open ? "flex" : "none";
            header.querySelector(".arrow").style.transform = open ? "rotate(90deg)" : "";
            if (open && !rendered) {
                rendered = true;
                requestAnimationFrame(function() {
                    matches.forEach(function(comp, mi) {
                        drawPair(
                            document.getElementById("left-"  + gi + "-" + mi),
                            document.getElementById("right-" + gi + "-" + mi),
                            comp.inchi_1, comp.inchi_2
                        );
                    });
                });
            }
        });
        matches.forEach(function(comp, mi) {
            content.appendChild(buildCard(comp, gi + "-" + mi, "cross", mapResults(comp.results || comp.matches || {})));
        });
        item.appendChild(header); item.appendChild(content); container.appendChild(item);
    });
}

function buildCard(comp, idx, mode, mapped) {
    var card = document.createElement("div"); card.className = "comparison-card";
    card.innerHTML =
        "<div class='mol-card-headers'>" +
            "<div class='molecule-header inchi-text'>" + comp.inchi_1 + "</div>" +
            "<div class='molecule-header inchi-text'>" + comp.inchi_2 + "</div>" +
        "</div>" +
        "<div class='mol-card-row'>" +
            "<div id='left-"  + idx + "' class='mol-cell'><div class='mol-loading'></div></div>" +
            "<div id='right-" + idx + "' class='mol-cell'><div class='mol-loading'></div></div>" +
        "</div>" +
        "<div class='layers-grid' id='layers-" + idx + "'></div>";

    var ld = card.querySelector("#layers-" + idx);
    Object.keys(layerLabels).forEach(function(key) {
        var val = mapped[key];
        if (mode === "pairwise" && val === null) return;
        if (mode === "cross"    && val !== true)  return;
        var d = document.createElement("div");
        d.className = "layer " + (val === true ? "match" : val === false ? "nomatch" : "");
        d.innerHTML = "<span class='layer-label'>" + layerLabels[key] + "</span>" +
            "<span class='badge " + (val === true ? "green" : val === false ? "red" : "") + "'>" +
            (val === true ? "EQUAL" : val === false ? "DIFF" : "N/A") + "</span>";
        ld.appendChild(d);
    });
    if (ld.children.length === 0)
        ld.innerHTML = "<div class='layer' style='grid-column:1/-1;color:var(--muted);font-size:12px'>No matching identity levels</div>";
    return card;
}

function group(comparisons) {
    return comparisons.reduce(function(acc, c) { (acc[c.inchi_1] = acc[c.inchi_1] || []).push(c); 
        return acc; 
    }, 
    {});
}

function readFile(file, cb) {
    var r = new FileReader();
    r.onload = function(e) { cb(e.target.result.split("\n").map(function(l) { 
        return l.trim(); 
    }).filter(Boolean)); 
    };
    r.readAsText(file);
}