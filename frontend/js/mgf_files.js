var mgfFile1Data = null;
var mgfFile2Data = null;

document.addEventListener("DOMContentLoaded", function() {

    var f1 = document.getElementById("mgf-file1");
    if (f1) {
        f1.addEventListener("change", function(e) {
            var file = e.target.files[0];
            if (!file) return;
            mgfFile1Data = file;
            var lbl = f1.closest("label");
            if (lbl) {
                lbl.classList.add("has-file");
                var n = lbl.querySelector(".file-btn-name");
                if (n) n.textContent = file.name;
            }
            showToast("Loaded " + file.name, "success");
        });
    }

    var f2 = document.getElementById("mgf-file2");
    if (f2) {
        f2.addEventListener("change", function(e) {
            var file = e.target.files[0];
            if (!file) return;
            mgfFile2Data = file;
            var lbl = f2.closest("label");
            if (lbl) {
                lbl.classList.add("has-file");
                var n = lbl.querySelector(".file-btn-name");
                if (n) n.textContent = file.name;
            }
            showToast("Loaded " + file.name, "success");
        });
    }

    var btn = document.getElementById("mgf-compare-btn");
    if (btn) {
        btn.addEventListener("click", function() {
            if (!mgfFile1Data || !mgfFile2Data) {
                showToast("Please upload both MGF files first", "error");
                return;
            }

            var level = document.getElementById("mgf-level").value;

            btn.disabled    = true;
            btn.textContent = "Comparing...";

            var formData = new FormData();
            formData.append("file1", mgfFile1Data);
            formData.append("file2", mgfFile2Data);
            formData.append("level", level);

            fetch("http://127.0.0.1:5000/api/compare_mgf_files", {
                method: "POST",
                body: formData
            })
            .then(function(res) {
                if (!res.ok) return res.json().then(function(e) { throw new Error(e.message || "Server error"); });
                return res.json();
            })
            .then(function(data) {
                renderResults(data);
            })
            .catch(function(err) {
                showToast(err.message || "Network error", "error");
            })
            .finally(function() {
                btn.disabled    = false;
                btn.textContent = "Compare Files";
            });
        });
    }
});

function renderResults(data) {
    var container = document.getElementById("mgf-results-container");
    container.innerHTML = "";

    var inputCounts  = data.input_counts  || {};
    var outputCount  = data.output_count  || 0;
    var changesCount = data.changes_count || 0;
    var breakdown    = data.changes_breakdown || {};
    var changesLog   = data.changes_log   || [];

    var totalInput = Object.values(inputCounts).reduce(function(a, b) { return a + b; }, 0);

    var summary = document.createElement("div");
    summary.className = "mgf-summary";
    summary.innerHTML =
        "<div class='mgf-stat'>" +
            "<div class='mgf-stat-value'>" + totalInput + "</div>" +
            "<div class='mgf-stat-label'>Total Input Entries</div>" +
        "</div>" +
        "<div class='mgf-stat'>" +
            "<div class='mgf-stat-value'>" + outputCount + "</div>" +
            "<div class='mgf-stat-label'>Total Output Entries</div>" +
        "</div>" +
        "<div class='mgf-stat'>" +
            "<div class='mgf-stat-value'>" + changesCount + "</div>" +
            "<div class='mgf-stat-label'>InChIs Normalized</div>" +
        "</div>";
    container.appendChild(summary);

    if (data.output_mgf_b64) {
        var outputName = (document.getElementById("mgf-output-name").value.trim() || "unified_output") + ".mgf";
        var downloadBar = document.createElement("div");
        downloadBar.className = "mgf-download-bar";
        downloadBar.innerHTML =
            "<div class='mgf-download-info'>" +
                "<span class='mgf-download-icon'>⬇</span>" +
                "<div>" +
                    "<div style='font-size:13px;font-weight:600;color:var(--text);'>" + outputName + "</div>" +
                    "<div style='font-size:11px;color:var(--muted);'>Normalized MGF with unified InChI strings</div>" +
                "</div>" +
            "</div>" +
            "<button class='mgf-download-btn' id='mgf-dl-btn'>Download</button>";
        container.appendChild(downloadBar);

        document.getElementById("mgf-dl-btn").addEventListener("click", function() {
            var bytes = atob(data.output_mgf_b64);
            var arr   = new Uint8Array(bytes.length);
            for (var i = 0; i < bytes.length; i++) arr[i] = bytes.charCodeAt(i);
            var blob  = new Blob([arr], { type: "chemical/x-mdl-molfile" });
            var url   = URL.createObjectURL(blob);
            var a     = document.createElement("a");
            a.href     = url;
            a.download = outputName;
            a.click();
            URL.revokeObjectURL(url);
        });
    }

    if (changesLog.length > 0) {
        var logSection = document.createElement("div");
        logSection.className = "mgf-log-section";

        var logHeader = document.createElement("div");
        logHeader.className = "mgf-log-header";
        logHeader.innerHTML =
            "<span class='mgf-log-title'>Normalization Log</span>" +
            "<span class='mgf-log-count'>" + changesLog.length + " change" + (changesLog.length !== 1 ? "s" : "") + "</span>";
        logSection.appendChild(logHeader);

        var logList = document.createElement("div");
        logList.className = "mgf-log-list";

        changesLog.forEach(function(change) {
            var item = document.createElement("div");
            item.className = "mgf-log-item";

            var steps = "";

            var typeLabel = change.structure_type === "SMILES" ? "SMILES" : "InChI";
            steps +=
                "<div class='mgf-log-step'>" +
                    "<span class='mgf-log-step-label'>Original (" + typeLabel + ")</span>" +
                    "<span class='mgf-log-inchi original'>" + change.original_structure + "</span>" +
                "</div>";

            if (change.smiles_to_inchi) {
                steps +=
                    "<div class='mgf-log-step'>" +
                        "<span class='mgf-log-step-label'>→ Converted to InChI</span>" +
                        "<span class='mgf-log-inchi'>" + change.smiles_to_inchi + "</span>" +
                    "</div>";
            }

            if (change.normalized_inchi && change.normalized_inchi !== change.smiles_to_inchi) {
                steps +=
                    "<div class='mgf-log-step'>" +
                        "<span class='mgf-log-step-label'> → Normalized (" + change.structure_type + ")</span>" +
                        "<span class='mgf-log-inchi'>" + change.normalized_inchi + "</span>" +
                    "</div>";
            }

            steps +=
                "<div class='mgf-log-step'>" +
                    "<span class='mgf-log-step-label'> → Canonical InChI</span>" +
                    "<span class='mgf-log-inchi canonical'>" + change.canonical_inchi + "</span>" +
                "</div>";

            item.innerHTML =
                "<div class='mgf-log-arrow'>→</div>" +
                "<div class='mgf-log-detail'>" + steps + "</div>";

            logList.appendChild(item);
        });

        logSection.appendChild(logList);
        container.appendChild(logSection);
    } else {
        var noChanges = document.createElement("div");
        noChanges.className = "mgf-empty-state";
        noChanges.innerHTML =
            "<div class='mgf-empty-icon' style='font-size:20px'>✓</div>" +
            "<div class='mgf-empty-text'>No normalization needed — all InChIs are already equivalent at this identity level</div>";
        container.appendChild(noChanges);
    }
}