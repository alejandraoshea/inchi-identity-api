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

    if (Object.keys(inputCounts).length > 0) {
        var countSection = document.createElement("div");
        countSection.className = "mgf-file-counts";

        Object.keys(inputCounts).forEach(function(fname) {
            var row = document.createElement("div");
            row.className = "mgf-file-count";

            var breakdownKey = Object.keys(breakdown).find(function(k) {
                return k.indexOf(fname) === 0 && k.indexOf("internal") !== -1;
            });
            var internalChanges = breakdownKey !== undefined ? breakdown[breakdownKey] : null;

            row.innerHTML =
                "<span class='mgf-file-count-name'>" + fname + "</span>" +
                "<span class='mgf-file-count-nums'>" +
                    inputCounts[fname] + " entries" +
                    (internalChanges !== null ? " · <strong>" + internalChanges + "</strong> InChIs normalized internally" : "") +
                "</span>";
            countSection.appendChild(row);
        });

        var crossKey = Object.keys(breakdown).find(function(k) { return k.indexOf("_to_") !== -1; });
        if (crossKey && breakdown[crossKey] > 0) {
            var crossRow = document.createElement("div");
            crossRow.className = "mgf-file-count mgf-cross-row";
            crossRow.innerHTML =
                "<span class='mgf-file-count-name'>Cross-file</span>" +
                "<span class='mgf-file-count-nums'><strong>" + breakdown[crossKey] + "</strong> entries unified across files</span>";
            countSection.appendChild(crossRow);
        }

        container.appendChild(countSection);
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
            item.innerHTML =
                "<div class='mgf-log-arrow'>→</div>" +
                "<div class='mgf-log-detail'>" +
                    "<div class='mgf-log-inchi original'>" + change.original_inchi + "</div>" +
                    "<div class='mgf-log-into'>unified into</div>" +
                    "<div class='mgf-log-inchi canonical'>" + change.canonical_inchi + "</div>" +
                "</div>";
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