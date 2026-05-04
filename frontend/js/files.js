var fileAEntries = [];
var fileBEntries = [];

document.addEventListener("DOMContentLoaded", function() {

    var fileInputA = document.getElementById("file1");
    if (fileInputA) fileInputA.addEventListener("change", function(e) {
        var selectedFile = e.target.files[0]; 
        if (!selectedFile) return;

        readFile(selectedFile, function(lines) {
            fileAEntries = lines;

            var labelWrapper = fileInputA.closest("label");
            if (labelWrapper) { 
                labelWrapper.classList.add("has-file"); 
                var fileNameDisplay = labelWrapper.querySelector(".file-btn-name"); 
                if (fileNameDisplay) fileNameDisplay.textContent = selectedFile.name; 
            }

            showToast("Loaded " + lines.length + " entries from File A", "success");
        });
    });

    var fileInputB = document.getElementById("file2");
    if (fileInputB) fileInputB.addEventListener("change", function(e) {
        var selectedFile = e.target.files[0]; 
        if (!selectedFile) return;

        readFile(selectedFile, function(lines) {
            fileBEntries = lines;

            var labelWrapper = fileInputB.closest("label");
            if (labelWrapper) { 
                labelWrapper.classList.add("has-file"); 
                var fileNameDisplay = labelWrapper.querySelector(".file-btn-name"); 
                if (fileNameDisplay) fileNameDisplay.textContent = selectedFile.name; 
            }

            showToast("Loaded " + lines.length + " entries from File B", "success");
        });
    });

    var compareButton = document.getElementById("compare-files-btn");
    if (compareButton) compareButton.addEventListener("click", function() {

        if (!fileAEntries.length || !fileBEntries.length) { 
            showToast("Please upload both files first", "error"); 
            return; 
        }

        var comparisonMode = (document.getElementById("mode-select") || {}).value || "pairwise";

        compareButton.disabled = true; 
        compareButton.textContent = "Comparing...";

        fetch("http://127.0.0.1:5000/api/compare_files", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ 
                list1: fileAEntries, 
                list2: fileBEntries, 
                mode: comparisonMode 
            })
        })
        .then(function(response) { 
            if (!response.ok) return response.json().then(function(errorData) { 
                throw new Error(errorData.message || "Error"); 
            }); 
            return response.json(); 
        })
        .then(function(data) { 
            comparisonMode === "pairwise" 
                ? renderPairwise(data.comparisons) 
                : renderCross(data.comparisons); 
        })
        .catch(function(error) { 
            showToast(error.message || "Network error", "error"); 
        })
        .finally(function() { 
            compareButton.disabled = false; 
            compareButton.textContent = "Compare Files"; 
        });
    });
});

function renderPairwise(comparisons) {
    var resultsContainer = document.getElementById("file-results-container");
    resultsContainer.innerHTML = "";

    comparisons.forEach(function(comparison, index) {
        var card = buildCard(
            comparison, 
            index, 
            "pairwise", 
            mapResults(comparison.results || comparison.matches || {})
        );

        resultsContainer.appendChild(card);

        requestAnimationFrame(function() {
            drawPair(
                document.getElementById("left-"  + index),
                document.getElementById("right-" + index),
                comparison.inchi_1, 
                comparison.inchi_2
            );
        });
    });
}

function renderCross(comparisons) {
    var resultsContainer = document.getElementById("file-results-container");
    resultsContainer.innerHTML = "";

    var groupedByFirstMolecule = group(comparisons);

    Object.keys(groupedByFirstMolecule).forEach(function(inchi1, groupIndex) {
        var matchList = groupedByFirstMolecule[inchi1];

        var groupWrapper = document.createElement("div");
        groupWrapper.className = "file-item";

        var header = document.createElement("div");
        header.className = "file-header";

        header.innerHTML =
            "<span class='arrow'>&#9658;</span>" +
            "<span class='inchi-text'>" + inchi1 + "</span>" +
            "<span class='badge' style='margin-left:auto;flex-shrink:0'>" +
            matchList.length + " match" + (matchList.length !== 1 ? "es" : "") +
            "</span>";

        var content = document.createElement("div");
        content.style.cssText = "display:none;flex-direction:column;gap:12px;padding:10px;";

        var hasRendered = false;

        header.addEventListener("click", function() {
            var isOpen = groupWrapper.classList.toggle("open");
            content.style.display = isOpen ? "flex" : "none";
            header.querySelector(".arrow").style.transform = isOpen ? "rotate(90deg)" : "";

            if (isOpen && !hasRendered) {
                hasRendered = true;

                requestAnimationFrame(function() {
                    matchList.forEach(function(comparison, matchIndex) {
                        drawPair(
                            document.getElementById("left-"  + groupIndex + "-" + matchIndex),
                            document.getElementById("right-" + groupIndex + "-" + matchIndex),
                            comparison.inchi_1, 
                            comparison.inchi_2
                        );
                    });
                });
            }
        });

        matchList.forEach(function(comparison, matchIndex) {
            content.appendChild(
                buildCard(
                    comparison, 
                    groupIndex + "-" + matchIndex, 
                    "cross", 
                    mapResults(comparison.results || comparison.matches || {})
                )
            );
        });

        groupWrapper.appendChild(header);
        groupWrapper.appendChild(content);
        resultsContainer.appendChild(groupWrapper);
    });
}

function buildCard(comparison, index, mode, mappedResults) {
    var card = document.createElement("div");
    card.className = "comparison-card";

    card.innerHTML =
        "<div class='mol-card-headers'>" +
            "<div class='molecule-header inchi-text'>" + comparison.inchi_1 + "</div>" +
            "<div class='molecule-header inchi-text'>" + comparison.inchi_2 + "</div>" +
        "</div>" +
        "<div class='mol-card-row'>" +
            "<div id='left-"  + index + "' class='mol-cell'><div class='mol-loading'></div></div>" +
            "<div id='right-" + index + "' class='mol-cell'><div class='mol-loading'></div></div>" +
        "</div>" +
        "<div class='layers-grid' id='layers-" + index + "'></div>";

    var layersContainer = card.querySelector("#layers-" + index);

    Object.keys(layerLabels).forEach(function(layerKey) {
        var value = mappedResults[layerKey];

        if (mode === "pairwise" && value === null) return;
        if (mode === "cross" && value !== true) return;

        var layerItem = document.createElement("div");
        layerItem.className = "layer " + (value === true ? "match" : value === false ? "nomatch" : "");

        layerItem.innerHTML =
            "<span class='layer-label'>" + layerLabels[layerKey] + "</span>" +
            "<span class='badge " + (value === true ? "green" : value === false ? "red" : "") + "'>" +
            (value === true ? "EQUAL" : value === false ? "DIFF" : "N/A") +
            "</span>";

        layersContainer.appendChild(layerItem);
    });

    if (layersContainer.children.length === 0) {
        layersContainer.innerHTML =
            "<div class='layer' style='grid-column:1/-1;color:var(--muted);font-size:12px'>" +
            "No matching identity levels</div>";
    }

    return card;
}

function group(comparisons) {
    return comparisons.reduce(function(accumulator, comparison) { 
        (accumulator[comparison.inchi_1] = accumulator[comparison.inchi_1] || []).push(comparison); 
        return accumulator; 
    }, {});
}

function readFile(file, callback) {
    var fileReader = new FileReader();

    fileReader.onload = function(event) {
        callback(
            event.target.result
                .split("\n")
                .map(function(line) { return line.trim(); })
                .filter(Boolean)
        );
    };

    fileReader.readAsText(file);
}