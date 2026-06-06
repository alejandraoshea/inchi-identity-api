function showToast(message, type) {
    type = type || "info";
    var container = document.getElementById("toast-container");
    if (!container) return;
    var toast = document.createElement("div");
    toast.classList.add("toast", "toast-" + type);
    toast.textContent = message;
    container.appendChild(toast);
    setTimeout(function() {
        toast.style.animation = "fadeOut 0.3s ease forwards";
        setTimeout(function() { toast.remove(); }, 300);
    }, 1800);
}

function setLoadingState(isLoading) {
    var btn = document.querySelector("button[data-compare]");
    if (btn) {
        btn.disabled    = isLoading;
        btn.textContent = isLoading ? "Comparing..." : "Compare";
    }
}

function autoResizeTextarea(el) {
    el.style.height = "auto";
    el.style.height = el.scrollHeight + "px";
}

function initTextareas() {
    document.querySelectorAll("textarea").forEach(function(ta) {
        ta.addEventListener("input", function() { autoResizeTextarea(ta); });
        autoResizeTextarea(ta);
    });
}

function markActiveNav() {
    var FILES_PAGES = ["files.html", "files-pairwise.html", "files-cross.html"];
    var page = location.pathname.split("/").pop();
    document.body.classList.toggle("allow-scroll", FILES_PAGES.indexOf(page) !== -1);
    document.querySelectorAll(".nav a, .dropdown-menu a").forEach(function(a) {
        a.classList.toggle("active", a.getAttribute("href") === page);
    });
}

function val(id) {
    var el = document.getElementById(id);
    return el ? el.value.trim() : "";
}