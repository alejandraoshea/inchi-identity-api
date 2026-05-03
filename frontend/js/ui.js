function showToast(message, type) {
    type = type || "info";
    var container = document.getElementById("toast-container");
    if (!container) return;

    var toast = document.createElement("div");
    toast.classList.add("toast", "toast-" + type);
    toast.textContent = message;
    container.appendChild(toast);

    setTimeout(function () {
        toast.style.animation = "fadeOut 0.3s ease forwards";
        setTimeout(function () { toast.remove(); }, 300);
    }, 1800);
}

function setLoadingState(isLoading) {
    var button = document.querySelector("button[data-compare]");
    if (button) {
        button.disabled    = isLoading;
        button.textContent = isLoading ? "Comparing…" : "Compare";
    }
}

function autoResizeTextarea(el) {
    el.style.height = "auto";
    el.style.height = el.scrollHeight + "px";
}

function initTextareas() {
    document.querySelectorAll("textarea").forEach(function (ta) {
        ta.addEventListener("input", function () { autoResizeTextarea(ta); });
        autoResizeTextarea(ta);
    });
}

function markActiveNav() {
    var page = location.pathname.split("/").pop();
    document.querySelectorAll(".nav a").forEach(function (a) {
        if (a.getAttribute("href") === page) a.classList.add("active");
        else a.classList.remove("active");
    });
}