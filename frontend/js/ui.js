function showToast(message, type) {
    type = type || "info";
    var container = document.getElementById("toast-container");
    if (!container) 
        return;
    var toast = document.createElement("div");
    toast.classList.add("toast", "toast-" + type);
    toast.textContent = message;
    container.appendChild(toast);
    setTimeout(function() {
        toast.style.animation = "fadeOut 0.3s ease forwards";
        setTimeout(function() { 
            toast.remove(); 
        }, 300);
    }, 1800);
}

function setLoadingState(isLoading) {
    var btn = document.querySelector("button[data-compare]");
    if (btn) { 
        btn.disabled = isLoading; 
        btn.textContent = isLoading ? "Comparing..." : "Compare"; 
    }
}

function autoResizeTextarea(element) {
    element.style.height = "auto";
    element.style.height = element.scrollHeight + "px";
}

function initTextareas() {
    document.querySelectorAll("textarea").forEach(function(textarea) {
        textarea.addEventListener("input", function() { autoResizeTextarea(textarea); });
        autoResizeTextarea(textarea);
    });
}

function markActiveNav() {
    var page = location.pathname.split("/").pop();
    document.querySelectorAll(".nav a").forEach(function(a) {
        a.classList.toggle("active", a.getAttribute("href") === page);
    });
}