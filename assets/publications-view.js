// Publications page: a two-view toggle.
//
// The HTML keeps the canonical "by type" listing (Book / Papers in journals /
// Refereed proceedings / etc.) untouched. Every <li> entry carries a
// data-theme="..." attribute. When the user switches to "By theme", we build
// (once, lazily) a parallel listing under #pub-by-theme by cloning each <li>
// into its theme bucket, then swap visibility between the two containers.
//
// Cloning (rather than moving) means anchors, titles, and styling all keep
// working under both views, and the data structure stays a single source of
// truth: edit a <li> once and both views update on next page load.

(function () {
    "use strict";

    const THEMES = [
        ["lattices",            "Lattices & Delaunay geometry"],
        ["polytopes",           "Polytopes"],
        ["tilings",             "Tilings"],
        ["plane-graphs",        "Plane graphs & polycycles"],
        ["cuts-metrics",        "Cuts, metrics & cones"],
        ["probability",         "Probability"],
        ["cohomology",          "Cohomology of arithmetic groups"],
        ["algebraic-geometry",  "Algebraic geometry & moduli"],
        ["math-physics",        "Mathematical physics"],
        ["geophysics",          "Geophysical modelling"],
    ];

    const byType = document.getElementById("pub-by-type");
    const byTheme = document.getElementById("pub-by-theme");
    let themeBuilt = false;

    function buildThemeView() {
        if (themeBuilt) return;
        themeBuilt = true;

        const buckets = {};
        byType.querySelectorAll("li[data-theme]").forEach((li) => {
            const t = li.dataset.theme;
            (buckets[t] || (buckets[t] = [])).push(li);
        });

        const frag = document.createDocumentFragment();
        for (const [key, label] of THEMES) {
            const items = buckets[key];
            if (!items || items.length === 0) continue;
            const section = document.createElement("section");
            const h3 = document.createElement("h3");
            h3.textContent = `${label} (${items.length})`;
            section.appendChild(h3);
            const ul = document.createElement("ul");
            for (const li of items) {
                ul.appendChild(li.cloneNode(true));
            }
            section.appendChild(ul);
            frag.appendChild(section);
        }
        byTheme.appendChild(frag);
    }

    document.querySelectorAll('input[name="pub-view"]').forEach((input) => {
        input.addEventListener("change", (e) => {
            if (e.target.value === "theme") {
                buildThemeView();
                byType.hidden = true;
                byTheme.hidden = false;
            } else {
                byType.hidden = false;
                byTheme.hidden = true;
            }
        });
    });
})();
