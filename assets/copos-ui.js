// Tools page glue: for each <section class="copos-tool" data-mode="..."> on
// the page, build an n×n input grid, call the corresponding WASM binding,
// and render the result. Each section is independent — its own dimension,
// matrix, presets, status, and result.
//
// Depends on /wasm/copos.js (Emscripten module factory) being loaded first.

(function () {
    "use strict";

    const MIN_N = 2;
    const MAX_N = 8;

    const PRESETS = {
        "identity": {
            label: "Identity I_n (copositive, completely positive)",
            fill: (n) => Array.from({ length: n * n }, (_, k) =>
                Math.floor(k / n) === (k % n) ? "1" : "0"),
        },
        "neg-diag": {
            label: "diag(-1, 1, ..., 1) (NOT copositive)",
            fill: (n) => Array.from({ length: n * n }, (_, k) => {
                const i = Math.floor(k / n), j = k % n;
                if (i !== j) return "0";
                return i === 0 ? "-1" : "1";
            }),
        },
        "horn": {
            label: "Horn matrix H_5 (copositive, NOT completely positive)",
            fill: () => [
                "1", "-1", "1", "1", "-1",
                "-1", "1", "-1", "1", "1",
                "1", "-1", "1", "-1", "1",
                "1", "1", "-1", "1", "-1",
                "-1", "1", "1", "-1", "1",
            ],
            n: 5,
        },
    };

    // Each entry describes one section's WASM call and its result phrasing.
    // `renderer` takes (result, ctx) and returns an HTML string for the
    // result pane. `ctx` exposes the helpers it might need (escapeHtml).
    const MODES = {
        copositivity: {
            method: "testCopositivity",
            renderer: (r, ctx) => renderCopositivity(r, "copositive",
                "v ≥ 0 with v<sup>T</sup>Av &lt; 0", ctx),
        },
        factorization: {
            method: "testCopositiveFactorization",
            renderer: (r, ctx) => renderFactorization(r, ctx),
        },
    };

    function renderCopositivity(r, label, witnessCondition, ctx) {
        const out = [];
        if (r.isCopositive) {
            out.push(`<p class="verdict verdict-yes">Matrix is ${label}.</p>`);
            out.push(`<p>Reason: <code>${ctx.escapeHtml(r.nature)}</code></p>`);
        } else {
            out.push(`<p class="verdict verdict-no">Matrix is NOT ${label}.</p>`);
            out.push(`<p>Nature of violation: <code>${ctx.escapeHtml(r.nature)}</code></p>`);
            if (r.witness && r.witness.size && r.witness.size() > 0) {
                const v = [];
                for (let i = 0; i < r.witness.size(); i++) v.push(r.witness.get(i));
                out.push(`<p>Witness vector ${witnessCondition}:</p>`);
                out.push(`<pre>v = (${v.map(ctx.escapeHtml).join(", ")})</pre>`);
            }
        }
        return out.join("");
    }

    function readVecFromFlat(flat, idx, n) {
        const v = [];
        for (let j = 0; j < n; j++) v.push(flat.get(idx * n + j));
        return v;
    }

    function renderFactorization(r, ctx) {
        const out = [];
        const n = r.dim;
        if (r.hasFactorization) {
            const nb = r.nBlocks;
            out.push(`<p class="verdict verdict-yes">Matrix admits a copositive factorization.</p>`);
            out.push(`<p>A = ${arraySum(nb, "α", "v")} with:</p>`);
            const rows = [];
            for (let i = 0; i < nb; i++) {
                const alpha = r.coefficients.get(i);
                const v = readVecFromFlat(r.vectorsFlat, i, n);
                rows.push(
                    `<tr><td>α<sub>${i + 1}</sub> = ${ctx.escapeHtml(alpha)}</td>` +
                    `<td>v<sub>${i + 1}</sub> = (${v.map(ctx.escapeHtml).join(", ")})</td></tr>`
                );
            }
            out.push(`<table class="copos-factor-table">${rows.join("")}</table>`);
        } else {
            out.push(`<p class="verdict verdict-no">Matrix does NOT admit a copositive factorization.</p>`);
            if (r.certificateFlat && r.certificateFlat.size && r.certificateFlat.size() === n * n) {
                out.push(`<p>Certificate: there is a copositive matrix C with ⟨C, A⟩ &lt; 0:</p>`);
                const lines = [];
                for (let i = 0; i < n; i++) {
                    const row = [];
                    for (let j = 0; j < n; j++) row.push(r.certificateFlat.get(i * n + j));
                    lines.push("  " + row.map(ctx.escapeHtml).join("  "));
                }
                out.push(`<pre>C =\n${lines.join("\n")}</pre>`);
            }
        }
        return out.join("");
    }

    function arraySum(nb, coefName, vecName) {
        const parts = [];
        const cap = Math.min(nb, 3);
        for (let i = 1; i <= cap; i++) {
            parts.push(`${coefName}<sub>${i}</sub> ${vecName}<sub>${i}</sub> ${vecName}<sub>${i}</sub><sup>T</sup>`);
        }
        if (nb > cap) parts.push("...");
        return parts.join(" + ");
    }

    // Singleton — the WASM module is shared between all sections on the page.
    let modulePromise = null;
    function loadModule() {
        if (modulePromise) return modulePromise;
        if (typeof createCoposModule !== "function") {
            return Promise.reject(new Error(
                "copos.js (Emscripten module) failed to load"));
        }
        // The cache-buster comes from window.__WASM_BUST, which is set inline
        // in Tools.html from Jekyll's site.time. Each site rebuild gets a new
        // value, so browsers re-fetch the wasm after we redeploy. Without
        // this, hard-refresh does not bust the wasm cache and embind init
        // can fail when a new copos.js meets a stale copos.wasm.
        const bust = window.__WASM_BUST ? "?v=" + window.__WASM_BUST : "";
        modulePromise = createCoposModule({
            locateFile: (path) => path.endsWith(".wasm") ? "/wasm/" + path + bust : path,
        });
        return modulePromise;
    }

    function escapeHtml(s) {
        return String(s).replace(/[&<>"']/g, (c) => ({
            "&": "&amp;", "<": "&lt;", ">": "&gt;",
            "\"": "&quot;", "'": "&#39;",
        }[c]));
    }

    function formatElapsed(ms) {
        if (ms < 1) return ms.toFixed(2) + " ms";
        if (ms < 1000) return ms.toFixed(1) + " ms";
        return (ms / 1000).toFixed(2) + " s";
    }

    // Parse the "nRows nCols\n<entries>" file format used by polyhedral_common.
    function parseMatrixFile(text) {
        const toks = text.split(/\s+/).filter(Boolean);
        if (toks.length < 2) throw new Error("File must start with: nRows nCols");
        const nRows = parseInt(toks[0], 10);
        const nCols = parseInt(toks[1], 10);
        if (!Number.isInteger(nRows) || !Number.isInteger(nCols) || nRows <= 0 || nCols <= 0) {
            throw new Error("Header (nRows nCols) must be two positive integers");
        }
        if (nRows !== nCols) {
            throw new Error(`Matrix must be square, got ${nRows} × ${nCols}`);
        }
        if (nRows < MIN_N || nRows > MAX_N) {
            throw new Error(`Dimension ${nRows} is out of range (${MIN_N}–${MAX_N})`);
        }
        const expected = nRows * nCols;
        const entries = toks.slice(2, 2 + expected);
        if (entries.length !== expected) {
            throw new Error(`Expected ${expected} entries after header, got ${entries.length}`);
        }
        return { n: nRows, entries };
    }

    // Wire up one <section class="copos-tool" data-mode="...">. Returns void;
    // each section keeps its DOM elements in a closure so they don't collide.
    function setupSection(root) {
        const modeKey = root.dataset.mode;
        const mode = MODES[modeKey];
        if (!mode) {
            console.warn("copos-ui: section has unknown data-mode:", modeKey, root);
            return;
        }

        const els = {
            dim: root.querySelector(".copos-dim"),
            grid: root.querySelector(".copos-grid"),
            run: root.querySelector(".copos-run"),
            symmetrize: root.querySelector(".copos-symmetrize"),
            presets: root.querySelector(".copos-presets"),
            file: root.querySelector(".copos-file"),
            status: root.querySelector(".copos-status"),
            result: root.querySelector(".copos-result"),
        };

        function rebuildGrid(n) {
            els.grid.innerHTML = "";
            els.grid.style.gridTemplateColumns = `repeat(${n}, minmax(3.5rem, 1fr))`;
            for (let i = 0; i < n; i++) {
                for (let j = 0; j < n; j++) {
                    const inp = document.createElement("input");
                    inp.type = "text";
                    inp.value = i === j ? "1" : "0";
                    inp.dataset.i = i;
                    inp.dataset.j = j;
                    inp.setAttribute("aria-label", `A[${i + 1},${j + 1}]`);
                    els.grid.appendChild(inp);
                }
            }
        }

        function readMatrix() {
            const n = parseInt(els.dim.value, 10);
            const entries = new Array(n * n);
            els.grid.querySelectorAll("input").forEach((inp) => {
                const i = +inp.dataset.i, j = +inp.dataset.j;
                entries[i * n + j] = inp.value.trim();
            });
            return { n, entries };
        }

        function writeMatrix(n, entries) {
            if (parseInt(els.dim.value, 10) !== n) {
                els.dim.value = String(n);
                rebuildGrid(n);
            }
            els.grid.querySelectorAll("input").forEach((inp) => {
                const i = +inp.dataset.i, j = +inp.dataset.j;
                inp.value = entries[i * n + j];
            });
        }

        function symmetrize() {
            const { n, entries } = readMatrix();
            for (let i = 0; i < n; i++) {
                for (let j = i + 1; j < n; j++) {
                    entries[j * n + i] = entries[i * n + j];
                }
            }
            writeMatrix(n, entries);
        }

        function setStatus(msg, kind) {
            els.status.textContent = msg;
            els.status.className = "copos-status" + (kind ? " " + kind : "");
        }

        function renderResult(r) {
            els.result.innerHTML = mode.renderer(r, { escapeHtml });
        }

        async function runTest() {
            els.result.innerHTML = "";
            setStatus("Loading module…");
            try {
                const m = await loadModule();
                const { n, entries } = readMatrix();
                setStatus("Computing…");
                const t0 = performance.now();
                const r = m[mode.method](n, entries);
                const elapsedMs = performance.now() - t0;
                if (r.error) {
                    els.result.innerHTML = "";
                    setStatus(r.error, "error");
                    return;
                }
                renderResult(r);
                setStatus(`Done in ${formatElapsed(elapsedMs)}.`, "ok");
            } catch (err) {
                setStatus("Error: " + (err && err.message ? err.message : String(err)), "error");
                console.error(err);
            }
        }

        function populatePresets() {
            for (const key of Object.keys(PRESETS)) {
                const opt = document.createElement("option");
                opt.value = key;
                opt.textContent = PRESETS[key].label;
                els.presets.appendChild(opt);
            }
        }

        function applyPreset(key) {
            if (!key) return;
            const p = PRESETS[key];
            const n = p.n || parseInt(els.dim.value, 10);
            if (n < MIN_N || n > MAX_N) {
                setStatus(`Preset requires n=${n} (out of range)`, "error");
                return;
            }
            writeMatrix(n, p.fill(n));
            setStatus("");
        }

        function loadFile(file) {
            const reader = new FileReader();
            reader.onload = () => {
                try {
                    const { n, entries } = parseMatrixFile(String(reader.result));
                    writeMatrix(n, entries);
                    els.result.innerHTML = "";
                    setStatus(`Loaded ${n} × ${n} matrix from ${file.name}`, "ok");
                } catch (err) {
                    setStatus("File error: " + err.message, "error");
                }
            };
            reader.onerror = () => setStatus("Could not read file", "error");
            reader.readAsText(file);
        }

        // Initial population.
        for (let n = MIN_N; n <= MAX_N; n++) {
            const opt = document.createElement("option");
            opt.value = n;
            opt.textContent = `${n} × ${n}`;
            if (n === 3) opt.selected = true;
            els.dim.appendChild(opt);
        }
        rebuildGrid(parseInt(els.dim.value, 10));
        populatePresets();

        els.dim.addEventListener("change", () => {
            rebuildGrid(parseInt(els.dim.value, 10));
            els.result.innerHTML = "";
            setStatus("");
        });
        els.symmetrize.addEventListener("click", symmetrize);
        els.run.addEventListener("click", runTest);
        els.presets.addEventListener("change", (e) => {
            applyPreset(e.target.value);
            e.target.value = "";
        });
        els.file.addEventListener("change", (e) => {
            const f = e.target.files && e.target.files[0];
            if (f) loadFile(f);
            e.target.value = "";
        });
    }

    function init() {
        document.querySelectorAll(".copos-tool").forEach(setupSection);
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();
