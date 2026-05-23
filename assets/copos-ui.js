// Tools page glue: for each <details class="copos-tool" data-mode="..."> on
// the page, build a matrix-input grid, call the corresponding WASM binding,
// and render the result. Each section is independent — its own dimensions,
// matrix, presets, status, and result.
//
// Shape (driven by data-shape on the section root):
//   "square"  — one n × n grid (default; used by copositivity-type tests).
//   "rect"    — separate rows / cols inputs (used by SHV inputs).
//
// Depends on /wasm/copos.js (Emscripten module factory) being loaded first.

(function () {
    "use strict";

    const MIN_N = 1;
    const MAX_N = 50;

    const PRESETS = {
        "identity": {
            label: "Identity I_n",
            shape: "square",
            fill: (rows, cols) => Array.from({ length: rows * cols }, (_, k) =>
                Math.floor(k / cols) === (k % cols) ? "1" : "0"),
        },
        "neg-diag": {
            label: "diag(-1, 1, ..., 1) (NOT copositive)",
            shape: "square",
            fill: (rows, cols) => Array.from({ length: rows * cols }, (_, k) => {
                const i = Math.floor(k / cols), j = k % cols;
                if (i !== j) return "0";
                return i === 0 ? "-1" : "1";
            }),
        },
        "horn": {
            label: "Horn matrix H_5 (copositive, NOT completely positive)",
            shape: "square",
            rows: 5, cols: 5,
            fill: () => [
                "1", "-1", "1", "1", "-1",
                "-1", "1", "-1", "1", "1",
                "1", "-1", "1", "-1", "1",
                "1", "1", "-1", "1", "-1",
                "-1", "1", "1", "-1", "1",
            ],
        },
    };

    // Per-mode: which WASM function, and how to render the result.
    // `shape` defaults to "square" if not specified.
    const MODES = {
        copositivity: {
            method: "testCopositivity",
            shape: "square",
            renderer: (r, ctx) => renderCopositivity(r, "copositive",
                "v ≥ 0 with v<sup>T</sup>Av &lt; 0", ctx),
        },
        factorization: {
            method: "testCopositiveFactorization",
            shape: "square",
            renderer: (r, ctx) => renderFactorization(r, ctx),
        },
        shvRealizability: {
            method: "testShortestVectorsRealizability",
            shape: "rect",
            integerOnly: true,
            renderer: (r, ctx) => renderRealizability(r, ctx),
        },
        shvAutomorphism: {
            method: "shortestVectorsAutomorphismGroup",
            shape: "rect",
            integerOnly: true,
            renderer: (r, ctx) => renderAutomorphism(r, ctx),
        },
        gramCanonical: {
            method: "gramCanonicalForm",
            shape: "square",
            renderer: (r, ctx) => renderCanonical(r, ctx),
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

    function renderMatrixPre(flat, rows, cols, ctx) {
        const lines = [];
        for (let i = 0; i < rows; i++) {
            const row = [];
            for (let j = 0; j < cols; j++) row.push(flat.get(i * cols + j));
            lines.push("  " + row.map(ctx.escapeHtml).join("  "));
        }
        return lines.join("\n");
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
                out.push(`<pre>C =\n${renderMatrixPre(r.certificateFlat, n, n, ctx)}</pre>`);
            }
        }
        return out.join("");
    }

    function renderRealizability(r, ctx) {
        const out = [];
        if (r.realizable) {
            out.push(`<p class="verdict verdict-yes">Configuration is realizable as a shortest-vector set.</p>`);
            out.push(`<p>A positive-definite Gram matrix realizing it:</p>`);
            out.push(`<pre>G =\n${renderMatrixPre(r.gramFlat, r.dim, r.dim, ctx)}</pre>`);
        } else {
            out.push(`<p class="verdict verdict-no">Configuration is NOT realizable as a shortest-vector set.</p>`);
        }
        return out.join("");
    }

    function renderAutomorphism(r, ctx) {
        const out = [];
        const n = r.dim;
        const nb = r.nGenerators;
        out.push(`<p class="verdict verdict-yes">Automorphism group computed.</p>`);
        out.push(`<p>${nb} generator${nb === 1 ? "" : "s"} (each is an n × n integer matrix acting on Z<sup>${n}</sup>):</p>`);
        for (let k = 0; k < nb; k++) {
            const lines = [];
            for (let i = 0; i < n; i++) {
                const row = [];
                for (let j = 0; j < n; j++) row.push(r.generatorsFlat.get(k * n * n + i * n + j));
                lines.push("  " + row.map(ctx.escapeHtml).join("  "));
            }
            out.push(`<pre>g<sub>${k + 1}</sub> =\n${lines.join("\n")}</pre>`);
        }
        return out.join("");
    }

    function renderCanonical(r, ctx) {
        const out = [];
        const n = r.dim;
        out.push(`<p class="verdict verdict-yes">Canonical form computed.</p>`);
        out.push(`<p>Canonical Gram matrix B · G · B<sup>T</sup>:</p>`);
        out.push(`<pre>${renderMatrixPre(r.canonicalFlat, n, n, ctx)}</pre>`);
        out.push(`<p>Canonicalizing integer basis B (rows of B give the new basis in the original coordinates):</p>`);
        out.push(`<pre>${renderMatrixPre(r.basisFlat, n, n, ctx)}</pre>`);
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

    // Parse the polyhedral_common matrix file format:
    //   nRows nCols
    //   A11 A12 ...
    //   ...
    function parseMatrixFile(text, shape) {
        const toks = text.split(/\s+/).filter(Boolean);
        if (toks.length < 2) throw new Error("File must start with: nRows nCols");
        const nRows = parseInt(toks[0], 10);
        const nCols = parseInt(toks[1], 10);
        if (!Number.isInteger(nRows) || !Number.isInteger(nCols) || nRows <= 0 || nCols <= 0) {
            throw new Error("Header (nRows nCols) must be two positive integers");
        }
        if (shape === "square" && nRows !== nCols) {
            throw new Error(`This tool needs a square matrix; got ${nRows} × ${nCols}`);
        }
        if (nRows > MAX_N || nCols > MAX_N) {
            throw new Error(`Dimensions exceed ${MAX_N} (got ${nRows} × ${nCols})`);
        }
        const expected = nRows * nCols;
        const entries = toks.slice(2, 2 + expected);
        if (entries.length !== expected) {
            throw new Error(`Expected ${expected} entries after header, got ${entries.length}`);
        }
        return { rows: nRows, cols: nCols, entries };
    }

    // Wire up one <details class="copos-tool" data-mode="...">. Returns void;
    // each section keeps its DOM elements + state in a closure.
    function setupSection(root) {
        const modeKey = root.dataset.mode;
        const mode = MODES[modeKey];
        if (!mode) {
            console.warn("copos-ui: section has unknown data-mode:", modeKey, root);
            return;
        }
        const shape = mode.shape || "square";

        const els = {
            rows: root.querySelector(".copos-rows"),
            cols: root.querySelector(".copos-cols"),
            grid: root.querySelector(".copos-grid"),
            run: root.querySelector(".copos-run"),
            symmetrize: root.querySelector(".copos-symmetrize"),
            presets: root.querySelector(".copos-presets"),
            file: root.querySelector(".copos-file"),
            status: root.querySelector(".copos-status"),
            result: root.querySelector(".copos-result"),
        };

        function readDims() {
            // For "square" the rows input drives both; cols input may be absent.
            const r = parseInt(els.rows.value, 10);
            if (shape === "square") return { rows: r, cols: r };
            return { rows: r, cols: parseInt(els.cols.value, 10) };
        }

        function setDims(rows, cols) {
            els.rows.value = String(rows);
            if (els.cols) els.cols.value = String(cols);
        }

        function clampDim(n) {
            if (!Number.isFinite(n) || n < MIN_N) return MIN_N;
            if (n > MAX_N) return MAX_N;
            return n;
        }

        function rebuildGrid(rows, cols) {
            els.grid.innerHTML = "";
            els.grid.style.gridTemplateColumns = `repeat(${cols}, minmax(3rem, 1fr))`;
            for (let i = 0; i < rows; i++) {
                for (let j = 0; j < cols; j++) {
                    const inp = document.createElement("input");
                    inp.type = "text";
                    inp.value = (shape === "square" && i === j) ? "1" : "0";
                    inp.dataset.i = i;
                    inp.dataset.j = j;
                    inp.setAttribute("aria-label", `M[${i + 1},${j + 1}]`);
                    els.grid.appendChild(inp);
                }
            }
        }

        function readMatrix() {
            const { rows, cols } = readDims();
            const entries = new Array(rows * cols);
            els.grid.querySelectorAll("input").forEach((inp) => {
                const i = +inp.dataset.i, j = +inp.dataset.j;
                entries[i * cols + j] = inp.value.trim();
            });
            return { rows, cols, entries };
        }

        function writeMatrix(rows, cols, entries) {
            const current = readDims();
            if (current.rows !== rows || current.cols !== cols) {
                setDims(rows, cols);
                rebuildGrid(rows, cols);
            }
            els.grid.querySelectorAll("input").forEach((inp) => {
                const i = +inp.dataset.i, j = +inp.dataset.j;
                inp.value = entries[i * cols + j];
            });
        }

        function symmetrize() {
            const { rows, cols, entries } = readMatrix();
            if (rows !== cols) {
                setStatus("Symmetrize only applies to a square matrix.", "error");
                return;
            }
            for (let i = 0; i < rows; i++) {
                for (let j = i + 1; j < cols; j++) {
                    entries[j * cols + i] = entries[i * cols + j];
                }
            }
            writeMatrix(rows, cols, entries);
        }

        function setStatus(msg, kind) {
            els.status.textContent = msg;
            els.status.className = "copos-status" + (kind ? " " + kind : "");
        }

        async function runTest() {
            els.result.innerHTML = "";
            setStatus("Loading module…");
            try {
                const m = await loadModule();
                const { rows, cols, entries } = readMatrix();
                setStatus("Computing…");
                const t0 = performance.now();
                const r = shape === "square"
                    ? m[mode.method](rows, entries)
                    : m[mode.method](rows, cols, entries);
                const elapsedMs = performance.now() - t0;
                if (r.error) {
                    els.result.innerHTML = "";
                    setStatus(r.error, "error");
                    return;
                }
                els.result.innerHTML = mode.renderer(r, { escapeHtml });
                setStatus(`Done in ${formatElapsed(elapsedMs)}.`, "ok");
            } catch (err) {
                setStatus("Error: " + (err && err.message ? err.message : String(err)), "error");
                console.error(err);
            }
        }

        function populatePresets() {
            if (!els.presets) return;
            for (const key of Object.keys(PRESETS)) {
                const p = PRESETS[key];
                // Only offer presets compatible with this section's shape.
                if (p.shape && p.shape !== shape && !(shape === "rect" && p.shape === "square")) {
                    continue;
                }
                const opt = document.createElement("option");
                opt.value = key;
                opt.textContent = p.label;
                els.presets.appendChild(opt);
            }
        }

        function applyPreset(key) {
            if (!key) return;
            const p = PRESETS[key];
            const cur = readDims();
            const rows = p.rows || cur.rows;
            const cols = p.cols || cur.cols;
            if (rows < MIN_N || cols < MIN_N || rows > MAX_N || cols > MAX_N) {
                setStatus(`Preset dimensions out of range`, "error");
                return;
            }
            writeMatrix(rows, cols, p.fill(rows, cols));
            setStatus("");
        }

        function loadFile(file) {
            const reader = new FileReader();
            reader.onload = () => {
                try {
                    const parsed = parseMatrixFile(String(reader.result), shape);
                    writeMatrix(parsed.rows, parsed.cols, parsed.entries);
                    els.result.innerHTML = "";
                    setStatus(`Loaded ${parsed.rows} × ${parsed.cols} matrix from ${file.name}`, "ok");
                } catch (err) {
                    setStatus("File error: " + err.message, "error");
                }
            };
            reader.onerror = () => setStatus("Could not read file", "error");
            reader.readAsText(file);
        }

        // Initial population.
        const defaultRows = parseInt(els.rows.value, 10) || 3;
        const defaultCols = (shape === "square")
            ? defaultRows
            : (parseInt(els.cols.value, 10) || defaultRows);
        setDims(defaultRows, defaultCols);
        rebuildGrid(defaultRows, defaultCols);
        populatePresets();

        const onDimChange = () => {
            const { rows, cols } = readDims();
            const r = clampDim(rows);
            const c = clampDim(cols);
            setDims(r, c);
            rebuildGrid(r, c);
            els.result.innerHTML = "";
            setStatus("");
        };
        els.rows.addEventListener("change", onDimChange);
        if (els.cols) els.cols.addEventListener("change", onDimChange);

        if (els.symmetrize) {
            els.symmetrize.addEventListener("click", symmetrize);
        }
        els.run.addEventListener("click", runTest);
        if (els.presets) {
            els.presets.addEventListener("change", (e) => {
                applyPreset(e.target.value);
                e.target.value = "";
            });
        }
        if (els.file) {
            els.file.addEventListener("change", (e) => {
                const f = e.target.files && e.target.files[0];
                if (f) loadFile(f);
                e.target.value = "";
            });
        }
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
