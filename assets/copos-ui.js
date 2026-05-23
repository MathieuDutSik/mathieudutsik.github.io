// Tools page glue: build the n×n input grid, call the WASM copositivity
// test, render the result. Depends on /wasm/copos.js (Emscripten module
// factory) being loaded first as a classic <script> tag.

(function () {
    "use strict";

    const MIN_N = 2;
    const MAX_N = 8;
    const PRESETS = {
        "identity": {
            label: "Identity I_n (copositive)",
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
            label: "Horn matrix H_5 (copositive but not PSD + nonneg)",
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

    const els = {
        dim: document.getElementById("copos-dim"),
        grid: document.getElementById("copos-grid"),
        run: document.getElementById("copos-run"),
        symmetrize: document.getElementById("copos-symmetrize"),
        presets: document.getElementById("copos-presets"),
        file: document.getElementById("copos-file"),
        status: document.getElementById("copos-status"),
        result: document.getElementById("copos-result"),
    };

    let modulePromise = null;
    function loadModule() {
        if (modulePromise) return modulePromise;
        if (typeof createCoposModule !== "function") {
            return Promise.reject(new Error(
                "copos.js (Emscripten module) failed to load"));
        }
        modulePromise = createCoposModule({
            locateFile: (path) => path.endsWith(".wasm")
                ? "/wasm/" + path
                : path,
        });
        return modulePromise;
    }

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
        const inputs = els.grid.querySelectorAll("input");
        inputs.forEach((inp) => {
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
        const inputs = els.grid.querySelectorAll("input");
        inputs.forEach((inp) => {
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
        const out = [];
        if (r.isCopositive) {
            out.push('<p class="verdict verdict-yes">Matrix is copositive.</p>');
            out.push(`<p>Reason: <code>${escapeHtml(r.nature)}</code></p>`);
        } else {
            out.push('<p class="verdict verdict-no">Matrix is NOT copositive.</p>');
            out.push(`<p>Nature of violation: <code>${escapeHtml(r.nature)}</code></p>`);
            if (r.witness && r.witness.size && r.witness.size() > 0) {
                const v = [];
                for (let i = 0; i < r.witness.size(); i++) v.push(r.witness.get(i));
                out.push(`<p>Witness vector v ≥ 0 with v<sup>T</sup>Av &lt; 0:</p>`);
                out.push(`<pre>v = (${v.map(escapeHtml).join(", ")})</pre>`);
            }
        }
        els.result.innerHTML = out.join("");
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

    async function runTest() {
        els.result.innerHTML = "";
        setStatus("Loading module…");
        try {
            const m = await loadModule();
            const { n, entries } = readMatrix();
            setStatus("Computing…");
            const t0 = performance.now();
            const r = m.testCopositivity(n, entries);
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

    // Parse the standard "nRows nCols\n<entries>" format used by polyhedral_common
    // command-line tools. Whitespace is any combination of spaces, tabs, newlines.
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

    function init() {
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

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();
