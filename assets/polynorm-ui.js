// Tools page glue for the polytope tools: the <details class="polynorm-tool"
// data-mode="..."> sections. In each, the polytope is given by M vectors v_i
// in R^n, read as its vertices, as the inequalities v_i . x <= 1, or as the
// generators of the zonotope v_1 [0,1] + ... + v_M [0,1]. The modes:
//   densities         packing and covering densities  (polytopeNormDensities)
//   isMinimalCovering minimal covering test           (isMinimalCovering)
//   minimalCovering   a minimal covering polytope     (minimalCoveringPolytope)
//
// The computations run in a Web Worker (assets/polynorm-worker.js) calling
// those bindings of /wasm/polyhedral.js, one worker per section.

(function () {
    "use strict";

    const MAX_ROWS = 50;
    const MAX_COLS = 10;

    // Presets for the densities tool.
    const DENSITY_PRESETS = {
        "square": {
            label: "Square [0,1]² (vertices)",
            inputType: "vertices", rows: 4, cols: 2,
            entries: ["0", "0", "1", "0", "0", "1", "1", "1"],
        },
        "triangle": {
            label: "Triangle conv(0, e₁, e₂) (vertices)",
            inputType: "vertices", rows: 3, cols: 2,
            entries: ["0", "0", "1", "0", "0", "1"],
        },
        "hexagon": {
            label: "Hexagon e₁[0,1] + e₂[0,1] + (e₁+e₂)[0,1] (zonotope)",
            inputType: "zonotope", rows: 3, cols: 2,
            entries: ["1", "0", "0", "1", "1", "1"],
        },
        "pentagon": {
            label: "Pentagon, 5 inequalities",
            inputType: "inequalities", rows: 5, cols: 2,
            entries: ["1", "0", "0", "1", "-1", "1", "-1", "-1", "1", "-1"],
        },
        "cube": {
            label: "Cube [-1,1]³ (inequalities ±xᵢ ≤ 1)",
            inputType: "inequalities", rows: 6, cols: 3,
            entries: ["1", "0", "0", "-1", "0", "0", "0", "1", "0",
                "0", "-1", "0", "0", "0", "1", "0", "0", "-1"],
        },
        "octahedron": {
            label: "Octahedron conv(±eᵢ) (vertices)",
            inputType: "vertices", rows: 6, cols: 3,
            entries: ["1", "0", "0", "-1", "0", "0", "0", "1", "0",
                "0", "-1", "0", "0", "0", "1", "0", "0", "-1"],
        },
        "octahedron-ineq": {
            label: "Octahedron |x|+|y|+|z| ≤ 1 (inequalities)",
            inputType: "inequalities", rows: 8, cols: 3,
            entries: ["1", "1", "1", "1", "1", "-1", "1", "-1", "1", "1", "-1", "-1",
                "-1", "1", "1", "-1", "1", "-1", "-1", "-1", "1", "-1", "-1", "-1"],
        },
        "tetrahedron": {
            label: "Tetrahedron conv(0, e₁, e₂, e₃) (vertices)",
            inputType: "vertices", rows: 4, cols: 3,
            entries: ["0", "0", "0", "1", "0", "0", "0", "1", "0", "0", "0", "1"],
        },
        "rhombic-dodecahedron": {
            label: "Rhombic dodecahedron (zonotope, 4 generators)",
            inputType: "zonotope", rows: 4, cols: 3,
            entries: ["1", "1", "1", "1", "1", "-1", "1", "-1", "1", "-1", "1", "1"],
        },
    };

    // Presets for the two minimal covering tools.
    const MINCOV_PRESETS = {
        "square": {
            label: "Square [0,1]² (minimal covering, a tile)",
            inputType: "vertices", rows: 4, cols: 2,
            entries: ["0", "0", "1", "0", "0", "1", "1", "1"],
        },
        "rectangle": {
            label: "Rectangle [0,3/2] × [0,1] (covering, not minimal)",
            inputType: "vertices", rows: 4, cols: 2,
            entries: ["0", "0", "3/2", "0", "0", "1", "3/2", "1"],
        },
        "triangle-1": {
            label: "Triangle conv(0, e₁, e₂) (not covering, μ = 2)",
            inputType: "vertices", rows: 3, cols: 2,
            entries: ["0", "0", "1", "0", "0", "1"],
        },
        "triangle-3": {
            label: "Triangle conv(0, 3e₁, 3e₂)",
            inputType: "vertices", rows: 3, cols: 2,
            entries: ["0", "0", "3", "0", "0", "3"],
        },
        "hexagon": {
            label: "Hexagon conv(±e₁, ±e₂, ±(e₁+e₂))",
            inputType: "vertices", rows: 6, cols: 2,
            entries: ["1", "0", "-1", "0", "0", "1", "0", "-1", "1", "1", "-1", "-1"],
        },
        "tetrahedron-3": {
            label: "Tetrahedron conv(0, 3e₁, 3e₂, 3e₃)",
            inputType: "vertices", rows: 4, cols: 3,
            entries: ["0", "0", "0", "3", "0", "0", "0", "3", "0", "0", "0", "3"],
        },
        "octahedron-3-2": {
            label: "Octahedron (3/2) conv(±eᵢ) (finding Q: about 20 s)",
            inputType: "vertices", rows: 6, cols: 3,
            entries: ["3/2", "0", "0", "-3/2", "0", "0", "0", "3/2", "0",
                "0", "-3/2", "0", "0", "0", "3/2", "0", "0", "-3/2"],
        },
        "octahedron-q": {
            label: "Minimal covering Q found in (3/2) conv(±eᵢ) (volume 5/4, not a tile)",
            inputType: "vertices", rows: 12, cols: 3,
            entries: ["-1/2", "1/2", "-1/2", "0", "-1/2", "0", "1/2", "0", "0",
                "-1", "1/2", "0", "-1/2", "1", "0", "-1/2", "-1/2", "1/2",
                "1/2", "-1/2", "1/2", "1/2", "1/2", "1/2", "0", "-1/2", "1",
                "-1/2", "0", "1", "1/2", "0", "1", "0", "1/2", "1"],
        },
    };

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

    // Decimal value of an exact rational "p/q" (or integer "p"), to `digits`
    // decimals, computed with BigInt so that large numerators stay exact.
    function rationalToDecimal(s, digits) {
        const parts = String(s).split("/");
        let p = BigInt(parts[0]);
        const q = parts.length > 1 ? BigInt(parts[1]) : 1n;
        const neg = p < 0n;
        if (neg) p = -p;
        const scale = 10n ** BigInt(digits);
        const scaled = (2n * p * scale + q) / (2n * q);  // rounded
        const intPart = scaled / scale;
        let frac = (scaled % scale).toString().padStart(digits, "0").replace(/0+$/, "");
        return (neg ? "-" : "") + intPart.toString() + (frac ? "." + frac : "");
    }

    function exactAndDecimal(s) {
        const exact = escapeHtml(s);
        if (!String(s).includes("/")) return `<code>${exact}</code>`;
        return `<code>${exact}</code> ≈ ${escapeHtml(rationalToDecimal(s, 6))}`;
    }

    function vectorList(flat, n, name) {
        const lines = [];
        for (let i = 0; i * n < flat.length; i++) {
            lines.push(`${name}${i + 1} = (${flat.slice(i * n, (i + 1) * n).map(escapeHtml).join(", ")})`);
        }
        return lines.join("\n");
    }

    function renderResult(r) {
        const n = r.dim;
        const out = [];
        out.push(`<table class="copos-factor-table polynorm-summary">`);
        out.push(`<tr><td>Dimension</td><td>${n}</td></tr>`);
        out.push(`<tr><td>Vertices / facets of P</td><td>${r.nVertices} / ${r.nFacets}</td></tr>`);
        out.push(`<tr><td>Volume of P</td><td>${exactAndDecimal(r.volume)}</td></tr>`);
        out.push(`<tr><td>Order of {A ∈ GL<sub>n</sub>(ℤ) : PA = P}</td><td>${escapeHtml(r.groupOrder)}</td></tr>`);
        if (r.hasPacking) {
            out.push(`<tr><td>Packing scalar α</td><td>${exactAndDecimal(r.alpha)}</td></tr>`);
            out.push(`<tr><td><strong>Packing density</strong> α<sup>n</sup> vol(P)</td><td>${exactAndDecimal(r.packingDensity)}</td></tr>`);
        }
        if (r.hasCovering) {
            out.push(`<tr><td>Covering radius μ</td><td>${exactAndDecimal(r.mu)}</td></tr>`);
            out.push(`<tr><td><strong>Covering density</strong> μ<sup>n</sup> vol(P)</td><td>${exactAndDecimal(r.coveringDensity)}</td></tr>`);
        }
        out.push(`</table>`);
        out.push(`<details><summary>Vertices of P (${r.nVertices})</summary>`);
        out.push(`<pre>${vectorList(r.verticesFlat, n, "w")}</pre></details>`);
        if (r.hasPacking) {
            out.push(`<details><summary>Contact vectors z, αP and αP + z touch (${r.nContacts})</summary>`);
            out.push(`<pre>${vectorList(r.contactsFlat, n, "z")}</pre></details>`);
        }
        if (r.hasCovering) {
            out.push(`<details><summary>Last covered point and the lattice points pinning p − μP (${r.nTight})</summary>`);
            out.push(`<pre>p = (${r.lastCoveredPoint.map(escapeHtml).join(", ")})\n\n` +
                `${vectorList(r.tightFlat, n, "z")}</pre></details>`);
        }
        return out.join("");
    }

    function verdict(ok, text) {
        return `<p class="verdict ${ok ? "verdict-yes" : "verdict-no"}">${text}</p>`;
    }

    function renderMinimalTest(r) {
        const n = r.dim;
        const out = [];
        if (!r.isCovering) {
            out.push(verdict(false, "P is NOT covering: P + ℤ<sup>n</sup> ≠ ℝ<sup>n</sup>."));
            out.push(`<p>Its covering radius is μ = ${exactAndDecimal(r.mu)} &gt; 1, so μP is covering; ` +
                `the minimal covering polytope tool works on μP.</p>`);
        } else if (r.isMinimal) {
            out.push(verdict(true, "P is minimal covering."));
            out.push(`<p>P + ℤ<sup>n</sup> = ℝ<sup>n</sup>, and no convex body strictly contained in P is covering: ` +
                `no vertex of P can be cut.</p>`);
        } else {
            out.push(verdict(false, "P is covering but NOT minimal."));
            out.push(`<p>${r.nCuttable} vert${r.nCuttable === 1 ? "ex" : "ices"} of P can be cut keeping the covering:</p>`);
            out.push(`<pre>${vectorList(r.cuttableFlat, n, "w")}</pre>`);
        }
        out.push(`<table class="copos-factor-table">`);
        out.push(`<tr><td>Vertices / facets of P</td><td>${r.nVertices} / ${r.nFacets}</td></tr>`);
        out.push(`<tr><td>Volume of P</td><td>${exactAndDecimal(r.volume)}</td></tr>`);
        out.push(`</table>`);
        out.push(`<details><summary>Vertices of P (${r.nVertices})</summary>`);
        out.push(`<pre>${vectorList(r.verticesFlat, n, "w")}</pre></details>`);
        return out.join("");
    }

    function renderMinimalCovering(r) {
        const n = r.dim;
        const out = [];
        if (!r.isCovering) {
            out.push(verdict(false, "P is NOT covering."));
            out.push(`<p>Its covering radius is μ = ${exactAndDecimal(r.muInput)} &gt; 1. ` +
                `Scale P by μ (or more) to get a covering polytope.</p>`);
            return out.join("");
        }
        if (r.isMinimal) {
            out.push(verdict(true, "Found a minimal covering polytope Q ⊂ μP." +
                (r.isTile ? " Q has volume 1: it is a lattice tile." : "")));
        } else {
            out.push(verdict(false, "Found a covering polytope Q ⊂ μP, NOT proved minimal."));
            out.push(`<p>Some vertex of Q can still be cut: try a finer grid, more rounds, or more reductions.</p>`);
        }
        out.push(`<table class="copos-factor-table">`);
        out.push(`<tr><td>Covering radius μ of P</td><td>${exactAndDecimal(r.muInput)}</td></tr>`);
        out.push(`<tr><td>Vertices / facets of Q</td><td>${r.nVertices} / ${r.nFacets}</td></tr>`);
        out.push(`<tr><td>Volume of Q</td><td>${exactAndDecimal(r.volume)}</td></tr>`);
        out.push(`</table>`);
        out.push(`<p>${escapeHtml(r.description)}</p>`);
        out.push(`<details open><summary>Vertices of Q (${r.nVertices})</summary>`);
        out.push(`<pre>${vectorList(r.verticesFlat, n, "w")}</pre></details>`);
        const facets = [];
        for (let i = 0; i < r.nFacets; i++) {
            const row = r.facetsFlat.slice(i * (n + 1), (i + 1) * (n + 1));
            const terms = row.slice(1).map((a, k) => `${a}·x${k + 1}`).join(" + ");
            facets.push(`${row[0]} + ${terms} ≥ 0`.replace(/\+ -/g, "− "));
        }
        out.push(`<details><summary>Facets of Q (${r.nFacets}), as b + a·x ≥ 0</summary>`);
        out.push(`<pre>${facets.map(escapeHtml).join("\n")}</pre></details>`);
        return out.join("");
    }

    function intInput(el, min, fallback) {
        const v = parseInt(el.value, 10);
        return Number.isInteger(v) && v >= min ? v : fallback;
    }

    // Per mode: the presets, the default preset, the binding and its
    // arguments after (inputType, rows, cols, entries), and the renderer.
    const MODES = {
        densities: {
            presets: DENSITY_PRESETS, initial: "square",
            method: "polytopeNormDensities",
            extraArgs: (root) => [root.querySelector(".polynorm-what").value],
            render: renderResult,
        },
        isMinimalCovering: {
            presets: MINCOV_PRESETS, initial: "rectangle",
            method: "isMinimalCovering",
            extraArgs: () => [],
            render: renderMinimalTest,
        },
        minimalCovering: {
            presets: MINCOV_PRESETS, initial: "rectangle",
            method: "minimalCoveringPolytope",
            extraArgs: (root) => {
                const method = root.querySelector(".polynorm-method").value;
                if (method === "grid") {
                    return [method, intInput(root.querySelector(".polynorm-grid-n"), 1, 1),
                        intInput(root.querySelector(".polynorm-rounds"), 1, 3)];
                }
                return [method, intInput(root.querySelector(".polynorm-max-step"), 1, 1000), 0];
            },
            render: renderMinimalCovering,
        },
    };

    // Matrix file: "M n" followed by the M vectors, one per row.
    function parseVectorFile(text) {
        const toks = text.split(/\s+/).filter(Boolean);
        if (toks.length < 2) throw new Error("File must start with: M n");
        const rows = parseInt(toks[0], 10);
        const cols = parseInt(toks[1], 10);
        if (!Number.isInteger(rows) || !Number.isInteger(cols) || rows <= 0 || cols <= 0) {
            throw new Error("Header (M n) must be two positive integers");
        }
        if (rows > MAX_ROWS || cols > MAX_COLS) {
            throw new Error(`At most ${MAX_ROWS} vectors in dimension at most ${MAX_COLS} (got ${rows} × ${cols})`);
        }
        const entries = toks.slice(2, 2 + rows * cols);
        if (entries.length !== rows * cols) {
            throw new Error(`Expected ${rows * cols} entries after header, got ${entries.length}`);
        }
        return { rows, cols, entries };
    }

    function setupSection(root) {
        const mode = MODES[root.dataset.mode];
        if (!mode) {
            console.warn("polynorm-ui: section has unknown data-mode:", root.dataset.mode, root);
            return;
        }
        const PRESETS = mode.presets;
        const els = {
            inputType: root.querySelector(".polynorm-input-type"),
            rows: root.querySelector(".copos-rows"),
            cols: root.querySelector(".copos-cols"),
            grid: root.querySelector(".copos-grid"),
            run: root.querySelector(".copos-run"),
            stop: root.querySelector(".polynorm-stop"),
            presets: root.querySelector(".copos-presets"),
            file: root.querySelector(".copos-file"),
            status: root.querySelector(".copos-status"),
            result: root.querySelector(".copos-result"),
        };
        let worker = null;

        function setStatus(msg, kind) {
            els.status.textContent = msg;
            els.status.className = "copos-status" + (kind ? " " + kind : "");
        }

        function clamp(v, max) {
            if (!Number.isFinite(v) || v < 1) return 1;
            return Math.min(v, max);
        }

        function readDims() {
            return {
                rows: clamp(parseInt(els.rows.value, 10), MAX_ROWS),
                cols: clamp(parseInt(els.cols.value, 10), MAX_COLS),
            };
        }

        function rebuildGrid(rows, cols, entries) {
            els.grid.innerHTML = "";
            els.grid.style.gridTemplateColumns = `auto repeat(${cols}, minmax(3rem, 1fr))`;
            for (let i = 0; i < rows; i++) {
                const label = document.createElement("span");
                label.className = "polynorm-row-label";
                label.innerHTML = `v<sub>${i + 1}</sub>`;
                els.grid.appendChild(label);
                for (let j = 0; j < cols; j++) {
                    const inp = document.createElement("input");
                    inp.type = "text";
                    inp.value = entries ? entries[i * cols + j] : (i === j ? "1" : "0");
                    inp.dataset.i = i;
                    inp.dataset.j = j;
                    inp.setAttribute("aria-label", `v${i + 1}[${j + 1}]`);
                    els.grid.appendChild(inp);
                }
            }
        }

        function readEntries(rows, cols) {
            const entries = new Array(rows * cols);
            els.grid.querySelectorAll("input").forEach((inp) => {
                entries[+inp.dataset.i * cols + +inp.dataset.j] = inp.value.trim();
            });
            return entries;
        }

        function load(rows, cols, entries) {
            els.rows.value = String(rows);
            els.cols.value = String(cols);
            rebuildGrid(rows, cols, entries);
            els.result.innerHTML = "";
        }

        function finish() {
            els.run.disabled = false;
            els.stop.disabled = true;
        }

        function run() {
            const { rows, cols } = readDims();
            const entries = readEntries(rows, cols);
            if (!worker) {
                worker = new Worker("/assets/polynorm-worker.js" +
                    (window.__WASM_BUST ? "?v=" + window.__WASM_BUST : ""));
            }
            els.result.innerHTML = "";
            els.run.disabled = true;
            els.stop.disabled = false;
            setStatus("Computing…");
            worker.onmessage = (e) => {
                finish();
                const msg = e.data;
                if (!msg.ok) {
                    setStatus("Error: " + msg.error, "error");
                    return;
                }
                if (msg.result.error) {
                    setStatus(msg.result.error, "error");
                    return;
                }
                els.result.innerHTML = mode.render(msg.result);
                setStatus(`Done in ${formatElapsed(msg.elapsedMs)}.`, "ok");
            };
            worker.onerror = (e) => {
                finish();
                setStatus("Error: " + (e.message || "the worker failed"), "error");
                worker.terminate();
                worker = null;
            };
            worker.postMessage({
                bust: window.__WASM_BUST ? "?v=" + window.__WASM_BUST : "",
                method: mode.method,
                args: [els.inputType.value, rows, cols, entries, ...mode.extraArgs(root)],
            });
        }

        function stop() {
            if (worker) {
                worker.terminate();
                worker = null;
            }
            finish();
            setStatus("Stopped.", "error");
        }

        for (const key of Object.keys(PRESETS)) {
            const opt = document.createElement("option");
            opt.value = key;
            opt.textContent = PRESETS[key].label;
            els.presets.appendChild(opt);
        }
        els.presets.addEventListener("change", (e) => {
            const p = PRESETS[e.target.value];
            e.target.value = "";
            if (!p) return;
            els.inputType.value = p.inputType;
            load(p.rows, p.cols, p.entries);
            setStatus("");
        });

        const onDimChange = () => {
            const { rows, cols } = readDims();
            load(rows, cols, null);
            setStatus("");
        };
        els.rows.addEventListener("change", onDimChange);
        els.cols.addEventListener("change", onDimChange);
        els.run.addEventListener("click", run);
        els.stop.addEventListener("click", stop);
        els.file.addEventListener("change", (e) => {
            const f = e.target.files && e.target.files[0];
            e.target.value = "";
            if (!f) return;
            const reader = new FileReader();
            reader.onload = () => {
                try {
                    const parsed = parseVectorFile(String(reader.result));
                    load(parsed.rows, parsed.cols, parsed.entries);
                    setStatus(`Loaded ${parsed.rows} vectors in dimension ${parsed.cols} from ${f.name}`, "ok");
                } catch (err) {
                    setStatus("File error: " + err.message, "error");
                }
            };
            reader.onerror = () => setStatus("Could not read file", "error");
            reader.readAsText(f);
        });

        // The grid method has parameters N and n_round, the reduction method
        // a maximal number of steps: show the ones of the chosen method.
        const methodSel = root.querySelector(".polynorm-method");
        if (methodSel) {
            const sync = () => {
                root.querySelectorAll("[data-method]").forEach((el) => {
                    el.hidden = el.dataset.method !== methodSel.value;
                });
            };
            methodSel.addEventListener("change", sync);
            sync();
        }

        const p = PRESETS[mode.initial];
        els.inputType.value = p.inputType;
        load(p.rows, p.cols, p.entries);
    }

    function init() {
        document.querySelectorAll(".polynorm-tool").forEach(setupSection);
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();
