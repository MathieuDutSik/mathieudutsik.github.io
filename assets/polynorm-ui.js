// Tools page glue for the polytope tools: the <details class="polynorm-tool"
// data-mode="..."> sections. In each, the polytope is given by M vectors v_i
// in R^n, read as its vertices, as the inequalities v_i . x <= 1, or as the
// generators of the zonotope v_1 [0,1] + ... + v_M [0,1]. The modes:
//   densities         packing and covering densities  (polytopeNormDensities)
//   isMinimalCovering minimal covering test           (isMinimalCovering)
//   minimalCovering   a minimal covering polytope     (minimalCoveringPolytope)
//   coveringPath      1/mu along a Minkowski path     (coveringPath), two
//                     polytopes K1, K2 in .polynorm-polytope blocks
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

    // Presets for the covering radius along a Minkowski path. The function
    // is f(lambda) = 1 / mu(lambda K1 + (1 - lambda) K2).
    const PATH_PRESETS = {
        "cfk": {
            label: "Two parallelogram tiles of Codenotti, Freyer, Krivokuca (4 pieces)",
            cols: 2,
            polytopes: [
                { inputType: "zonotope", rows: 2, entries: ["1", "1", "-3/2", "-1/2"] },
                { inputType: "zonotope", rows: 2, entries: ["1", "0", "1/2", "1"] },
            ],
        },
        "boxes": {
            label: "Boxes [0,1/2] × [0,2] and [0,2] × [0,1/2]",
            cols: 2,
            polytopes: [
                { inputType: "zonotope", rows: 2, entries: ["1/2", "0", "0", "2"] },
                { inputType: "zonotope", rows: 2, entries: ["2", "0", "0", "1/2"] },
            ],
        },
        "square-cross": {
            label: "Square |x|, |y| ≤ 1 and cross-polytope |x| + |y| ≤ 1 (inequalities)",
            cols: 2,
            polytopes: [
                { inputType: "inequalities", rows: 4, entries: ["1", "0", "-1", "0", "0", "1", "0", "-1"] },
                { inputType: "inequalities", rows: 4, entries: ["1", "1", "1", "-1", "-1", "1", "-1", "-1"] },
            ],
        },
        "boxes3": {
            label: "Boxes [0,1/2] × [0,2] × [0,1] and [0,2] × [0,1/2] × [0,1] (about 45 s)",
            cols: 3,
            polytopes: [
                { inputType: "zonotope", rows: 3, entries: ["1/2", "0", "0", "0", "2", "0", "0", "0", "1"] },
                { inputType: "zonotope", rows: 3, entries: ["2", "0", "0", "0", "1/2", "0", "0", "0", "1"] },
            ],
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

    // Exact rationals {p, q} over BigInt, q > 0, for the values of the
    // covering path at its break points.
    function bigAbs(a) { return a < 0n ? -a : a; }
    function bigGcd(a, b) {
        a = bigAbs(a); b = bigAbs(b);
        while (b) [a, b] = [b, a % b];
        return a;
    }
    function rat(p, q) {
        if (q < 0n) { p = -p; q = -q; }
        const g = bigGcd(p, q) || 1n;
        return { p: p / g, q: q / g };
    }
    function ratParse(s) {
        const parts = String(s).split("/");
        return rat(BigInt(parts[0]), parts.length > 1 ? BigInt(parts[1]) : 1n);
    }
    const ratAdd = (a, b) => rat(a.p * b.q + b.p * a.q, a.q * b.q);
    const ratMul = (a, b) => rat(a.p * b.p, a.q * b.q);
    const ratInv = (a) => rat(a.q, a.p);
    const ratCmp = (a, b) => (a.p * b.q < b.p * a.q ? -1 : a.p * b.q > b.p * a.q ? 1 : 0);
    const ratStr = (a) => a.q === 1n ? a.p.toString() : `${a.p}/${a.q}`;
    const ratNum = (a) => Number(a.p) / Number(a.q);

    // The pieces of the path, and its break points with their exact values.
    function pathData(r) {
        const pieces = [];
        for (let k = 0; k < r.nPieces; k++) {
            const [alpha, beta, s0, s1] = r.piecesFlat.slice(4 * k, 4 * k + 4).map(ratParse);
            pieces.push({ alpha, beta, s0, s1 });
        }
        const value = (P, l) => ratAdd(P.s0, ratMul(P.s1, l));
        const points = pieces.map((P) => ({ lambda: P.alpha, f: value(P, P.alpha) }));
        const last = pieces[pieces.length - 1];
        points.push({ lambda: last.beta, f: value(last, last.beta) });
        return { pieces, points, value };
    }

    function renderCoveringPath(r) {
        const { pieces, points, value } = pathData(r);
        let iMax = 0, iMin = 0;
        points.forEach((pt, i) => {
            if (ratCmp(pt.f, points[iMax].f) > 0) iMax = i;
            if (ratCmp(pt.f, points[iMin].f) < 0) iMin = i;
        });
        const exact = (a) => exactAndDecimal(ratStr(a));
        const out = [];
        out.push(verdict(true, `f(λ) = 1/μ(λK<sub>1</sub> + (1−λ)K<sub>2</sub>) has ` +
            `${pieces.length} affine piece${pieces.length === 1 ? "" : "s"} on [0,1], each one proved.`));
        out.push(`<table class="copos-factor-table">`);
        out.push(`<tr><td>μ(K<sub>2</sub>), at λ = 0</td><td>${exact(ratInv(points[0].f))}</td></tr>`);
        out.push(`<tr><td>μ(K<sub>1</sub>), at λ = 1</td><td>${exact(ratInv(points[points.length - 1].f))}</td></tr>`);
        out.push(`<tr><td>Largest μ (smallest 1/μ)</td><td>${exact(ratInv(points[iMin].f))} at λ = ${escapeHtml(ratStr(points[iMin].lambda))}</td></tr>`);
        out.push(`<tr><td>Smallest μ (largest 1/μ)</td><td>${exact(ratInv(points[iMax].f))} at λ = ${escapeHtml(ratStr(points[iMax].lambda))}</td></tr>`);
        out.push(`</table>`);
        out.push(`<div class="polynorm-chart" data-imax="${iMax}" data-imin="${iMin}"></div>`);
        const rows = pieces.map((P) => {
            const sign = P.s1.p < 0n ? "−" : "+";
            const s1abs = rat(bigAbs(P.s1.p), P.s1.q);
            return `<tr><td>[${escapeHtml(ratStr(P.alpha))}, ${escapeHtml(ratStr(P.beta))}]</td>` +
                `<td>${escapeHtml(ratStr(P.s0))} ${sign} ${escapeHtml(ratStr(s1abs))} λ</td>` +
                `<td>${escapeHtml(ratStr(ratInv(value(P, P.alpha))))} to ${escapeHtml(ratStr(ratInv(value(P, P.beta))))}</td></tr>`;
        });
        out.push(`<table class="copos-factor-table polynorm-pieces">` +
            `<tr><th>λ</th><th>1/μ</th><th>μ</th></tr>${rows.join("")}</table>`);
        out.push(`<p>${escapeHtml(r.description)}</p>`);
        return out.join("");
    }

    // The plot of f = 1/mu over [0,1]: the polyline through the break
    // points (f is affine in between), a marker at each break point, the
    // extreme values labeled, and a crosshair with a tooltip on hover.
    const SVG_NS = "http://www.w3.org/2000/svg";
    function svgEl(tag, attrs, text) {
        const el = document.createElementNS(SVG_NS, tag);
        for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
        if (text !== undefined) el.textContent = text;
        return el;
    }

    // A round step (1, 2, 2.5 or 5 times a power of 10) giving about
    // `target` intervals on [lo, hi].
    function niceStep(lo, hi, target) {
        const raw = (hi - lo) / target;
        const mag = Math.pow(10, Math.floor(Math.log10(raw)));
        return [1, 2, 2.5, 5, 10].map((m) => m * mag).find((st) => st >= raw * (1 - 1e-9));
    }

    function wirePathChart(container, r) {
        const host = container.querySelector(".polynorm-chart");
        const { pieces, points, value } = pathData(r);
        const W = 640, H = 300, ML = 56, MR = 20, MT = 24, MB = 44;
        const fs = points.map((pt) => ratNum(pt.f));
        let lo = Math.min(...fs), hi = Math.max(...fs);
        const pad = hi > lo ? (hi - lo) * 0.15 : 0.05;
        lo -= pad; hi += pad;
        // The range is widened to whole tick steps, so that both ends get a tick.
        const step = niceStep(lo, hi, 5);
        const k0 = Math.floor(lo / step + 1e-9), k1 = Math.ceil(hi / step - 1e-9);
        lo = k0 * step;
        hi = k1 * step;
        const ticks = [];
        for (let k = k0; k <= k1; k++) ticks.push(k * step);
        const X = (l) => ML + l * (W - ML - MR);
        const Y = (f) => MT + (hi - f) / (hi - lo) * (H - MT - MB);
        // As many decimals as the step has: 0.05 -> 2, 0.025 -> 3.
        const digits = (String(+step.toPrecision(12)).split(".")[1] || "").length;

        const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
            "aria-label": "Plot of 1/mu(lambda K1 + (1 - lambda) K2) for lambda in [0,1]; the values are in the table below" });
        for (const t of ticks) {
            svg.appendChild(svgEl("line", { x1: ML, x2: W - MR, y1: Y(t), y2: Y(t), class: "grid" }));
            svg.appendChild(svgEl("text", { x: ML - 8, y: Y(t) + 4, "text-anchor": "end", class: "tick" }, t.toFixed(digits)));
        }
        for (const l of [0, 0.25, 0.5, 0.75, 1]) {
            svg.appendChild(svgEl("line", { x1: X(l), x2: X(l), y1: MT, y2: H - MB, class: "grid" }));
            svg.appendChild(svgEl("text", { x: X(l), y: H - MB + 18, "text-anchor": "middle", class: "tick" }, String(l)));
        }
        svg.appendChild(svgEl("text", { x: (ML + W - MR) / 2, y: H - 6, "text-anchor": "middle", class: "axis-label" }, "λ"));
        svg.appendChild(svgEl("text", { x: 14, y: MT + (H - MT - MB) / 2, "text-anchor": "middle", class: "axis-label",
            transform: `rotate(-90 14 ${MT + (H - MT - MB) / 2})` }, "1/μ"));
        svg.appendChild(svgEl("polyline", { class: "series",
            points: points.map((pt) => `${X(ratNum(pt.lambda))},${Y(ratNum(pt.f))}`).join(" ") }));
        // Every break point attaining the largest or the smallest value is
        // labeled: the extremes can be attained several times (the two
        // local maxima of the parallelograms).
        const fMax = points[+host.dataset.imax].f, fMin = points[+host.dataset.imin].f;
        points.forEach((pt) => {
            svg.appendChild(svgEl("circle", { cx: X(ratNum(pt.lambda)), cy: Y(ratNum(pt.f)), r: 4, class: "marker" }));
            const isMax = ratCmp(pt.f, fMax) === 0, isMin = ratCmp(pt.f, fMin) === 0;
            if ((isMax || isMin) && ratCmp(fMax, fMin) !== 0) {
                const above = isMax;
                const x = X(ratNum(pt.lambda));
                const anchor = x < ML + 40 ? "start" : x > W - MR - 40 ? "end" : "middle";
                svg.appendChild(svgEl("text", { x, y: Y(ratNum(pt.f)) + (above ? -10 : 20),
                    "text-anchor": anchor, class: "direct-label" }, ratStr(pt.f)));
            }
        });

        // Hover: the crosshair snaps to a break point within 8 px, showing
        // exact values there, and shows decimals elsewhere.
        const cross = svgEl("line", { y1: MT, y2: H - MB, class: "crosshair", visibility: "hidden" });
        const dot = svgEl("circle", { r: 5, class: "marker hover-dot", visibility: "hidden" });
        svg.appendChild(cross);
        svg.appendChild(dot);
        const hit = svgEl("rect", { x: ML, y: MT, width: W - ML - MR, height: H - MT - MB, class: "hit" });
        svg.appendChild(hit);
        host.appendChild(svg);
        const tip = document.createElement("div");
        tip.className = "polynorm-tooltip";
        tip.hidden = true;
        host.appendChild(tip);

        function show(evt) {
            const box = svg.getBoundingClientRect();
            const scale = W / box.width;
            const sx = (evt.clientX - box.left) * scale;
            let l = Math.min(1, Math.max(0, (sx - ML) / (W - ML - MR)));
            let html;
            const near = points.find((pt) => Math.abs(X(ratNum(pt.lambda)) - sx) <= 8);
            let fNum;
            if (near) {
                l = ratNum(near.lambda);
                fNum = ratNum(near.f);
                html = `λ = ${escapeHtml(ratStr(near.lambda))}<br>1/μ = ${escapeHtml(ratStr(near.f))}<br>` +
                    `μ = ${escapeHtml(ratStr(ratInv(near.f)))}`;
            } else {
                const P = pieces.find((Q) => ratNum(Q.beta) >= l) || pieces[pieces.length - 1];
                fNum = ratNum(P.s0) + ratNum(P.s1) * l;
                html = `λ ≈ ${l.toFixed(4)}<br>1/μ ≈ ${fNum.toFixed(6)}<br>μ ≈ ${(1 / fNum).toFixed(6)}`;
            }
            cross.setAttribute("x1", X(l)); cross.setAttribute("x2", X(l));
            dot.setAttribute("cx", X(l)); dot.setAttribute("cy", Y(fNum));
            cross.setAttribute("visibility", "visible");
            dot.setAttribute("visibility", "visible");
            tip.innerHTML = html;
            tip.hidden = false;
            const px = X(l) / scale, py = Y(fNum) / scale;
            tip.style.left = (px + 12 + tip.offsetWidth > box.width ? px - 12 - tip.offsetWidth : px + 12) + "px";
            tip.style.top = Math.max(0, py - tip.offsetHeight - 8) + "px";
        }
        function hide() {
            cross.setAttribute("visibility", "hidden");
            dot.setAttribute("visibility", "hidden");
            tip.hidden = true;
        }
        hit.addEventListener("pointermove", show);
        hit.addEventListener("pointerdown", show);
        hit.addEventListener("pointerleave", hide);
    }

    function intInput(el, min, fallback) {
        const v = parseInt(el.value, 10);
        return Number.isInteger(v) && v >= min ? v : fallback;
    }

    // The arguments (inputType, rows, cols, entries, ...extra) of the
    // bindings that take one polytope.
    const singleArgs = (extraArgs) => (polys, cols, root) =>
        [polys[0].inputType, polys[0].rows, cols, polys[0].entries, ...extraArgs(root)];

    // Per mode: the presets, the default preset, the binding and its
    // arguments, the renderer and what to wire after rendering.
    const MODES = {
        densities: {
            presets: DENSITY_PRESETS, initial: "square",
            method: "polytopeNormDensities",
            buildArgs: singleArgs((root) => [root.querySelector(".polynorm-what").value]),
            render: renderResult,
        },
        isMinimalCovering: {
            presets: MINCOV_PRESETS, initial: "rectangle",
            method: "isMinimalCovering",
            buildArgs: singleArgs(() => []),
            render: renderMinimalTest,
        },
        minimalCovering: {
            presets: MINCOV_PRESETS, initial: "rectangle",
            method: "minimalCoveringPolytope",
            buildArgs: singleArgs((root) => {
                const method = root.querySelector(".polynorm-method").value;
                if (method === "grid") {
                    return [method, intInput(root.querySelector(".polynorm-grid-n"), 1, 1),
                        intInput(root.querySelector(".polynorm-rounds"), 1, 3)];
                }
                return [method, intInput(root.querySelector(".polynorm-max-step"), 1, 1000), 0];
            }),
            render: renderMinimalCovering,
        },
        coveringPath: {
            presets: PATH_PRESETS, initial: "cfk",
            method: "coveringPath",
            buildArgs: (polys, cols) => [polys[0].inputType, polys[0].rows, polys[0].entries,
                polys[1].inputType, polys[1].rows, polys[1].entries, cols],
            render: renderCoveringPath,
            afterRender: wirePathChart,
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
            cols: root.querySelector(".copos-cols"),
            run: root.querySelector(".copos-run"),
            stop: root.querySelector(".polynorm-stop"),
            presets: root.querySelector(".copos-presets"),
            status: root.querySelector(".copos-status"),
            result: root.querySelector(".copos-result"),
        };
        // The polytope inputs: one block per .polynorm-polytope, or the
        // section itself when it takes a single polytope. Each block has its
        // input type, number of vectors, file input and grid; the dimension
        // is shared.
        const blockRoots = root.querySelectorAll(".polynorm-polytope");
        const blocks = Array.from(blockRoots.length ? blockRoots : [root], (b) => ({
            name: b.dataset.name || "",
            inputType: b.querySelector(".polynorm-input-type"),
            rows: b.querySelector(".copos-rows"),
            file: b.querySelector(".copos-file"),
            grid: b.querySelector(".copos-grid"),
        }));
        let worker = null;

        function setStatus(msg, kind) {
            els.status.textContent = msg;
            els.status.className = "copos-status" + (kind ? " " + kind : "");
        }

        function clamp(v, max) {
            if (!Number.isFinite(v) || v < 1) return 1;
            return Math.min(v, max);
        }

        function readCols() {
            return clamp(parseInt(els.cols.value, 10), MAX_COLS);
        }

        function readRows(blk) {
            return clamp(parseInt(blk.rows.value, 10), MAX_ROWS);
        }

        function rebuildGrid(blk, rows, cols, entries) {
            blk.grid.innerHTML = "";
            blk.grid.style.gridTemplateColumns = `auto repeat(${cols}, minmax(3rem, 1fr))`;
            for (let i = 0; i < rows; i++) {
                const label = document.createElement("span");
                label.className = "polynorm-row-label";
                label.innerHTML = `v<sub>${i + 1}</sub>`;
                blk.grid.appendChild(label);
                for (let j = 0; j < cols; j++) {
                    const inp = document.createElement("input");
                    inp.type = "text";
                    inp.value = entries ? entries[i * cols + j] : (i === j ? "1" : "0");
                    inp.dataset.i = i;
                    inp.dataset.j = j;
                    inp.setAttribute("aria-label", `${blk.name}v${i + 1}[${j + 1}]`);
                    blk.grid.appendChild(inp);
                }
            }
        }

        function readPolytope(blk, cols) {
            const rows = readRows(blk);
            const entries = new Array(rows * cols);
            blk.grid.querySelectorAll("input").forEach((inp) => {
                entries[+inp.dataset.i * cols + +inp.dataset.j] = inp.value.trim();
            });
            return { inputType: blk.inputType.value, rows, entries };
        }

        // Load one block; its grid takes the current dimension.
        function loadBlock(blk, rows, entries) {
            blk.rows.value = String(rows);
            rebuildGrid(blk, rows, readCols(), entries);
            els.result.innerHTML = "";
        }

        // A preset is one polytope {inputType, rows, cols, entries}, or
        // {cols, polytopes: [...]} for a section of several blocks.
        function applyPreset(p) {
            els.cols.value = String(p.cols);
            (p.polytopes || [p]).forEach((q, k) => {
                blocks[k].inputType.value = q.inputType;
                loadBlock(blocks[k], q.rows, q.entries);
            });
        }

        function finish() {
            els.run.disabled = false;
            els.stop.disabled = true;
        }

        function run() {
            const cols = readCols();
            const polys = blocks.map((blk) => readPolytope(blk, cols));
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
                if (mode.afterRender) mode.afterRender(els.result, msg.result);
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
                args: mode.buildArgs(polys, cols, root),
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
            applyPreset(p);
            setStatus("");
        });

        els.cols.addEventListener("change", () => {
            els.cols.value = String(readCols());
            blocks.forEach((blk) => loadBlock(blk, readRows(blk), null));
            setStatus("");
        });
        for (const blk of blocks) {
            blk.rows.addEventListener("change", () => {
                loadBlock(blk, readRows(blk), null);
                setStatus("");
            });
            blk.file.addEventListener("change", (e) => {
                const f = e.target.files && e.target.files[0];
                e.target.value = "";
                if (!f) return;
                const reader = new FileReader();
                reader.onload = () => {
                    try {
                        const parsed = parseVectorFile(String(reader.result));
                        if (parsed.cols !== readCols()) {
                            // A new dimension resets the other polytopes.
                            els.cols.value = String(parsed.cols);
                            blocks.forEach((b) => { if (b !== blk) loadBlock(b, readRows(b), null); });
                        }
                        loadBlock(blk, parsed.rows, parsed.entries);
                        setStatus(`Loaded ${parsed.rows} vectors in dimension ${parsed.cols} from ${f.name}`, "ok");
                    } catch (err) {
                        setStatus("File error: " + err.message, "error");
                    }
                };
                reader.onerror = () => setStatus("Could not read file", "error");
                reader.readAsText(f);
            });
        }
        els.run.addEventListener("click", run);
        els.stop.addEventListener("click", stop);

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

        applyPreset(PRESETS[mode.initial]);
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
