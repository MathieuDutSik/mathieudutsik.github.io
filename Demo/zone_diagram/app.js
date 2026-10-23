// Interactive zone diagram: drawing and interaction around the wasm sweep.
"use strict";

(function () {
  const ZONE_COLORS = [
    "#2f6db5", "#d1542f", "#2f9a6a", "#a5479c", "#c8961a", "#4aa3b8",
    "#7a5bd1", "#b8425f", "#6b8f2c", "#d27a2e", "#3d7f8f", "#8c6b4f",
  ];
  const KIND_TIP = 2, KIND_CORNER = 3, KIND_CONTACT = 4;

  const canvas = document.getElementById("view");
  const ctx = canvas.getContext("2d");
  const ui = {
    preset: document.getElementById("preset"),
    reload: document.getElementById("reload"),
    clear: document.getElementById("clear"),
    play: document.getElementById("play"),
    level: document.getElementById("level"),
    levelValue: document.getElementById("levelValue"),
    showFacing: document.getElementById("showFacing"),
    showMarks: document.getElementById("showMarks"),
    showVoronoi: document.getElementById("showVoronoi"),
    showShrunk: document.getElementById("showShrunk"),
    alpha: document.getElementById("alpha"),
    alphaValue: document.getElementById("alphaValue"),
    alphaVoronoi: document.getElementById("alphaVoronoi"),
    alphaZone: document.getElementById("alphaZone"),
    levelNote: document.getElementById("levelNote"),
    showChain: document.getElementById("showChain"),
    rays: document.getElementById("rays"),
    skeleton: document.getElementById("skeleton"),
    stats: document.getElementById("stats"),
    check: document.getElementById("check"),
    checkResult: document.getElementById("checkResult"),
    readout: document.getElementById("readout"),
  };

  let Z = null;           // the wasm module
  let sites = [];
  let diagram = null;
  let voronoi = null;
  let view = { cx: 0, cy: 0, scale: 100 };
  let level = Infinity;   // sweep level shown
  let hover = null;       // { x, y } in world coordinates
  let chainPts = null;
  let anim = null;

  const color = (i) => ZONE_COLORS[i % ZONE_COLORS.length];
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  // ---- examples --------------------------------------------------------

  function preset(name) {
    switch (name) {
      case "two": return [[0, 1], [0, -1]];
      case "third": return [[0, 1], [0, -1], [1.7, -1.1]];
      case "collinear": return [[-1, 0], [0, 0], [1, 0]];
      case "pair": return [[0, 0], [2, 0], [2, 0.35]];
      case "hexagon": {
        const out = [[0, 0]];
        for (let k = 0; k < 6; k++) out.push([Math.cos(k * Math.PI / 3), Math.sin(k * Math.PI / 3)]);
        return out;
      }
      default: {
        // A fixed pseudo-random sample, so that the example is reproducible.
        let s = 20261023;
        const rnd = () => ((s = (s * 1103515245 + 12345) % 2147483648) / 2147483648);
        const out = [];
        while (out.length < 12) {
          const p = [rnd() * 2 - 1, rnd() * 2 - 1];
          if (out.every((q) => Math.hypot(p[0] - q[0], p[1] - q[1]) > 0.25)) out.push(p);
        }
        return out;
      }
    }
  }

  // ---- view ------------------------------------------------------------

  // The page may load while hidden, for example in a closed <details> of the
  // host page; the canvas then has no size and fitting waits for one.
  let needsFit = false;

  function resize() {
    const r = canvas.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.max(1, Math.round(r.width * dpr));
    canvas.height = Math.max(1, Math.round(r.height * dpr));
    if (needsFit && r.width > 0 && r.height > 0) {
      fit();
      recompute();
      return;
    }
    draw();
  }

  function fit() {
    const r = canvas.getBoundingClientRect();
    needsFit = !(r.width > 0 && r.height > 0);
    if (needsFit) return;
    if (!sites.length) { view = { cx: 0, cy: 0, scale: Math.min(r.width, r.height) / 4 }; return; }
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of sites) { x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y); }
    const span = Math.max(x1 - x0, y1 - y0, 1);
    view = { cx: (x0 + x1) / 2, cy: (y0 + y1) / 2, scale: Math.min(r.width, r.height) / (span * 2.2) };
  }

  const toScreen = (x, y) => {
    const r = canvas.getBoundingClientRect();
    return [r.width / 2 + (x - view.cx) * view.scale, r.height / 2 - (y - view.cy) * view.scale];
  };
  const toWorld = (sx, sy) => {
    const r = canvas.getBoundingClientRect();
    return [view.cx + (sx - r.width / 2) / view.scale, view.cy - (sy - r.height / 2) / view.scale];
  };

  // ---- computation -----------------------------------------------------

  function recompute(rays) {
    chainPts = null;
    ui.checkResult.textContent = "";
    if (!Z || sites.length === 0) { diagram = null; voronoi = null; draw(); renderPanel(); return; }
    const flat = sites.flat();
    diagram = Z.compute(flat, rays || Number(ui.rays.value), 0, Number(ui.alpha.value));
    const L = diagram.maxLevel * 2;
    voronoi = Z.voronoi(flat, view.cx - L, view.cy - L, view.cx + L, view.cy + L);
    updateChain();
    draw();
    renderPanel();
  }

  // ---- drawing ---------------------------------------------------------

  // Points of the boundary of R_i cut at level s, in angular order. Each
  // point carries the label of the ray it ends, for colouring.
  function zoneOutline(z, s) {
    const pts = [];
    const m = z.rho.length;
    let c = 0;
    const nc = z.corners.length / 3;
    for (let r = 0; r < m; r++) {
      const rho = Math.min(z.rho[r], s);
      const reached = z.rho[r] <= s && z.kind[r] !== 0;
      pts.push({ x: z.x + rho * Math.cos(z.angle[r]), y: z.y + rho * Math.sin(z.angle[r]),
                 facing: reached ? z.facing[r] : -1, kind: reached ? z.kind[r] : 0, ray: r });
      while (c < nc && z.corners[3 * c + 2] === r) {
        const cx = z.corners[3 * c], cy = z.corners[3 * c + 1];
        if (Math.hypot(cx - z.x, cy - z.y) <= s && z.rho[r] <= s && z.rho[(r + 1) % m] <= s) {
          pts.push({ x: cx, y: cy, corner: true });
        }
        c++;
      }
    }
    return pts;
  }

  function draw() {
    const dpr = window.devicePixelRatio || 1;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const r = canvas.getBoundingClientRect();
    ctx.clearRect(0, 0, r.width, r.height);
    ctx.fillStyle = css("--canvas");
    ctx.fillRect(0, 0, r.width, r.height);
    const ink = css("--ink"), muted = css("--muted");
    const alpha = parseFloat(css("--zone-alpha")) || 0.2;

    if (ui.showVoronoi.checked && voronoi) {
      ctx.save();
      ctx.strokeStyle = muted;
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 4]);
      for (const cell of voronoi) {
        ctx.beginPath();
        for (let k = 0; k < cell.length; k += 2) {
          const [sx, sy] = toScreen(cell[k], cell[k + 1]);
          k ? ctx.lineTo(sx, sy) : ctx.moveTo(sx, sy);
        }
        ctx.closePath();
        ctx.stroke();
      }
      ctx.restore();
    }

    if (diagram) {
      const s = Math.min(level, diagram.maxLevel);
      const outlines = diagram.zones.map((z) => zoneOutline(z, s));
      // Fills.
      outlines.forEach((pts, i) => {
        ctx.beginPath();
        pts.forEach((p, k) => { const [sx, sy] = toScreen(p.x, p.y); k ? ctx.lineTo(sx, sy) : ctx.moveTo(sx, sy); });
        ctx.closePath();
        ctx.globalAlpha = alpha;
        ctx.fillStyle = color(i);
        ctx.fill();
        ctx.globalAlpha = 1;
      });
      // Boundaries: each arc in the colour of the zone it faces.
      outlines.forEach((pts, i) => drawBoundary(pts, i, muted));
      if (ui.showShrunk.checked) drawShrunk(outlines, ink);
      if (ui.showMarks.checked) drawMarks(s, ink);
    }

    if (chainPts && ui.showChain.checked) drawChain(ink);

    // Sites.
    sites.forEach(([x, y], i) => {
      const [sx, sy] = toScreen(x, y);
      ctx.beginPath();
      ctx.arc(sx, sy, 5, 0, 2 * Math.PI);
      ctx.fillStyle = color(i);
      ctx.fill();
      ctx.lineWidth = 1.5;
      ctx.strokeStyle = css("--site-stroke");
      ctx.stroke();
      ctx.fillStyle = ink;
      ctx.font = "500 11px " + css("--font-mono");
      ctx.fillText(String(i + 1), sx + 7, sy - 7);
    });
  }

  function drawBoundary(pts, owner, muted) {
    const n = pts.length;
    if (!n) return;
    const facingColour = ui.showFacing.checked;
    // Segment k goes from pts[k] to pts[k+1] and takes the label of the
    // endpoint that is not a corner.
    const pieces = [];
    for (let k = 0; k < n; k++) {
      const a = pts[k], b = pts[(k + 1) % n];
      if (a.corner || b.corner || a.facing === b.facing) {
        pieces.push([a, b, a.corner ? b : a]);
      } else {
        // The facing zone changes between two rays but no corner point was
        // found: split the segment between the two labels.
        const m = { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 };
        pieces.push([a, m, a], [m, b, b]);
      }
    }
    for (const [a, b, lab] of pieces) {
      const open = lab.facing < 0;
      ctx.beginPath();
      const [ax, ay] = toScreen(a.x, a.y), [bx, by] = toScreen(b.x, b.y);
      ctx.moveTo(ax, ay);
      ctx.lineTo(bx, by);
      if (open) {
        ctx.strokeStyle = muted;
        ctx.lineWidth = 1;
        ctx.setLineDash([2, 3]);
      } else {
        ctx.strokeStyle = facingColour ? color(lab.facing) : color(owner);
        ctx.lineWidth = 2.5;
        ctx.setLineDash(lab.kind === KIND_CORNER ? [5, 3] : []);
      }
      ctx.stroke();
    }
    ctx.setLineDash([]);
  }

  // The competitors S_j = p_j + alpha (R_j - p_j).
  function drawShrunk(outlines, ink) {
    const a = diagram.alpha;
    ctx.save();
    ctx.strokeStyle = ink;
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    outlines.forEach((pts, i) => {
      const z = diagram.zones[i];
      ctx.beginPath();
      pts.forEach((p, k) => {
        const [sx, sy] = toScreen(z.x + a * (p.x - z.x), z.y + a * (p.y - z.y));
        k ? ctx.lineTo(sx, sy) : ctx.moveTo(sx, sy);
      });
      ctx.closePath();
      ctx.stroke();
    });
    ctx.restore();
  }

  function drawMarks(s, ink) {
    diagram.zones.forEach((z) => {
      const m = z.rho.length;
      let c = 0;
      for (let r = 0; r < m; r++) {
        if (z.kind[r] !== KIND_CONTACT || z.rho[r] > s) continue;
        const [sx, sy] = toScreen(z.x + z.rho[r] * Math.cos(z.angle[r]), z.y + z.rho[r] * Math.sin(z.angle[r]));
        const [tx, ty] = toScreen(z.contacts[2 * c], z.contacts[2 * c + 1]);
        c++;
        ctx.strokeStyle = ink;
        ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(sx, sy); ctx.lineTo(tx, ty); ctx.stroke();
        ctx.fillStyle = ink;
        ctx.beginPath();
        ctx.moveTo(sx, sy - 5); ctx.lineTo(sx + 5, sy); ctx.lineTo(sx, sy + 5); ctx.lineTo(sx - 5, sy);
        ctx.closePath(); ctx.fill();
      }
      const nc = z.corners.length / 3;
      for (let k = 0; k < nc; k++) {
        const cx = z.corners[3 * k], cy = z.corners[3 * k + 1];
        if (Math.hypot(cx - z.x, cy - z.y) > s) continue;
        const [sx, sy] = toScreen(cx, cy);
        ctx.strokeStyle = ink;
        ctx.lineWidth = 1.5;
        ctx.beginPath(); ctx.arc(sx, sy, 3.5, 0, 2 * Math.PI); ctx.stroke();
      }
    });
  }

  function drawChain(ink) {
    ctx.save();
    ctx.strokeStyle = ink;
    ctx.fillStyle = ink;
    ctx.lineWidth = 1.25;
    const pts = [];
    for (let k = 0; k < chainPts.length; k += 2) pts.push(toScreen(chainPts[k], chainPts[k + 1]));
    for (let k = 0; k + 1 < pts.length; k++) {
      const [ax, ay] = pts[k], [bx, by] = pts[k + 1];
      ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, by); ctx.stroke();
      const ang = Math.atan2(by - ay, bx - ax), len = Math.hypot(bx - ax, by - ay);
      if (len > 14) {
        const mx = ax + (bx - ax) * 0.55, my = ay + (by - ay) * 0.55;
        ctx.beginPath();
        ctx.moveTo(mx + 6 * Math.cos(ang), my + 6 * Math.sin(ang));
        ctx.lineTo(mx - 4 * Math.cos(ang) - 4 * Math.sin(ang), my - 4 * Math.sin(ang) + 4 * Math.cos(ang));
        ctx.lineTo(mx - 4 * Math.cos(ang) + 4 * Math.sin(ang), my - 4 * Math.sin(ang) - 4 * Math.cos(ang));
        ctx.closePath(); ctx.fill();
      }
    }
    pts.forEach(([x, y], k) => { ctx.beginPath(); ctx.arc(x, y, k === 0 ? 4 : 2.5, 0, 2 * Math.PI); ctx.fill(); });
    ctx.restore();
  }

  // ---- dependency chain under the pointer --------------------------------

  // The boundary ray nearest to the pointer, within a few pixels.
  function pickRay(x, y) {
    if (!diagram) return null;
    const s = Math.min(level, diagram.maxLevel);
    let best = null;
    diagram.zones.forEach((z, i) => {
      let a = Math.atan2(y - z.y, x - z.x);
      if (a < 0) a += 2 * Math.PI;
      const ang = z.angle;
      let lo = 0, hi = ang.length;
      while (lo < hi) { const mid = (lo + hi) >> 1; if (ang[mid] < a) lo = mid + 1; else hi = mid; }
      for (const r of [lo % ang.length, (lo + ang.length - 1) % ang.length]) {
        if (z.kind[r] < KIND_TIP || z.rho[r] > s) continue;
        const px = z.x + z.rho[r] * Math.cos(ang[r]), py = z.y + z.rho[r] * Math.sin(ang[r]);
        const d = Math.hypot(px - x, py - y) * view.scale;
        if (d < 10 && (!best || d < best.d)) best = { d, zone: i, ray: r };
      }
    });
    return best;
  }

  function updateChain() {
    chainPts = null;
    // For alpha < 1 influence also runs inwards, so there are no chains.
    if (!hover || !ui.showChain.checked || !Z || !diagram || diagram.alpha < 1) return;
    const pick = pickRay(hover.x, hover.y);
    if (pick) chainPts = Z.chain(pick.zone, pick.ray);
  }

  // ---- side panel ------------------------------------------------------

  function renderPanel() {
    ui.skeleton.textContent = "";
    ui.stats.textContent = "";
    if (!diagram) return;
    diagram.zones.forEach((z, i) => {
      const row = document.createElement("div");
      row.className = "zone";
      const owner = document.createElement("span");
      owner.className = "owner";
      owner.style.color = color(i);
      owner.textContent = "R" + (i + 1);
      row.append(owner);
      for (const a of z.arcs) {
        const chip = document.createElement("span");
        chip.className = "chip" + (a.straight ? " straight" : "") + (a.facing < 0 ? " open" : "");
        if (a.facing >= 0) {
          const sw = document.createElement("i");
          sw.style.background = color(a.facing);
          chip.append(sw);
        }
        chip.append(a.facing < 0 ? "open" : (a.straight ? "corner of R" : "R") + (a.facing + 1) + (a.contact ? " ◆" : ""));
        row.append(chip);
      }
      ui.skeleton.append(row);
    });
    const alpha = diagram.alpha;
    ui.levelNote.textContent = alpha >= 1
      ? "At level s the picture shows R\u1d62 \u2229 B(p\u1d62, s). Everything below s is final: what lies further out never changes it."
      : "At level s the picture shows R\u1d62 \u2229 B(p\u1d62, s). For \u03b1 < 1 a shrunk far part of a zone can reach back, so the sweep is followed by fixed-point passes, and the dependency chains are not shown.";
    const rows = [
      ["sites", sites.length],
      ["α", alpha.toFixed(2)],
      ["fixed-point passes", alpha > 0 && alpha < 1 ? diagram.iterations : "none needed"],
      ["rays per zone", ui.rays.value],
      ["time (wasm)", diagram.ms.toFixed(1) + " ms"],
      ["rays fixed", diagram.pops.toLocaleString("en")],
      ["relaxations", diagram.relaxations.toLocaleString("en")],
      ["contacts", diagram.contacts],
      ["corners", diagram.corners],
      ["sweep stops at", diagram.maxLevel.toFixed(2)],
    ];
    for (const [k, v] of rows) {
      const dt = document.createElement("dt"); dt.textContent = k;
      const dd = document.createElement("dd"); dd.textContent = v;
      ui.stats.append(dt, dd);
    }
  }

  function setLevelFromSlider() {
    const v = Number(ui.level.value);
    if (!diagram || v >= 1000) { level = Infinity; ui.levelValue.textContent = "all"; }
    else { level = (v / 1000) * levelRange(); ui.levelValue.textContent = level.toFixed(3); }
    updateChain();
    draw();
  }

  // The slider runs up to the level where the last bounded ray stops.
  function levelRange() {
    let top = 0;
    for (const z of diagram.zones) for (let r = 0; r < z.rho.length; r++) if (z.kind[r] !== 0) top = Math.max(top, z.rho[r]);
    return Math.max(top * 1.15, 1e-9);
  }

  // ---- interaction -----------------------------------------------------

  let drag = null;

  function siteAt(sx, sy) {
    let best = -1, bd = 12;
    sites.forEach(([x, y], i) => { const [px, py] = toScreen(x, y); const d = Math.hypot(px - sx, py - sy); if (d < bd) { bd = d; best = i; } });
    return best;
  }

  let pending = false;
  function scheduleCompute(rays) {
    if (pending) return;
    pending = true;
    requestAnimationFrame(() => { pending = false; recompute(rays); });
  }

  canvas.addEventListener("pointerdown", (e) => {
    const r = canvas.getBoundingClientRect();
    const sx = e.clientX - r.left, sy = e.clientY - r.top;
    const i = siteAt(sx, sy);
    canvas.setPointerCapture(e.pointerId);
    drag = { site: i, sx, sy, moved: false, cx: view.cx, cy: view.cy };
  });
  canvas.addEventListener("pointermove", (e) => {
    const r = canvas.getBoundingClientRect();
    const sx = e.clientX - r.left, sy = e.clientY - r.top;
    const [wx, wy] = toWorld(sx, sy);
    ui.readout.textContent = "x " + wx.toFixed(3) + "   y " + wy.toFixed(3);
    if (drag) {
      if (Math.hypot(sx - drag.sx, sy - drag.sy) > 3) drag.moved = true;
      if (!drag.moved) return;
      if (drag.site >= 0) {
        sites[drag.site] = [wx, wy];
        scheduleCompute(Math.min(360, Number(ui.rays.value)));
      } else {
        view.cx = drag.cx - (sx - drag.sx) / view.scale;
        view.cy = drag.cy + (sy - drag.sy) / view.scale;
        draw();
      }
      return;
    }
    hover = { x: wx, y: wy };
    updateChain();
    draw();
  });
  canvas.addEventListener("pointerup", (e) => {
    const r = canvas.getBoundingClientRect();
    const sx = e.clientX - r.left, sy = e.clientY - r.top;
    if (drag && !drag.moved && drag.site < 0) {
      sites.push(toWorld(sx, sy));
      recompute();
    } else if (drag && drag.moved && drag.site >= 0) {
      recompute();
    }
    drag = null;
  });
  canvas.addEventListener("pointerleave", () => { hover = null; chainPts = null; draw(); });
  canvas.addEventListener("dblclick", (e) => {
    const r = canvas.getBoundingClientRect();
    const i = siteAt(e.clientX - r.left, e.clientY - r.top);
    if (i >= 0) { sites.splice(i, 1); recompute(); }
  });
  canvas.addEventListener("wheel", (e) => {
    e.preventDefault();
    const r = canvas.getBoundingClientRect();
    const sx = e.clientX - r.left, sy = e.clientY - r.top;
    const [wx, wy] = toWorld(sx, sy);
    view.scale *= Math.exp(-e.deltaY * 0.0015);
    const [nx, ny] = toWorld(sx, sy);
    view.cx += wx - nx;
    view.cy += wy - ny;
    draw();
  }, { passive: false });

  function load(name) {
    sites = preset(name);
    fit();
    recompute();
  }

  ui.reload.addEventListener("click", () => load(ui.preset.value));
  ui.preset.addEventListener("change", () => load(ui.preset.value));
  ui.clear.addEventListener("click", () => { sites = []; recompute(); });
  ui.rays.addEventListener("change", () => recompute());
  const showAlpha = () => { ui.alphaValue.textContent = Number(ui.alpha.value).toFixed(2); };
  ui.alpha.addEventListener("input", () => { showAlpha(); scheduleCompute(Math.min(360, Number(ui.rays.value))); });
  ui.alpha.addEventListener("change", () => { showAlpha(); recompute(); });
  for (const [button, value] of [[ui.alphaVoronoi, 0], [ui.alphaZone, 1]]) {
    button.addEventListener("click", () => { ui.alpha.value = value; showAlpha(); recompute(); });
  }
  ui.level.addEventListener("input", () => { stopAnim(); setLevelFromSlider(); });
  for (const box of [ui.showFacing, ui.showMarks, ui.showVoronoi, ui.showShrunk, ui.showChain]) {
    box.addEventListener("change", () => { updateChain(); draw(); });
  }
  ui.check.addEventListener("click", () => {
    if (!diagram) return;
    const res = Z.residual(1);
    ui.checkResult.textContent = "residual " + res.toExponential(2);
  });

  function stopAnim() {
    if (anim) cancelAnimationFrame(anim);
    anim = null;
    ui.play.textContent = "Grow";
  }
  ui.play.addEventListener("click", () => {
    if (anim) { stopAnim(); return; }
    if (!diagram) return;
    const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    if (reduced) { ui.level.value = 1000; setLevelFromSlider(); return; }
    ui.play.textContent = "Pause";
    let v = Number(ui.level.value) >= 1000 ? 0 : Number(ui.level.value);
    let last = performance.now();
    const step = (t) => {
      v = Math.min(1000, v + (t - last) * 0.12);
      last = t;
      ui.level.value = v;
      setLevelFromSlider();
      if (v >= 1000) { stopAnim(); return; }
      anim = requestAnimationFrame(step);
    };
    anim = requestAnimationFrame(step);
  });

  new ResizeObserver(resize).observe(canvas);
  new MutationObserver(draw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);

  createZoneModule().then((module) => {
    Z = module;
    showAlpha();
    resize();
    load(ui.preset.value);
  });
})();
