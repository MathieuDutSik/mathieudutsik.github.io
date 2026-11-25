// Web Worker running the packing / covering computation of the Tools page
// (assets/polynorm-ui.js). The covering search can take a long time, so it
// runs here rather than on the page's main thread: the page stays responsive
// and "Stop" terminates the worker.
//
// Message in : { bust, inputType, rows, cols, entries, what }
// Message out: { ok: true, result, elapsedMs } or { ok: false, error }
// The result is converted to plain arrays, since embind vectors cannot be
// posted between threads.

"use strict";

let modulePromise = null;

function loadModule(bust) {
    if (modulePromise) return modulePromise;
    importScripts("/wasm/polyhedral.js" + bust);
    modulePromise = createPolyhedralModule({
        locateFile: (path) => path.endsWith(".wasm") ? "/wasm/" + path + bust : path,
    });
    return modulePromise;
}

function toArray(vec) {
    const a = [];
    for (let i = 0; i < vec.size(); i++) a.push(vec.get(i));
    vec.delete();
    return a;
}

self.onmessage = async (e) => {
    const msg = e.data;
    try {
        const m = await loadModule(msg.bust);
        const t0 = performance.now();
        const r = m.polytopeNormDensities(msg.inputType, msg.rows, msg.cols,
            msg.entries, msg.what);
        const elapsedMs = performance.now() - t0;
        const result = {};
        for (const key of Object.keys(r)) {
            const v = r[key];
            result[key] = (v && typeof v.size === "function") ? toArray(v) : v;
        }
        self.postMessage({ ok: true, result, elapsedMs });
    } catch (err) {
        self.postMessage({ ok: false, error: err && err.message ? err.message : String(err) });
    }
};
