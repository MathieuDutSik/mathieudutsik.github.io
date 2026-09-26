import { useState, useMemo, useRef, useCallback, useEffect } from "react";

function dist(a, b) {
  const dx = a.x - b.x;
  const dy = a.y - b.y;
  return Math.sqrt(dx * dx + dy * dy);
}

// Seeded PRNG (mulberry32) for deterministic results at a given N + vectors
function mulberry32(seed) {
  return () => {
    seed |= 0;
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function App() {
  const svgRef = useRef(null);
  const [radius, setRadius] = useState(6);
  const [sampleCount, setSampleCount] = useState(1000);
  const [k, setK] = useState(1);
  const [v1End, setV1End] = useState({ x: 80, y: 0 });
  const [v2End, setV2End] = useState({ x: 0, y: -80 });
  const [dragging, setDragging] = useState(null);
  const [svgSize, setSvgSize] = useState({ w: 800, h: 600 });

  const containerRef = useRef(null);

  useEffect(() => {
    const update = () => {
      if (containerRef.current) {
        const rect = containerRef.current.getBoundingClientRect();
        setSvgSize({ w: rect.width, h: rect.height });
      }
    };
    update();
    window.addEventListener("resize", update);
    return () => window.removeEventListener("resize", update);
  }, []);

  const cx = svgSize.w / 2;
  const cy = svgSize.h / 2;

  const toSvg = useCallback(
    (e) => {
      const svg = svgRef.current;
      const rect = svg.getBoundingClientRect();
      return {
        x: e.clientX - rect.left - cx,
        y: e.clientY - rect.top - cy,
      };
    },
    [cx, cy],
  );

  const onMouseDown = useCallback(
    (which) => (e) => {
      e.preventDefault();
      setDragging(which);
    },
    [],
  );

  const onMouseMove = useCallback(
    (e) => {
      if (!dragging) return;
      const pos = toSvg(e);
      if (dragging === "v1") setV1End(pos);
      else setV2End(pos);
    },
    [dragging, toSvg],
  );

  const onMouseUp = useCallback(() => {
    setDragging(null);
  }, []);

  // All lattice points (including origin) for nearest-neighbor search
  const allLattice = useMemo(() => {
    const pts = [];
    for (let n1 = -5; n1 <= 5; n1++) {
      for (let n2 = -5; n2 <= 5; n2++) {
        pts.push({
          x: n1 * v1End.x + n2 * v2End.x,
          y: n1 * v1End.y + n2 * v2End.y,
        });
      }
    }
    return pts;
  }, [v1End, v2End]);

  // Compute k-covering radius r (distance to k-th nearest lattice point)
  const coveringRadius = useMemo(() => {
    const rng = mulberry32(42);
    let maxDist = 0;
    for (let i = 0; i < sampleCount; i++) {
      const x1 = rng();
      const x2 = rng();
      const px = x1 * v1End.x + x2 * v2End.x;
      const py = x1 * v1End.y + x2 * v2End.y;
      const pt = { x: px, y: py };
      // Collect distances to all lattice points and pick k-th smallest
      const dists = allLattice.map((lp) => dist(pt, lp));
      dists.sort((a, b) => a - b);
      const kthDist = dists[k - 1] ?? dists[dists.length - 1];
      if (kthDist > maxDist) maxDist = kthDist;
    }
    return maxDist;
  }, [sampleCount, k, v1End, v2End, allLattice]);

  // Lattice points visible in the rectangle (excluding origin)
  const latticePoints = [];
  for (const lp of allLattice) {
    if (lp.x === 0 && lp.y === 0) continue;
    const margin = radius;
    if (
      lp.x >= -cx + margin &&
      lp.x <= cx - margin &&
      lp.y >= -cy + margin &&
      lp.y <= cy - margin
    ) {
      latticePoints.push(lp);
    }
  }

  // All lattice points visible (including origin) for drawing covering circles
  const visibleAllLattice = allLattice.filter((lp) => {
    const margin = -coveringRadius; // allow circles that partially overlap the viewport
    return (
      lp.x >= -cx + margin &&
      lp.x <= cx - margin &&
      lp.y >= -cy + margin &&
      lp.y <= cy - margin
    );
  });

  return (
    <div className="app">
      {/* Controls */}
      <div className="toolbar">
        <label className="control">
          Point radius:
          <input
            type="range"
            min={2}
            max={20}
            value={radius}
            onChange={(e) => setRadius(Number(e.target.value))}
            className="range"
          />
          <span className="range-value">{radius}</span>
        </label>
        <label className="control">
          N:
          <input
            type="text"
            value={sampleCount}
            onChange={(e) => {
              const v = parseInt(e.target.value, 10);
              if (!isNaN(v) && v > 0) setSampleCount(v);
            }}
            className="field field-wide"
          />
        </label>
        <label className="control">
          k:
          <input
            type="number"
            min={1}
            value={k}
            onChange={(e) => {
              const v = parseInt(e.target.value, 10);
              if (!isNaN(v) && v >= 1) setK(v);
            }}
            className="field"
          />
        </label>
        <span className="note">
          {k}-covering radius r = {coveringRadius.toFixed(1)}px
        </span>
        <span className="note hint">
          Drag the arrow tips to move v1 / v2
        </span>
      </div>

      {/* SVG canvas */}
      <div ref={containerRef} className="canvas">
        <svg
          ref={svgRef}
          width={svgSize.w}
          height={svgSize.h}
          className="plane"
          onMouseMove={onMouseMove}
          onMouseUp={onMouseUp}
          onMouseLeave={onMouseUp}
        >
          <defs>
            <marker
              id="arrow-v1"
              markerWidth="10"
              markerHeight="7"
              refX="9"
              refY="3.5"
              orient="auto"
            >
              <polygon points="0 0, 10 3.5, 0 7" fill="#22c55e" />
            </marker>
            <marker
              id="arrow-v2"
              markerWidth="10"
              markerHeight="7"
              refX="9"
              refY="3.5"
              orient="auto"
            >
              <polygon points="0 0, 10 3.5, 0 7" fill="#a855f7" />
            </marker>
          </defs>

          {/* Covering circles (black, drawn first so they're behind everything) */}
          {visibleAllLattice.map((pt, i) => (
            <circle
              key={`cov-${i}`}
              cx={cx + pt.x}
              cy={cy + pt.y}
              r={coveringRadius}
              fill="none"
              stroke="red"
              strokeWidth={1.5}
            />
          ))}

          {/* Lattice points (blue) */}
          {latticePoints.map((pt, i) => (
            <circle
              key={i}
              cx={cx + pt.x}
              cy={cy + pt.y}
              r={radius}
              fill="#3b82f6"
              opacity={0.8}
            />
          ))}

          {/* Vector v1 (green) */}
          <line
            x1={cx}
            y1={cy}
            x2={cx + v1End.x}
            y2={cy + v1End.y}
            stroke="#22c55e"
            strokeWidth={2}
            markerEnd="url(#arrow-v1)"
          />
          <text
            x={cx + v1End.x + 10}
            y={cy + v1End.y - 10}
            fill="#22c55e"
            fontSize={14}
            fontWeight="bold"
          >
            v1
          </text>
          <circle
            cx={cx + v1End.x}
            cy={cy + v1End.y}
            r={8}
            fill="#22c55e"
            stroke="white"
            strokeWidth={2}
            className="handle"
            onMouseDown={onMouseDown("v1")}
          />

          {/* Vector v2 (purple) */}
          <line
            x1={cx}
            y1={cy}
            x2={cx + v2End.x}
            y2={cy + v2End.y}
            stroke="#a855f7"
            strokeWidth={2}
            markerEnd="url(#arrow-v2)"
          />
          <text
            x={cx + v2End.x + 10}
            y={cy + v2End.y - 10}
            fill="#a855f7"
            fontSize={14}
            fontWeight="bold"
          >
            v2
          </text>
          <circle
            cx={cx + v2End.x}
            cy={cy + v2End.y}
            r={8}
            fill="#a855f7"
            stroke="white"
            strokeWidth={2}
            className="handle"
            onMouseDown={onMouseDown("v2")}
          />

          {/* Center point p (red, always on top) */}
          <circle cx={cx} cy={cy} r={radius} fill="#ef4444" />
          <text
            x={cx + radius + 4}
            y={cy - radius - 4}
            fill="#ef4444"
            fontSize={14}
            fontWeight="bold"
          >
            p
          </text>
        </svg>
      </div>
    </div>
  );
}

export default App;
