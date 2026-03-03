// generative-model-viz.js
//
// 3D "lift" visualization of the hidden manifold model:
//   — latent points c^μ sit on a ghost plane at the bottom (z = Z_FLOOR)
//   — teacher direction θ⁰ is drawn as an arrow on that plane
//   — projected/observed points x^μ float above in 3D (after σ(Fc))
//   — thin lift lines connect each latent point to its observed copy
//   — everything colored by label value with a custom muted-pastel scale

// ── Palette ──────────────────────────────────────────────────────────────────

// Purple-to-teal diverging scale (matches scientific software post palette)
// negative label → light purple, zero → neutral, positive → teal
const PASTEL_SCALE = [
  [0.00, "#e8c8f0"],   // light purple
  [0.25, "#c49bd4"],   // medium purple
  [0.50, "#f5f5f5"],   // neutral
  [0.75, "#a0d4d8"],   // light teal
  [1.00, "#0CBABA"],   // teal
];

const THETA_COLOR  = "#861388";   // purple — main accent
const PLANE_COLOR  = "#e0e0e0";   // clean light gray
const LIFT_COLOR   = "#aaaaaa";   // neutral gray
const FONT         = "'Inter', system-ui, sans-serif";
const BG           = "#ffffff";   // clean white
const AXIS_COLOR   = "#666666";

// ── Math helpers ──────────────────────────────────────────────────────────────

function rng() {
  let u, v;
  do { u = Math.random(); } while (!u);
  do { v = Math.random(); } while (!v);
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function applyNl(x, type) {
  if (type === "tanh")  return Math.tanh(x);
  if (type === "relu")  return Math.max(0, x);
  if (type === "sign")  return x >= 0 ? 1 : -1;
  return x;  // identity
}

// ── Main export ───────────────────────────────────────────────────────────────

export function generativeModelViz({ sigmaType, f0Type, thetaAngle, rerollTrigger, d3, Plotly }) {

  const N = 90;   // number of points — dense enough to read the cloud

  // Ground-truth direction θ (2D unit vector)
  const theta = [Math.cos(thetaAngle), Math.sin(thetaAngle)];

  // Random projection matrix F  (2 × 3, each column is a 3D direction)
  const F = [
    [rng(), rng(), rng()],
    [rng(), rng(), rng()],
  ];

  // Generate data
  const pts = Array.from({ length: N }, () => {
    const c = [rng(), rng()];
    const nu = (c[0] * theta[0] + c[1] * theta[1]);   // c · θ  (d=2, so no 1/√d needed for vis)
    const y  = applyNl(nu, f0Type);
    // x = σ( F^T c )  →  3D vector
    const xRaw = [
      F[0][0] * c[0] + F[1][0] * c[1],
      F[0][1] * c[0] + F[1][1] * c[1],
      F[0][2] * c[0] + F[1][2] * c[1],
    ];
    const x = xRaw.map(v => applyNl(v, sigmaType));
    return { c, x, y };
  });

  // ── Z floor for latent plane ──────────────────────────────────────────────
  const allZ = pts.map(p => p.x[2]);
  const zMin = Math.min(...allZ);
  const zMax = Math.max(...allZ);
  const Z_FLOOR = zMin - 1.4;   // a little below the lowest projected point

  // ── Plotly traces ─────────────────────────────────────────────────────────
  const traces = [];

  // 1. Ghost plane (latent space floor) — mesh3d
  const PS = 3.2;  // half-size
  traces.push({
    type: "mesh3d",
    x: [-PS,  PS,  PS, -PS],
    y: [-PS, -PS,  PS,  PS],
    z: [Z_FLOOR, Z_FLOOR, Z_FLOOR, Z_FLOOR],
    i: [0, 0], j: [1, 2], k: [2, 3],
    color: PLANE_COLOR,
    opacity: 0.18,
    showlegend: false,
    hoverinfo: "skip",
    flatshading: true,
    lighting: { diffuse: 0.5, ambient: 1 },
  });

  // 2. Lift lines: latent → projected  (one scatter3d per point, hidden from legend)
  pts.forEach(p => {
    traces.push({
      type: "scatter3d", mode: "lines",
      x: [p.c[0], p.x[0]],
      y: [p.c[1], p.x[1]],
      z: [Z_FLOOR, p.x[2]],
      line: { color: LIFT_COLOR, width: 1 },
      opacity: 0.35,
      showlegend: false,
      hoverinfo: "skip",
    });
  });

  // 3. Latent points on the floor plane — colored by label
  traces.push({
    type: "scatter3d", mode: "markers",
    x: pts.map(p => p.c[0]),
    y: pts.map(p => p.c[1]),
    z: pts.map(() => Z_FLOOR),
    marker: {
      size: 4.5,
      color: pts.map(p => p.y),
      colorscale: PASTEL_SCALE,
      line: { color: "white", width: 0.5 },
      opacity: 1,
    },
    name: "c^μ  — latent vectors (ℝᵈ)",
    hovertemplate: "c^μ<br>label = %{marker.color:.2f}<extra></extra>",
  });

  // 4. Teacher direction θ⁰ — arrow on the latent plane
  const ARROW_LEN = 1.8;
  traces.push({
    type: "scatter3d", mode: "lines+markers",
    x: [0, theta[0] * ARROW_LEN],
    y: [0, theta[1] * ARROW_LEN],
    z: [Z_FLOOR, Z_FLOOR],
    line: { color: THETA_COLOR, width: 4 },
    marker: {
      size: [0, 9],
      color: THETA_COLOR,
      symbol: ["circle", "diamond"],
    },
    name: "θ⁰  (teacher direction)",
    hoverinfo: "skip",
  });

  // 4b. θ⁰ text label at arrowhead — placed as a text trace so it scales with the scene
  traces.push({
    type: "scatter3d", mode: "text",
    x: [theta[0] * ARROW_LEN * 1.28],
    y: [theta[1] * ARROW_LEN * 1.28],
    z: [Z_FLOOR],
    text: ["θ⁰"],
    textfont: { size: 13, color: THETA_COLOR, family: FONT },
    textposition: "middle center",
    showlegend: false,
    hoverinfo: "skip",
  });

  // 5. Projected / observed points in 3D — same colorscale, slightly transparent
  traces.push({
    type: "scatter3d", mode: "markers",
    x: pts.map(p => p.x[0]),
    y: pts.map(p => p.x[1]),
    z: pts.map(p => p.x[2]),
    marker: {
      size: 5.5,
      color: pts.map(p => p.y),
      colorscale: PASTEL_SCALE,
      cmin: Math.min(...pts.map(p => p.y)),
      cmax: Math.max(...pts.map(p => p.y)),
      colorbar: {
        title: { text: "label", font: { size: 11, family: FONT, color: AXIS_COLOR } },
        tickfont: { size: 10, family: FONT, color: AXIS_COLOR },
        thickness: 10,
        len: 0.65,
        x: 1.01,
        outlinewidth: 0,
      },
      line: { color: "white", width: 0.6 },
      opacity: 0.88,
    },
    name: "x^μ = σ(Fc^μ)  — observed (ℝᵖ)",
    hovertemplate: "x^μ<br>label = %{marker.color:.2f}<extra></extra>",
  });

  // ── Layout ────────────────────────────────────────────────────────────────
  const axisStyle = (title) => ({
    title: { text: title, font: { size: 11, color: AXIS_COLOR, family: FONT } },
    tickfont: { size: 9, color: AXIS_COLOR, family: FONT },
    showgrid: false,
    showline: false,
    zeroline: false,
    showticklabels: false,   // coordinates are arbitrary — hide for clarity
    backgroundcolor: BG,
    gridcolor: "#ddd",
  });

  const layout = {
    paper_bgcolor: BG,
    font: { family: FONT },
    margin: { l: 0, r: 60, t: 16, b: 0 },
    scene: {
      xaxis: axisStyle("x₁"),
      yaxis: axisStyle("x₂"),
      zaxis: axisStyle("x₃"),
      bgcolor: BG,
      camera: { eye: { x: 1.6, y: -1.8, z: 1.2 } },
      aspectmode: "cube",
      annotations: [
        // ── Latent plane label (back-left corner of the ghost plane) ──────
        {
          x: -PS * 0.72, y: PS * 0.72, z: Z_FLOOR,
          text: "Latent space ℝ<sup>d</sup>",
          showarrow: false,
          font: { size: 11, color: AXIS_COLOR, family: FONT },
          bgcolor: "rgba(255,255,255,0.92)",
          bordercolor: PLANE_COLOR,
          borderwidth: 1,
          borderpad: 4,
          xanchor: "center",
        },
        // ── Observed cloud label (above the cloud centroid) ───────────────
        {
          x: 0, y: 0, z: zMax + 0.9,
          text: "Observed space ℝ<sup>p</sup>",
          showarrow: false,
          font: { size: 11, color: AXIS_COLOR, family: FONT },
          bgcolor: "rgba(255,255,255,0.92)",
          bordercolor: "#ddd",
          borderwidth: 1,
          borderpad: 4,
          xanchor: "center",
        },
        // ── Lift-line label (labels the σ(Fc) mapping) ────────────────────
        {
          // position near the midpoint between one latent & observed point
          x: pts[0].c[0] * 0.5 + pts[0].x[0] * 0.5,
          y: pts[0].c[1] * 0.5 + pts[0].x[1] * 0.5,
          z: Z_FLOOR + (pts[0].x[2] - Z_FLOOR) * 0.5,
          text: "σ(Fc)",
          showarrow: true,
          ax: 30, ay: -20,
          font: { size: 10, color: LIFT_COLOR, family: FONT },
          bgcolor: "rgba(0,0,0,0)",
          borderwidth: 0,
          arrowcolor: LIFT_COLOR,
          arrowsize: 0.8,
          arrowwidth: 1,
        },
      ],
    },
    legend: {
      x: 0.02, y: 0.98,
      bgcolor: "rgba(255,255,255,0.90)",
      bordercolor: "#ddd", borderwidth: 1,
      font: { size: 11, family: FONT, color: "#555" },
    },
    hovermode: "closest",
  };

  const config = {
    responsive: true,
    displayModeBar: false,   // no toolbar
    scrollZoom: false,
  };

  // ── DOM assembly ──────────────────────────────────────────────────────────
  const wrap = d3.create("div")
    .style("width", "100%")
    .style("height", "560px")
    .style("border-radius", "8px")
    .style("overflow", "hidden")
    .style("background", BG);

  Plotly.newPlot(wrap.node(), traces, layout, config);

  return wrap.node();
}
