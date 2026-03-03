// interactive-plots.js
// Clean, minimal interactive generalization error explorer.
// Controls: λ slider, nonlinearity selector, task selector, Δ slider (regression).
// Plots:    (1) gen error vs p/n  (2) 3D surface n/d × p/n × error

// ── Palette (matches scientific software post) ────────────────────────────────
const P = {
  bg:       "#ffffff",   // clean white
  panelBg:  "#f8f8f8",   // very light gray for control panel
  grid:     "#e8e8e8",
  border:   "#ddd",
  text:     "#333",
  muted:    "#666",
  blue:     "#2F6690",   // teal-blue (software post H3 color)
  spike:    "#861388",   // purple — marks the double-descent spike
  font:     "'Inter', system-ui, sans-serif",
};

// Purple-to-teal colorscale for the 3D surface
const SURFACE_SCALE = [
  [0.00, "#e8c8f0"],   // light purple
  [0.25, "#c49bd4"],   // medium purple
  [0.50, "#f5f5f5"],   // neutral
  [0.75, "#a0d4d8"],   // light teal
  [1.00, "#0CBABA"],   // teal
];

// ── Minimal Plotly config ────────────────────────────────────────────────────
const CFG_2D = {
  responsive: true,
  displayModeBar: false,   // clean — no toolbar
};

const CFG_3D = {
  responsive: true,
  displayModeBar: true,
  displaylogo: false,
  modeBarButtonsToRemove: [
    "zoom3d", "pan3d", "orbitRotation", "tableRotation",
    "handleDrag3d", "resetCameraDefault3d", "resetCameraLastSave3d",
    "hoverClosest3d",
  ],
};

// ── Shared Plotly axis style ──────────────────────────────────────────────────
function axStyle(title) {
  return {
    title: { text: title, font: { size: 12, color: P.muted, family: P.font } },
    tickfont: { size: 10, color: P.muted, family: P.font },
    gridcolor: P.grid,
    linecolor: P.border,
    zerolinecolor: P.grid,
    showgrid: true,
    zeroline: false,
  };
}

// ── CSV loader ────────────────────────────────────────────────────────────────
async function loadCSV(filename) {
  const res = await fetch(`data/${filename}`);
  if (!res.ok) throw new Error(`Cannot load ${filename}`);
  const text = await res.text();
  const [header, ...rows] = text.trim().split("\n");
  const keys = header.split(",");
  return rows.map(line => {
    const vals = line.split(",");
    const row = {};
    keys.forEach((k, i) => {
      row[k] = isNaN(vals[i]) ? vals[i] : parseFloat(vals[i]);
    });
    return row;
  });
}

// ── CSS injection (one-time) ──────────────────────────────────────────────────
function injectStyles() {
  if (document.getElementById("gen-err-styles")) return;
  const el = document.createElement("style");
  el.id = "gen-err-styles";
  el.textContent = `
    .gen-ctrl {
      display: flex;
      flex-wrap: wrap;
      gap: 24px;
      align-items: flex-end;
      padding: 14px 20px;
      background: ${P.panelBg};
      border: 1px solid ${P.border};
      border-radius: 8px;
      margin-bottom: 20px;
      font-family: ${P.font};
    }
    .gen-ctrl-group {
      display: flex;
      flex-direction: column;
      gap: 6px;
      min-width: 120px;
    }
    .gen-ctrl label {
      font-size: 11px;
      font-weight: 600;
      letter-spacing: 0.04em;
      text-transform: uppercase;
      color: ${P.muted};
    }
    .gen-ctrl .val-badge {
      font-size: 12px;
      font-weight: 500;
      color: ${P.blue};
      margin-top: -2px;
    }
    .gen-ctrl input[type=range] {
      -webkit-appearance: none;
      width: 140px;
      height: 4px;
      border-radius: 2px;
      background: ${P.grid};
      outline: none;
      cursor: pointer;
    }
    .gen-ctrl input[type=range]::-webkit-slider-thumb {
      -webkit-appearance: none;
      width: 14px; height: 14px;
      border-radius: 50%;
      background: ${P.blue};
      border: 2px solid ${P.bg};
      box-shadow: 0 1px 4px rgba(0,0,0,.12);
      cursor: pointer;
    }
    .gen-ctrl select {
      padding: 5px 10px;
      border: 1px solid ${P.border};
      border-radius: 6px;
      background: white;
      font-size: 12px;
      color: ${P.text};
      font-family: ${P.font};
      cursor: pointer;
      outline: none;
    }
    .gen-ctrl select:focus { border-color: ${P.blue}; }
    .gen-plot-wrap {
      border: 1px solid ${P.border};
      border-radius: 8px;
      overflow: hidden;
      background: white;
    }
  `;
  document.head.appendChild(el);
}

// ── Main export ───────────────────────────────────────────────────────────────
export function createInteractivePlots({ d3, Plotly }) {
  injectStyles();

  // ── Lambda + Delta value maps ────────────────────────────────────────────
  const LAMBDAS = [1e-4, 1e-3, 1e-2, 1e-1, 1.0];
  const DELTAS  = [0.0, 0.1, 0.5, 1.0];

  // ── Build UI ─────────────────────────────────────────────────────────────
  const root = d3.create("div").style("font-family", P.font).style("max-width", "900px").style("margin", "0 auto");

  // Control panel
  const ctrl = root.append("div").attr("class", "gen-ctrl");

  // λ slider
  const gLambda = ctrl.append("div").attr("class", "gen-ctrl-group");
  gLambda.append("label").text("Regularization λ");
  const lambdaBadge = gLambda.append("div").attr("class", "val-badge").text(`λ = ${LAMBDAS[1]}`);
  const lambdaIn = gLambda.append("input")
    .attr("type", "range").attr("min", 0).attr("max", 4).attr("step", 1).attr("value", 1);

  // Nonlinearity selector
  const gNl = ctrl.append("div").attr("class", "gen-ctrl-group");
  gNl.append("label").text("Nonlinearity σ");
  const nlSel = gNl.append("select");
  nlSel.selectAll("option").data(["sign", "erf"]).enter().append("option")
    .attr("value", d => d).text(d => d === "sign" ? "Sign" : "Erf");

  // Task selector
  const gTask = ctrl.append("div").attr("class", "gen-ctrl-group");
  gTask.append("label").text("Task");
  const taskSel = gTask.append("select");
  taskSel.selectAll("option").data(["classification", "regression"]).enter().append("option")
    .attr("value", d => d).text(d => d[0].toUpperCase() + d.slice(1));

  // Δ (noise) slider — only shown for regression
  const gDelta = ctrl.append("div").attr("class", "gen-ctrl-group");
  gDelta.append("label").text("Noise Δ");
  const deltaBadge = gDelta.append("div").attr("class", "val-badge").text("Δ = 0.0");
  const deltaIn = gDelta.append("input")
    .attr("type", "range").attr("min", 0).attr("max", 3).attr("step", 1).attr("value", 0);

  // ── Plot containers ──────────────────────────────────────────────────────
  const wrap1d = root.append("div").attr("class", "gen-plot-wrap")
    .style("height", "380px").style("margin-bottom", "16px");
  wrap1d.node().id = "gep-1d";

  const wrap3d = root.append("div").attr("class", "gen-plot-wrap")
    .style("height", "460px");
  wrap3d.node().id = "gep-3d";

  // ── State helpers ────────────────────────────────────────────────────────
  const getLambda = () => LAMBDAS[parseInt(lambdaIn.property("value"))];
  const getDelta  = () => DELTAS[parseInt(deltaIn.property("value"))];

  function updateBadges() {
    lambdaBadge.text(`λ = ${getLambda()}`);
    deltaBadge.text(`Δ = ${getDelta()}`);
  }

  function toggleDelta() {
    const show = taskSel.property("value") === "regression";
    gDelta.style("display", show ? "flex" : "none");
  }

  // ── Plot: 1D generalization error vs p/n ────────────────────────────────
  async function render1D() {
    const task = taskSel.property("value");
    const nl   = nlSel.property("value");
    const lam  = getLambda();
    const del  = getDelta();

    let rows;
    try {
      if (task === "classification") {
        rows = (await loadCSV("classification_1d.csv"))
          .filter(r => r.nonlinearity === nl && Math.abs(r.lambda - lam) < 1e-12);
      } else {
        rows = (await loadCSV("regression_1d.csv"))
          .filter(r =>
            r.nonlinearity === nl &&
            Math.abs(r.lambda - lam) < 1e-12 &&
            Math.abs(r.Delta - del) < 1e-12
          );
      }
    } catch { rows = []; }

    if (!rows.length) return;
    rows.sort((a, b) => a.p_n - b.p_n);

    // Find spike if λ is small
    const annotations = [];
    if (lam <= 1e-2) {
      const near = rows.filter(r => r.p_n >= 0.8 && r.p_n <= 1.2);
      if (near.length) {
        const pk = near.reduce((m, r) => r.gen_error > m.gen_error ? r : m);
        annotations.push({
          x: pk.p_n, y: pk.gen_error,
          text: "interpolation<br>threshold",
          showarrow: true, arrowhead: 2, arrowwidth: 1.5,
          arrowcolor: P.spike, ax: 28, ay: -34,
          font: { color: P.spike, size: 11, family: P.font },
        });
      }
    }

    const trace = {
      x: rows.map(r => r.p_n),
      y: rows.map(r => r.gen_error),
      type: "scatter", mode: "lines",
      line: { color: P.blue, width: 2.5 },
      hovertemplate: "p/n = %{x:.2f}<br>error = %{y:.4f}<extra></extra>",
    };

    const layout = {
      paper_bgcolor: P.bg,
      plot_bgcolor:  P.bg,
      font: { family: P.font, color: P.text },
      margin: { l: 56, r: 24, t: 20, b: 52 },
      xaxis: { ...axStyle("p / n  (features per sample)"), range: [0.25, 2.6] },
      yaxis: axStyle("Generalization Error"),
      annotations,
      shapes: [
        {
          type: "line", x0: 1, x1: 1, y0: 0, y1: 1, yref: "paper",
          line: { color: P.border, width: 1, dash: "dot" },
        },
      ],
      hovermode: "closest",
    };

    Plotly.react("gep-1d", [trace], layout, CFG_2D);
  }

  // ── Plot: 3D surface n/d × p/n × error ──────────────────────────────────
  async function render3D() {
    const nl  = nlSel.property("value");
    const lam = getLambda();

    let rows;
    try {
      rows = (await loadCSV("surface_3d.csv"))
        .filter(r => r.nonlinearity === nl && Math.abs(r.lambda - lam) < 1e-12);
    } catch { rows = []; }

    if (!rows.length) return;

    // Build a grid for surface plot
    const p_ns = [...new Set(rows.map(r => r.p_n))].sort((a, b) => a - b);
    const n_ds = [...new Set(rows.map(r => r.n_d))].sort((a, b) => a - b);

    const zGrid = n_ds.map(nd =>
      p_ns.map(pn => {
        const hit = rows.find(r => Math.abs(r.p_n - pn) < 1e-9 && Math.abs(r.n_d - nd) < 1e-9);
        return hit ? hit.gen_error : null;
      })
    );

    const trace = {
      x: p_ns, y: n_ds, z: zGrid,
      type: "surface",
      colorscale: SURFACE_SCALE,
      reversescale: false,
      colorbar: {
        title: { text: "Error", font: { size: 11, family: P.font } },
        tickfont: { size: 10, family: P.font },
        thickness: 12,
        len: 0.75,
      },
      contours: {
        z: { show: true, usecolormap: true, highlightcolor: "#fff", project: { z: false } },
      },
      hovertemplate: "p/n = %{x:.2f}<br>n/d = %{y:.2f}<br>error = %{z:.4f}<extra></extra>",
    };

    const layout = {
      paper_bgcolor: P.bg,
      font: { family: P.font, color: P.text },
      margin: { l: 0, r: 0, t: 16, b: 0 },
      scene: {
        xaxis: { title: "p / n", titlefont: { size: 11 }, tickfont: { size: 9 }, gridcolor: P.grid },
        yaxis: { title: "n / d", titlefont: { size: 11 }, tickfont: { size: 9 }, gridcolor: P.grid },
        zaxis: { title: "Error",  titlefont: { size: 11 }, tickfont: { size: 9 }, gridcolor: P.grid },
        bgcolor: P.bg,
        camera: { eye: { x: 1.6, y: 1.6, z: 1.2 } },
        aspectmode: "manual",
        aspectratio: { x: 1.2, y: 1, z: 0.75 },
      },
    };

    Plotly.react("gep-3d", [trace], layout, CFG_3D);
  }

  // ── Wire events ──────────────────────────────────────────────────────────
  async function update() {
    await Promise.all([render1D(), render3D()]);
  }

  lambdaIn.on("input", () => { updateBadges(); update(); });
  deltaIn.on("input",  () => { updateBadges(); update(); });
  nlSel.on("change",   () => update());
  taskSel.on("change", () => { toggleDelta(); update(); });

  // ── Init ─────────────────────────────────────────────────────────────────
  updateBadges();
  toggleDelta();

  // Defer initial load so the DOM nodes exist
  setTimeout(() => update().catch(console.error), 200);

  return root.node();
}
