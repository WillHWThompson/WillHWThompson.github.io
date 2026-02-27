// Helper function for random Gaussian
function randomGaussian() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

// Generate random data
function generateData(nPoints = 50, thetaAngle = Math.PI/4) {
  // Generate latent points c^μ ~ N(0, I_2)
  const cPoints = [];
  for (let i = 0; i < nPoints; i++) {
    cPoints.push([randomGaussian(), randomGaussian()]);
  }
  
  // Ground truth vector θ on unit circle
  const theta = [Math.cos(thetaAngle), Math.sin(thetaAngle)];
  
  // Random projection matrix F (3x2) - projects from 2D to 3D
  // Increase spread by removing division and multiplying by 2
  const F = [];
  for (let i = 0; i < 3; i++) {
    F[i] = [];
    for (let j = 0; j < 2; j++) {
      F[i][j] = randomGaussian() * 2;
    }
  }
  return { cPoints, theta, F };
}

// Main visualization function
// Helper to apply nonlinearity
function applyNonlinearity(x, type) {
  if (type === "identity") return x;
  if (type === "tanh") return Math.tanh(x);
  if (type === "relu") return Math.max(0, x);
  if (type === "sign") return Math.sign(x);
  return x;
}

// Compute projections and labels
function computeProjections(data, f0Type, sigmaType) {
  const { cPoints, theta, F } = data;
  return cPoints.map(c => {
    // Project to 3D: x = sigma(F^T c)
    const xRaw = [
      F[0][0]*c[0] + F[0][1]*c[1],
      F[1][0]*c[0] + F[1][1]*c[1],
      F[2][0]*c[0] + F[2][1]*c[1]
    ];
    const x = xRaw.map(v => applyNonlinearity(v, sigmaType));
    // Label: y = f0((1/d) c · theta)
    const d = c.length;
    const dot = c[0]*theta[0] + c[1]*theta[1];
    const y = applyNonlinearity(dot/d, f0Type);
    return { c, x, y };
  });
}
// Only one valid generativeModelViz definition, no stray references to height or export
function generativeModelViz({
  sigmaType = "identity",
  f0Type = "identity", 
  thetaAngle = Math.PI/4,
  rerollTrigger = 0,
  d3,
  Plotly
}) {
  
  // Generate or regenerate data when button is clicked
  const data = generateData(50, thetaAngle);
  const projectedData = computeProjections(data, f0Type, sigmaType);
  const { theta } = data;
  
  // Create color scale based on labels
  const yValues = projectedData.map(d => d.y);
  const yExtent = d3.extent(yValues);
  const colorScale = d3.scaleSequential(d3.interpolateViridis)
    .domain(yExtent);
  
  // Calculate the z-range from projected points to place latent points at bottom
  const zValues = projectedData.map(d => d.x[2]);
  const zExtent = d3.extent(zValues);
  const latentZ = zExtent[0] - 1; // Place latent points 1 unit below minimum projected z
  
  // Set latent points at the lowest z value shown in the plot
  const minZ = -4; // Match Plotly z axis minimum

  // Prepare data for Plotly
  const traces = [];

  traces.push({
    x: projectedData.map(d => d.c[0]),
    y: projectedData.map(d => d.c[1]),
    z: projectedData.map(() => minZ),
    mode: 'markers',
    type: 'scatter3d',
    marker: {
      size: 5,
      color: projectedData.map(d => d.y),
      colorscale: 'Cividis', // Sleek colormap
      colorbar: {
        title: 'Label',
        x: 1.18
      },
      line: {
        color: '#222',
        width: 0.5
      }
    },
    name: 'Latent',
    text: projectedData.map(d => `Label: ${d.y.toFixed(2)}`),
    hovertemplate: '%{text}<extra></extra>'
  });
  
  // 2. Ground truth vector θ EXACTLY on x-y plane (z=0)
  // 2. Ground truth vector θ on latent plane (z=-4)
  traces.push({
    x: [0, theta[0] * 2],
    y: [0, theta[1] * 2],
    z: [minZ, minZ],
    mode: 'lines+markers',
    type: 'scatter3d',
    line: {
      color: '#0077b6',
      width: 2
    },
    marker: {
      size: [3, 8],
      color: '#0077b6',
      symbol: ['circle', 'diamond']
    },
    name: 'θ',
    hoverinfo: 'skip'
  });
  
  // 3. Projected points (x^μ) in 3D space
  traces.push({
    x: projectedData.map(d => d.x[0]),
    y: projectedData.map(d => d.x[1]),
    z: projectedData.map(d => d.x[2]),
    mode: 'markers',
    type: 'scatter3d',
    marker: {
      size: 5,
      color: projectedData.map(d => d.y),
      colorscale: 'Turbo', // Sleek colormap
      opacity: 0.7,
      line: {
        color: '#222',
        width: 0.5
      }
    },
    name: 'Projected',
    text: projectedData.map(d => `Label: ${d.y.toFixed(2)}`),
    hovertemplate: '%{text}<extra></extra>'
  });
  
  // 4. Reference plane at z=0 (latent space) to make it visually clear
  const planeSize = 4;
  const planeResolution = 2;
  const planeX = [];
  const planeY = [];
  const planeZ = [];
  
  for (let i = 0; i <= planeResolution; i++) {
    for (let j = 0; j <= planeResolution; j++) {
      planeX.push(-planeSize + (i * 2 * planeSize) / planeResolution);
      planeY.push(-planeSize + (j * 2 * planeSize) / planeResolution);  
      planeZ.push(0); // Always at z=0
    }
  }
  
  traces.push({
    x: planeX,
    y: planeY,
    z: planeZ,
    mode: 'markers',
    type: 'scatter3d',
    marker: {
      size: 1.5,
      color: '#eee',
      opacity: 0.15
    },
    showlegend: false,
    hoverinfo: 'skip',
    name: 'Latent plane'
  });
  
  // 5. Clear dotted lines connecting each latent point (z=0) to its projected point (3D)
  projectedData.forEach((d, i) => {
    traces.push({
      x: [d.x[0], d.c[0]],
      y: [d.x[1], d.c[1]],
      z: [d.x[2], minZ],
      mode: 'lines',
      type: 'scatter3d',
      line: {
        color: '#888',
        width: 1,
        dash: 'dot'
      },
      showlegend: false,
      hoverinfo: 'skip',
      opacity: 0.3
    });
  });
  
  // Layout configuration
  const layout = {
    title: {
      text: 'Generative Model Visualization',
      font: { size: 15, color: '#222' }
    },
    scene: {
      xaxis: { 
        title: 'Latent Dim 1',
        range: [-4, 4],
        showgrid: false,
        zeroline: false,
        color: '#222',
        backgroundcolor: '#fff'
      },
      yaxis: { 
        title: 'Latent Dim 2',
        range: [-4, 4],
        showgrid: false,
        zeroline: false,
        color: '#222',
        backgroundcolor: '#fff'
      },
      zaxis: { 
        title: 'Projected Dim',
        range: [-4, 4],
        showgrid: false,
        zeroline: false,
        color: '#222',
        backgroundcolor: '#fff'
      },
      camera: {
        eye: { x: 1.8, y: 1.8, z: 1.5 }
      },
      aspectmode: 'cube',
      bgcolor: '#fff'
    },
    width: 800,
    height: 600,
    margin: { l: 0, r: 0, t: 40, b: 0 },
    legend: {
      x: 0.85,
      y: 0.95,
      bgcolor: '#fff',
      bordercolor: '#eee',
      borderwidth: 1,
      font: { size: 12, color: '#222' }
    },
    hovermode: false // Remove abscissa/crosshair hover lines
  };
  
  // Create container
  const container = d3.create("div")
    .style("width", "800px")
    .style("height", "600px")
    .style("background", "#fff");
  
  // Plot with Plotly
  Plotly.newPlot(container.node(), traces, layout, {
    responsive: true,
    displayModeBar: true
  });
  
  // Add info panel
  // Minimal info panel
  const infoPanel = d3.create("div")
    .style("margin-top", "8px")
    .style("padding", "6px")
    .style("background", "#fff")
    .style("border-radius", "3px")
    .style("font-size", "13px")
    .style("color", "#222");

  infoPanel.append("div")
    .style("font-weight", "bold")
    .text("Parameters:");

  infoPanel.append("div")
    .text(`f₀: ${f0Type}, σ: ${sigmaType}, θ: ${(thetaAngle * 180 / Math.PI).toFixed(1)}°`);

  // Combine plot and info
  const fullContainer = d3.create("div");
  fullContainer.node().appendChild(container.node());
  fullContainer.node().appendChild(infoPanel.node());

  return fullContainer.node();
}
// End of generativeModelViz function