// Interactive Generalization Error Plots
// This module creates interactive visualizations for the generalization error analysis

function createInteractivePlots({ d3, Plotly }) {
  
  // Load and parse CSV data
  async function loadData(filename) {
    try {
      const response = await fetch(`data/${filename}`);
      if (!response.ok) throw new Error(`Failed to load ${filename}`);
      const text = await response.text();
      
      const lines = text.trim().split('\n');
      const headers = lines[0].split(',');
      const data = lines.slice(1).map(line => {
        const values = line.split(',');
        const row = {};
        headers.forEach((header, i) => {
          const value = values[i];
          // Try to parse as number, fallback to string
          row[header] = isNaN(value) ? value : parseFloat(value);
        });
        return row;
      });
      return data;
    } catch (error) {
      console.error(`Error loading data: ${error}`);
      return [];
    }
  }

  // Create the main container
  const mainDiv = d3.create('div')
    .style('width', '100%')
    .style('max-width', '1200px')
    .style('margin', 'auto');

  // Add title
  mainDiv.append('h3')
    .text('Interactive Generalization Error Analysis')
    .style('text-align', 'center')
    .style('margin-bottom', '20px');

  // Control panel
  const controlsDiv = mainDiv.append('div')
    .style('background', '#f8f9fa')
    .style('padding', '15px')
    .style('border-radius', '5px')
    .style('margin-bottom', '20px')
    .style('display', 'flex')
    .style('flex-wrap', 'wrap')
    .style('gap', '20px')
    .style('justify-content', 'center');

  // Controls
  const lambdaSlider = controlsDiv.append('div')
    .style('display', 'flex')
    .style('flex-direction', 'column')
    .style('align-items', 'center');
  
  lambdaSlider.append('label')
    .text('Regularization λ')
    .style('margin-bottom', '5px')
    .style('font-weight', 'bold');
  
  const lambdaInput = lambdaSlider.append('input')
    .attr('type', 'range')
    .attr('min', 0)
    .attr('max', 4)
    .attr('step', 1)
    .attr('value', 1)
    .style('width', '150px');
    
  const lambdaValue = lambdaSlider.append('div')
    .style('margin-top', '5px')
    .style('font-size', '12px');

  // Nonlinearity selector
  const nonlinDiv = controlsDiv.append('div')
    .style('display', 'flex')
    .style('flex-direction', 'column')
    .style('align-items', 'center');

  nonlinDiv.append('label')
    .text('Nonlinearity')
    .style('margin-bottom', '5px')
    .style('font-weight', 'bold');

  const nonlinSelect = nonlinDiv.append('select')
    .style('padding', '5px');

  nonlinSelect.selectAll('option')
    .data(['sign', 'erf'])
    .enter().append('option')
    .attr('value', d => d)
    .text(d => d);

  // Task selector
  const taskDiv = controlsDiv.append('div')
    .style('display', 'flex')
    .style('flex-direction', 'column')
    .style('align-items', 'center');

  taskDiv.append('label')
    .text('Task')
    .style('margin-bottom', '5px')
    .style('font-weight', 'bold');

  const taskSelect = taskDiv.append('select')
    .style('padding', '5px');

  taskSelect.selectAll('option')
    .data(['classification', 'regression'])
    .enter().append('option')
    .attr('value', d => d)
    .text(d => d);

  // Delta slider (for regression)
  const deltaSlider = controlsDiv.append('div')
    .style('display', 'flex')
    .style('flex-direction', 'column')
    .style('align-items', 'center');
  
  deltaSlider.append('label')
    .text('Noise Δ')
    .style('margin-bottom', '5px')
    .style('font-weight', 'bold');
  
  const deltaInput = deltaSlider.append('input')
    .attr('type', 'range')
    .attr('min', 0)
    .attr('max', 3)
    .attr('step', 1)
    .attr('value', 0)
    .style('width', '150px');
    
  const deltaValue = deltaSlider.append('div')
    .style('margin-top', '5px')
    .style('font-size', '12px');

  // Plot containers
  const plotsDiv = mainDiv.append('div')
    .style('display', 'flex')
    .style('flex-direction', 'column')
    .style('gap', '30px');

  // 1D plot container
  const plot1dDiv = plotsDiv.append('div')
    .attr('id', 'plot1d')
    .style('height', '400px')
    .style('border', '1px solid #ddd')
    .style('border-radius', '5px');

  // 3D plot container
  const plot3dDiv = plotsDiv.append('div')
    .attr('id', 'plot3d')
    .style('height', '500px')
    .style('border', '1px solid #ddd')
    .style('border-radius', '5px');

  plot3dDiv.append('h4')
    .text('3D Surface: n/d vs p/n vs Generalization Error')
    .style('text-align', 'center')
    .style('margin', '10px 0');

  // Lambda values mapping
  const lambdaValues = [1e-4, 1e-3, 1e-2, 1e-1, 1.0];
  const deltaValues = [0.0, 0.1, 0.5, 1.0];

  // Update lambda value display
  function updateLambdaValue() {
    const idx = parseInt(lambdaInput.property('value'));
    lambdaValue.text(`λ = ${lambdaValues[idx]}`);
  }

  // Update delta value display  
  function updateDeltaValue() {
    const idx = parseInt(deltaInput.property('value'));
    deltaValue.text(`Δ = ${deltaValues[idx]}`);
  }

  // Update visibility of delta slider
  function updateDeltaVisibility() {
    const task = taskSelect.property('value');
    deltaSlider.style('display', task === 'regression' ? 'flex' : 'none');
  }

  // Plot 1D data
  async function plot1D() {
    const task = taskSelect.property('value');
    const nonlinearity = nonlinSelect.property('value');
    const lambdaIdx = parseInt(lambdaInput.property('value'));
    const deltaIdx = parseInt(deltaInput.property('value'));
    
    const lambda = lambdaValues[lambdaIdx];
    const delta = deltaValues[deltaIdx];

    try {
      let data;
      if (task === 'classification') {
        data = await loadData('classification_1d.csv');
        // Filter data
        data = data.filter(d => 
          d.nonlinearity === nonlinearity && 
          Math.abs(d.lambda - lambda) < 1e-10
        );
      } else {
        data = await loadData('regression_1d.csv');
        // Filter data
        data = data.filter(d => 
          d.nonlinearity === nonlinearity && 
          Math.abs(d.lambda - lambda) < 1e-10 && 
          Math.abs(d.Delta - delta) < 1e-10
        );
      }

      if (data.length === 0) {
        console.log('No data found for current selection');
        return;
      }

      // Sort by p_n
      data.sort((a, b) => a.p_n - b.p_n);

      const trace = {
        x: data.map(d => d.p_n),
        y: data.map(d => d.gen_error),
        type: 'scatter',
        mode: 'lines+markers',
        name: task === 'classification' ? 
          `${nonlinearity} (λ=${lambda})` : 
          `${nonlinearity} (λ=${lambda}, Δ=${delta})`,
        line: { width: 3 },
        marker: { size: 6 }
      };

      const layout = {
        title: {
          text: `${task.charAt(0).toUpperCase() + task.slice(1)} Generalization Error`,
          font: { size: 16 }
        },
        xaxis: { 
          title: 'p/n (features per sample)',
          titlefont: { size: 14 }
        },
        yaxis: { 
          title: 'Generalization Error',
          titlefont: { size: 14 }
        },
        margin: { l: 60, r: 50, t: 50, b: 50 },
        hovermode: 'closest',
        annotations: []
      };

      // Add annotation for double descent spike if visible
      const spikeData = data.filter(d => d.p_n >= 0.8 && d.p_n <= 1.2);
      if (spikeData.length > 0 && lambda <= 1e-2) {
        const maxPoint = spikeData.reduce((max, d) => d.gen_error > max.gen_error ? d : max);
        layout.annotations.push({
          x: maxPoint.p_n,
          y: maxPoint.gen_error,
          text: "Double Descent<br>Spike",
          showarrow: true,
          arrowhead: 2,
          arrowsize: 1,
          arrowwidth: 2,
          arrowcolor: "#ff0000",
          ax: 20,
          ay: -30,
          font: { color: "#ff0000", size: 12 }
        });
      }

      const config = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['pan2d', 'select2d', 'lasso2d', 'resetScale2d']
      };

      await Plotly.newPlot('plot1d', [trace], layout, config);

    } catch (error) {
      console.error('Error in plot1D:', error);
    }
  }

  // Plot 3D surface
  async function plot3D() {
    const nonlinearity = nonlinSelect.property('value');
    const lambdaIdx = parseInt(lambdaInput.property('value'));
    const lambda = lambdaValues[lambdaIdx];

    try {
      const data = await loadData('surface_3d.csv');
      
      // Filter data
      const filteredData = data.filter(d => 
        d.nonlinearity === nonlinearity && 
        Math.abs(d.lambda - lambda) < 1e-10
      );

      if (filteredData.length === 0) {
        console.log('No 3D data found for current selection');
        return;
      }

      const trace = {
        x: filteredData.map(d => d.p_n),
        y: filteredData.map(d => d.n_d),
        z: filteredData.map(d => d.gen_error),
        type: 'scatter3d',
        mode: 'markers',
        marker: {
          size: 4,
          color: filteredData.map(d => d.gen_error),
          colorscale: 'Viridis',
          showscale: true,
          colorbar: { title: 'Generalization<br>Error' }
        },
        text: filteredData.map(d => 
          `p/n: ${d.p_n.toFixed(2)}<br>n/d: ${d.n_d.toFixed(2)}<br>Error: ${d.gen_error.toFixed(4)}`
        ),
        hovertemplate: '%{text}<extra></extra>'
      };

      const layout = {
        scene: {
          xaxis: { title: 'p/n (features per sample)' },
          yaxis: { title: 'n/d (samples per latent dimension)' },
          zaxis: { title: 'Generalization Error' },
          camera: {
            eye: { x: 1.5, y: 1.5, z: 1.5 }
          }
        },
        title: `3D Surface: ${nonlinearity} nonlinearity (λ=${lambda})`,
        margin: { l: 0, r: 0, t: 50, b: 0 }
      };

      const config = {
        responsive: true,
        displayModeBar: true
      };

      await Plotly.newPlot('plot3d', [trace], layout, config);

    } catch (error) {
      console.error('Error in plot3D:', error);
    }
  }

  // Update all plots
  async function updatePlots() {
    await plot1D();
    await plot3D();
  }

  // Event listeners with error handling
  lambdaInput.on('input', function() {
    updateLambdaValue();
    updatePlots().catch(console.error);
  });

  deltaInput.on('input', function() {
    updateDeltaValue();
    updatePlots().catch(console.error);
  });

  nonlinSelect.on('change', () => updatePlots().catch(console.error));
  taskSelect.on('change', function() {
    updateDeltaVisibility();
    updatePlots().catch(console.error);
  });

  // Initialize
  updateLambdaValue();
  updateDeltaValue();
  updateDeltaVisibility();
  
  // Add a loading message instead of immediate plot loading
  const loadingDiv = mainDiv.append('div')
    .style('text-align', 'center')
    .style('padding', '40px')
    .style('color', '#666')
    .text('Loading interactive plots...');
  
  // Load plots after a delay to prevent blocking
  setTimeout(() => {
    loadingDiv.remove();
    updatePlots().catch(error => {
      console.error('Error loading plots:', error);
      mainDiv.append('div')
        .style('color', 'red')
        .style('padding', '20px')
        .style('text-align', 'center')
        .text(`Error loading interactive plots: ${error.message}`);
    });
  }, 500);

  // Return the main container node like the working generativeModelViz
  return mainDiv.node();
}

// Export for OJS compatibility
export { createInteractivePlots };