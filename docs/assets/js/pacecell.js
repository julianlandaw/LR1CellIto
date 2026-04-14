(function () {
  'use strict';

  const DEFAULTS = {
    inafac: 1.0,
    itofac: 0.0,
    icalfac: 1.0,
    ikfac: 1.0,
    ikifac: 1.0,
    tauXfac: 1.0,
    yshift: 0.0,
    beats: 10,
    pcl: 500,
  };

  const PRESETS = {
    default: { ...DEFAULTS },
    chaos: {
      inafac: 1.0,
      itofac: 1.05,
      pcl: 378,
      beats: 50,
      tauXfac: 5.0,
      icalfac: 1.15,
      ikifac: 2.2,
      ikfac: 1.0,
      yshift: 8.0,
    },
    eads: {
      inafac: 1.0,
      itofac: 0.0,
      icalfac: 1.0,
      ikfac: 1.0,
      ikifac: 1.0,
      tauXfac: 10.0,
      yshift: 0.0,
      beats: 15,
      pcl: 1575,
    },
  };

  const INPUT_IDS = ['inafac', 'itofac', 'icalfac', 'ikfac', 'ikifac', 'tauXfac', 'yshift', 'beats', 'pcl'];
  const SAMPLE_INTERVAL_MS = 5;
  const WARMUP_MS = 10000;
  const AUTO_UPDATE_DEBOUNCE_MS = 450;

  let ui = {};
  let currentPreset = 'default';
  let isRunning = false;
  let pendingRun = false;
  let debouncedTimer = null;
  let lastResult = null;

  function init() {
    ui = {
      form: document.getElementById('controlsForm'),
      plot: document.getElementById('myDiv1'),
      status: document.getElementById('statusIndicator'),
      statusText: document.getElementById('statusText'),
      autoUpdate: document.getElementById('autoUpdate'),
      runBtn: document.getElementById('runBtn'),
      resetBtn: document.getElementById('resetBtn'),
      exportBtn: document.getElementById('exportBtn'),
      presetButtons: {
        default: document.getElementById('presetDefaultBtn'),
        chaos: document.getElementById('presetChaosBtn'),
        eads: document.getElementById('presetEadsBtn'),
      },
      metrics: {
        apdMin: document.getElementById('metricApdMin'),
        apdMax: document.getElementById('metricApdMax'),
        apdRange: document.getElementById('metricApdRange'),
        beatCount: document.getElementById('metricBeatCount'),
        preset: document.getElementById('metricPreset'),
        presetNote: document.getElementById('metricPresetNote'),
      },
      inputs: {},
    };

    INPUT_IDS.forEach((id) => {
      ui.inputs[id] = document.getElementById(id);
    });

    setFormValues(DEFAULTS);
    setActivePreset('default');
    bindEvents();
    setStatus('ready', 'Ready to simulate');
    startSimulation();
  }

  function bindEvents() {
    ui.form.addEventListener('submit', (event) => {
      event.preventDefault();
      queueRun(true);
    });

    ui.resetBtn.addEventListener('click', () => {
      setFormValues(DEFAULTS);
      setActivePreset('default');
      queueRun(true);
    });

    ui.exportBtn.addEventListener('click', exportPlotImage);

    Object.entries(ui.presetButtons).forEach(([presetName, button]) => {
      button.addEventListener('click', () => {
        applyPreset(presetName);
      });
    });

    INPUT_IDS.forEach((id) => {
      const input = ui.inputs[id];
      ['input', 'change'].forEach((eventName) => {
        input.addEventListener(eventName, handleInputChange);
      });
    });

    window.addEventListener('resize', debounce(() => {
      if (ui.plot && window.Plotly && ui.plot.data) {
        Plotly.Plots.resize(ui.plot);
      }
    }, 150));
  }

  function handleInputChange() {
    detectPresetFromCurrentValues();

    if (ui.autoUpdate.checked) {
      queueRun();
    } else {
      setStatus('ready', 'Parameter updated. Click “Run simulation” to refresh plots.');
    }
  }

  function queueRun(immediate = false) {
    if (debouncedTimer) {
      clearTimeout(debouncedTimer);
      debouncedTimer = null;
    }

    if (immediate) {
      startSimulation();
      return;
    }

    debouncedTimer = setTimeout(() => {
      startSimulation();
    }, AUTO_UPDATE_DEBOUNCE_MS);
  }

  function startSimulation() {
    if (isRunning) {
      pendingRun = true;
      return;
    }

    isRunning = true;
    updateRunningState(true);
    setStatus('running', 'Running simulation…');

    window.requestAnimationFrame(() => {
      try {
        const inputValues = getInputValues();
        validateInputs(inputValues);

        const result = runSimulation(inputValues);
        lastResult = result;
        renderMetrics(result.metrics);
        renderPlots(result);

        setStatus('ready', `Simulation complete in ${result.runtimeMs.toFixed(0)} ms.`);
      } catch (error) {
        console.error(error);
        setStatus('error', error && error.message ? error.message : 'Simulation failed.');
      } finally {
        isRunning = false;
        updateRunningState(false);

        if (pendingRun) {
          pendingRun = false;
          startSimulation();
        }
      }
    });
  }

  function runSimulation(inputValues) {
    if (typeof LR1CellIto !== 'function') {
      throw new Error('LR1CellIto.js is not loaded. Check the script path near the bottom of index.html.');
    }

    if (typeof Plotly === 'undefined') {
      throw new Error('Plotly is not loaded. Update the local Plotly script path in index.html.');
    }

    const startedAt = performance.now();
    const cell = new LR1CellIto();

    cell.inafac = inputValues.inafac;
    cell.itofac = inputValues.itofac;
    cell.icalfac = inputValues.icalfac;
    cell.ikfac = inputValues.ikfac;
    cell.ikifac = inputValues.ikifac;
    cell.tauXfac = inputValues.tauXfac;
    cell.yshift = inputValues.yshift;

    const warmupSteps = Math.ceil(WARMUP_MS / cell.dt);
    for (let i = 0; i < warmupSteps; i += 1) {
      cell.stepdt(0);
    }

    let t = 0;
    let nextStim = 100;

    for (t = cell.dt; t < inputValues.pcl * 10; t += cell.dt) {
      if (t > nextStim - cell.dt / 2 && t < nextStim + cell.stimduration - cell.dt / 2) {
        cell.stepdt(cell.stimulus);
      } else {
        cell.stepdt(0);
      }

      if (t > nextStim + cell.stimduration + cell.dt / 2) {
        nextStim += inputValues.pcl;
      }
    }

    const ts = [0];
    const vs = [cell.v];
    const xs = [cell.xr];
    const apds = [];
    const xsInit = [];
    const apdCounter = [];

    let oldV = cell.v;
    let newV = cell.v;
    let startApdTime = null;
    let minApd = Infinity;
    let maxApd = -Infinity;
    let minX = Infinity;
    let maxX = -Infinity;
    let lastSampleTime = SAMPLE_INTERVAL_MS;
    nextStim = 100;

    for (t = cell.dt; t < inputValues.pcl * inputValues.beats; t += cell.dt) {
      oldV = cell.v;

      if (t > nextStim - cell.dt / 2 && t < nextStim + cell.stimduration - cell.dt / 2) {
        cell.stepdt(cell.stimulus);
      } else {
        cell.stepdt(0);
      }

      newV = cell.v;

      if (newV > cell.threshold && oldV < cell.threshold) {
        startApdTime = t;
        xsInit.push(cell.xr);
        minX = Math.min(minX, cell.xr);
        maxX = Math.max(maxX, cell.xr);
      }

      if (newV < cell.threshold && oldV > cell.threshold && startApdTime !== null) {
        const thisApd = t - startApdTime;
        apds.push(thisApd);
        apdCounter.push(apds.length);
        minApd = Math.min(minApd, thisApd);
        maxApd = Math.max(maxApd, thisApd);
        startApdTime = null;
      }

      if (t > nextStim + cell.stimduration + cell.dt / 2) {
        nextStim += inputValues.pcl;
      }

      if (t >= lastSampleTime - cell.dt / 2) {
        ts.push(t / 1000);
        vs.push(cell.v);
        xs.push(cell.xr);
        lastSampleTime += SAMPLE_INTERVAL_MS;
      }
    }

    if (!Number.isFinite(minApd)) {
      minApd = NaN;
      maxApd = NaN;
    }

    if (!Number.isFinite(minX)) {
      minX = NaN;
      maxX = NaN;
    }

    const apdRange = Number.isFinite(minApd) && Number.isFinite(maxApd) ? maxApd - minApd : NaN;

    return {
      ts,
      vs,
      xs,
      apds,
      xsInit,
      apdCounter,
      runtimeMs: performance.now() - startedAt,
      inputs: inputValues,
      metrics: {
        minApd,
        maxApd,
        apdRange,
        beatCount: apds.length,
        minX,
        maxX,
        preset: currentPreset,
      },
    };
  }

  function renderPlots(result) {
    const yRange = computeApdRange(result.apds);
    const xScatterRange = computeXScatterRange(result.metrics.minX, result.metrics.maxX);

    const traces = [
      {
        x: result.ts,
        y: result.vs,
        type: 'scatter',
        mode: 'lines',
        name: 'Action potential',
        line: { color: '#2563eb', width: 2.25 },
        hovertemplate: 'Time: %{x:.3f} s<br>Voltage: %{y:.2f} mV<extra></extra>',
      },
      {
        x: result.ts,
        y: result.xs,
        type: 'scatter',
        mode: 'lines',
        name: 'X gate',
        xaxis: 'x2',
        yaxis: 'y2',
        line: { color: '#7c3aed', width: 2.0 },
        hovertemplate: 'Time: %{x:.3f} s<br>X: %{y:.4f}<extra></extra>',
      },
      {
        x: result.apdCounter,
        y: result.apds,
        type: 'scatter',
        mode: 'markers+lines',
        name: 'APD by beat',
        xaxis: 'x3',
        yaxis: 'y3',
        marker: { color: '#dc2626', size: 7 },
        line: { color: '#fca5a5', width: 1.5 },
        hovertemplate: 'Beat: %{x}<br>APD: %{y:.2f} ms<extra></extra>',
      },
      {
        x: result.xsInit,
        y: result.apds,
        type: 'scatter',
        mode: 'markers',
        name: 'APD vs initial X',
        xaxis: 'x4',
        yaxis: 'y4',
        marker: { color: '#059669', size: 7, opacity: 0.85 },
        hovertemplate: 'Initial X: %{x:.4f}<br>APD: %{y:.2f} ms<extra></extra>',
      },
    ];

    const layout = {
      title: {
        text: 'LR1 simulation output',
        font: { size: 22 },
        x: 0.02,
        xanchor: 'left',
      },
      paper_bgcolor: '#ffffff',
      plot_bgcolor: '#fbfdff',
      margin: { l: 64, r: 24, t: 58, b: 48 },
      showlegend: false,
      height: 820,
      grid: {
        rows: 4,
        columns: 1,
        pattern: 'independent',
        roworder: 'top to bottom',
        ygap: 0.08,
      },
      xaxis: {
        title: 'Time (s)',
        showgrid: true,
        gridcolor: '#e2e8f0',
        zeroline: false,
      },
      yaxis: {
        title: 'Voltage (mV)',
        showgrid: true,
        gridcolor: '#e2e8f0',
        zeroline: false,
      },
      xaxis2: {
        title: 'Time (s)',
        showgrid: true,
        gridcolor: '#e2e8f0',
        zeroline: false,
      },
      yaxis2: {
        title: 'X',
        showgrid: true,
        gridcolor: '#e2e8f0',
        zeroline: false,
      },
      xaxis3: {
        title: 'Beat number',
        showgrid: true,
        gridcolor: '#e2e8f0',
        zeroline: false,
      },
      yaxis3: {
        title: 'APD (ms)',
        showgrid: true,
        gridcolor: '#e2e8f0',
        zeroline: false,
        range: yRange,
      },
      xaxis4: {
        title: 'Initial X',
        showgrid: true,
        gridcolor: '#e2e8f0',
        zeroline: false,
        range: xScatterRange,
      },
      yaxis4: {
        title: 'APD (ms)',
        showgrid: true,
        gridcolor: '#e2e8f0',
        zeroline: false,
        range: yRange,
      },
      annotations: [
        {
          text: `Preset: ${formatPresetName(currentPreset)} · PCL: ${formatNumber(result.inputs.pcl, 0)} ms · Beats: ${formatNumber(result.inputs.beats, 0)}`,
          xref: 'paper',
          yref: 'paper',
          x: 1,
          y: 1.08,
          xanchor: 'right',
          showarrow: false,
          font: { color: '#475569', size: 12 },
        },
      ],
    };

    const config = {
      responsive: true,
      displaylogo: false,
      toImageButtonOptions: {
        format: 'png',
        filename: 'lr1-simulation',
        scale: 2,
      },
    };

    if (ui.plot && ui.plot.data) {
      Plotly.react(ui.plot, traces, layout, config);
    } else {
      Plotly.newPlot(ui.plot, traces, layout, config);
    }
  }

  function renderMetrics(metrics) {
    ui.metrics.apdMin.textContent = formatMetric(metrics.minApd, ' ms');
    ui.metrics.apdMax.textContent = formatMetric(metrics.maxApd, ' ms');
    ui.metrics.apdRange.textContent = formatMetric(metrics.apdRange, ' ms');
    ui.metrics.beatCount.textContent = Number.isFinite(metrics.beatCount) ? String(metrics.beatCount) : '—';
    ui.metrics.preset.textContent = formatPresetName(metrics.preset);

    if (Number.isFinite(metrics.minX) && Number.isFinite(metrics.maxX)) {
      ui.metrics.presetNote.textContent = `Initial X range: ${metrics.minX.toFixed(4)} to ${metrics.maxX.toFixed(4)}`;
    } else {
      ui.metrics.presetNote.textContent = 'No complete APDs detected in this run';
    }
  }

  function setFormValues(values) {
    INPUT_IDS.forEach((id) => {
      const value = values[id];
      if (!ui.inputs[id]) return;

      ui.inputs[id].value = Number.isFinite(value)
        ? (Number.isInteger(value) ? String(value) : value.toFixed(2))
        : '';
    });
  }

  function getInputValues() {
    const values = {};
    INPUT_IDS.forEach((id) => {
      values[id] = Number.parseFloat(ui.inputs[id].value);
    });

    values.beats = Math.round(values.beats);
    values.pcl = Math.round(values.pcl);
    return values;
  }

  function validateInputs(values) {
    const ranges = {
      inafac: [0, 5],
      itofac: [0, 5],
      icalfac: [0, 5],
      ikfac: [0, 5],
      ikifac: [0, 5],
      tauXfac: [0.1, 25],
      yshift: [-40, 40],
      beats: [1, 250],
      pcl: [50, 5000],
    };

    for (const [key, [min, max]] of Object.entries(ranges)) {
      const value = values[key];
      if (!Number.isFinite(value)) {
        throw new Error(`Please enter a valid numeric value for ${key}.`);
      }
      if (value < min || value > max) {
        throw new Error(`${key} must be between ${min} and ${max}.`);
      }
    }
  }

  function applyPreset(presetName) {
    const preset = PRESETS[presetName];
    if (!preset) return;

    setFormValues(preset);
    setActivePreset(presetName);
    queueRun(true);
  }

  function setActivePreset(presetName) {
    currentPreset = presetName;

    Object.entries(ui.presetButtons).forEach(([name, button]) => {
      button.classList.remove('btn-pill-active', 'btn-secondary', 'btn-ghost');

      if (name === presetName) {
        button.classList.add('btn-pill-active');
      } else if (name === 'default') {
        button.classList.add('btn-secondary');
      } else {
        button.classList.add('btn-ghost');
      }
    });
  }

  function detectPresetFromCurrentValues() {
    const currentValues = getInputValues();

    const matchedPreset = Object.entries(PRESETS).find(([, preset]) => {
      return isSameConfig(currentValues, preset);
    });

    if (matchedPreset) {
      setActivePreset(matchedPreset[0]);
      return;
    }

    currentPreset = 'custom';
    Object.entries(ui.presetButtons).forEach(([name, button]) => {
      button.classList.remove('btn-pill-active', 'btn-secondary', 'btn-ghost');
      button.classList.add(name === 'default' ? 'btn-secondary' : 'btn-ghost');
    });

    ui.metrics.preset.textContent = 'Custom';
    ui.metrics.presetNote.textContent = 'Manual parameter edits';
  }

  function exportPlotImage() {
    if (!ui.plot || !ui.plot.data || typeof Plotly === 'undefined') {
      setStatus('error', 'Run a simulation before exporting the plot.');
      return;
    }

    Plotly.downloadImage(ui.plot, {
      format: 'png',
      filename: 'lr1-simulation',
      width: 1600,
      height: 900,
      scale: 2,
    });
  }

  function updateRunningState(running) {
    ui.runBtn.disabled = running;
    ui.resetBtn.disabled = running;
    ui.exportBtn.disabled = running && !lastResult;

    Object.values(ui.presetButtons).forEach((button) => {
      button.disabled = running;
    });
  }

  function setStatus(kind, message) {
    ui.status.classList.remove('ready', 'running', 'error');
    ui.status.classList.add(kind);
    ui.statusText.textContent = message;
  }

  function computeApdRange(apds) {
    if (!Array.isArray(apds) || apds.length === 0) return undefined;

    const min = Math.min(...apds);
    const max = Math.max(...apds);
    const pad = Math.max(25, (max - min) * 0.2);
    return [Math.max(0, min - pad), max + pad];
  }

  function computeXScatterRange(minX, maxX) {
    if (!Number.isFinite(minX) || !Number.isFinite(maxX)) return undefined;

    const pad = Math.max(0.01, (maxX - minX) * 0.2);
    return [minX - pad, maxX + pad];
  }

  function formatMetric(value, suffix = '') {
    return Number.isFinite(value) ? `${value.toFixed(1)}${suffix}` : '—';
  }

  function formatNumber(value, decimals) {
    return Number.isFinite(value) ? value.toFixed(decimals) : '—';
  }

  function formatPresetName(preset) {
    switch (preset) {
      case 'default':
        return 'Default';
      case 'chaos':
        return 'Chaos';
      case 'eads':
        return 'EADs';
      case 'custom':
        return 'Custom';
      default:
        return 'Default';
    }
  }

  function isSameConfig(a, b) {
    return INPUT_IDS.every((key) => {
      return Math.abs(Number(a[key]) - Number(b[key])) < 1e-9;
    });
  }

  function debounce(fn, delay) {
    let timerId = null;
    return function debounced(...args) {
      if (timerId) clearTimeout(timerId);
      timerId = setTimeout(() => fn.apply(this, args), delay);
    };
  }

  window.addEventListener('DOMContentLoaded', init);
})();