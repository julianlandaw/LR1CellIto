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
  const FIELD_LABELS = {
    inafac: 'I Na multiplier', itofac: 'I to multiplier', icalfac: 'I Ca,L multiplier',
    ikfac: 'I Kr multiplier', ikifac: 'I K1 multiplier', tauXfac: 'X gate scaling',
    yshift: 'y shift', beats: 'Beats', pcl: 'PCL',
  };
  const AUTO_UPDATE_DEBOUNCE_MS = 450;

  let ui = {};
  let currentPreset = 'default';
  let debouncedTimer = null;
  let activeWorker = null;
  let runId = 0;
  let lastResult = null;
  let resultsStale = true;
  let plotQueue = Promise.resolve();

  function init() {
    ui = {
      form: document.getElementById('controlsForm'),
      plot: document.getElementById('myDiv1'),
      status: document.getElementById('statusIndicator'),
      statusText: document.getElementById('statusText'),
      inputError: document.getElementById('inputError'),
      resultState: document.getElementById('resultState'),
      apdTableBody: document.getElementById('apdTableBody'),
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
        beatNote: document.getElementById('metricBeatNote'),
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
      setResultsStale(true);
      queueRun(true);
    });

    ui.exportBtn.addEventListener('click', exportPlotImage);
    ui.autoUpdate.addEventListener('change', () => {
      if (ui.autoUpdate.checked && resultsStale) queueRun(true);
    });

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
        Plotly.relayout(ui.plot, { height: plotHeight() }).catch(console.error);
      }
    }, 150));
  }

  function handleInputChange() {
    clearInputError();
    detectPresetFromCurrentValues();
    setResultsStale(true);

    if (ui.autoUpdate.checked) {
      queueRun();
    } else {
      setStatus('ready', 'Parameters changed. Run the simulation to refresh results.');
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
      debouncedTimer = null;
      startSimulation();
    }, AUTO_UPDATE_DEBOUNCE_MS);
  }

  function startSimulation() {
    let inputValues;
    try {
      inputValues = getInputValues(true);
      clearInputError();
    } catch (error) {
      cancelActiveWorker();
      setResultsStale(true);
      showInputError(error);
      updateRunningState(false);
      return;
    }

    cancelActiveWorker();
    const id = ++runId;
    const preset = currentPreset;
    let worker;
    try {
      worker = new Worker('assets/js/simulation-worker.js');
    } catch (error) {
      setStatus('error', `Could not start the simulation worker: ${error.message}`);
      return;
    }
    activeWorker = worker;
    updateRunningState(true);
    setStatus('running', 'Running simulation…');
    let shownProgress = 0;

    worker.onmessage = async ({ data }) => {
      if (!data || data.id !== runId || worker !== activeWorker) return;
      if (data.kind === 'progress') {
        const progress = Math.floor(data.percent / 10) * 10;
        if (progress > shownProgress) {
          shownProgress = progress;
          setStatus('running', `Running simulation… ${progress}%`);
        }
        return;
      }
      if (data.kind === 'error') {
        finishWithError(data.message);
        return;
      }
      if (data.kind !== 'result') return;

      worker.terminate();
      activeWorker = null;
      const result = data.result;
      result.metrics.preset = preset;
      try {
        plotQueue = plotQueue.catch(() => {}).then(() => renderPlots(result));
        await plotQueue;
        lastResult = result;
        renderMetrics(result.metrics);
        renderApdTable(result);
        const fresh = id === runId && isSameConfig(getInputValues(false), result.inputs);
        setResultsStale(!fresh);
        if (id !== runId) return;
        setStatus('ready', fresh
          ? `Simulation complete in ${result.runtimeMs.toFixed(0)} ms.`
          : 'Simulation complete. Parameters changed while it was running.');
      } catch (error) {
        if (id === runId) finishWithError(error.message || 'Plot rendering failed.');
      } finally {
        if (id === runId) updateRunningState(false);
      }
    };

    worker.onerror = (event) => {
      if (id !== runId) return;
      console.error(event);
      finishWithError('The simulation worker failed to load or run.');
    };

    worker.postMessage({ kind: 'run', id, inputs: inputValues });
  }

  function cancelActiveWorker() {
    runId += 1;
    if (activeWorker) activeWorker.terminate();
    activeWorker = null;
  }

  function finishWithError(message) {
    cancelActiveWorker();
    setResultsStale(true);
    updateRunningState(false);
    setStatus('error', message);
  }

  function renderPlots(result) {
    if (typeof Plotly === 'undefined') {
      throw new Error('Plotly is not loaded from assets/js/plotly.js.');
    }
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
        name: 'APD by detected event',
        xaxis: 'x3',
        yaxis: 'y3',
        marker: { color: '#dc2626', size: 7 },
        line: { color: '#fca5a5', width: 1.5 },
        hovertemplate: 'Detected APD: %{x}<br>APD: %{y:.2f} ms<extra></extra>',
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
      margin: { l: 64, r: 24, t: 100, b: 100 },
      showlegend: false,
      height: plotHeight(),
      grid: {
        rows: 4,
        columns: 1,
        pattern: 'independent',
        roworder: 'top to bottom',
        ygap: 0.15,
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
        title: 'Detected APD number',
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
          text: `Preset: ${formatPresetName(result.metrics.preset)} · PCL: ${formatNumber(result.inputs.pcl, 0)} ms · Beats: ${formatNumber(result.inputs.beats, 0)}`,
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
      return Plotly.react(ui.plot, traces, layout, config);
    } else {
      return Plotly.newPlot(ui.plot, traces, layout, config);
    }
  }

  function plotHeight() {
    return window.matchMedia('(max-width: 720px)').matches ? 1240 : 960;
  }

  function renderMetrics(metrics) {
    ui.metrics.apdMin.textContent = formatMetric(metrics.minApd, ' ms');
    ui.metrics.apdMax.textContent = formatMetric(metrics.maxApd, ' ms');
    ui.metrics.apdRange.textContent = formatMetric(metrics.apdRange, ' ms');
    ui.metrics.beatCount.textContent = Number.isFinite(metrics.beatCount) ? String(metrics.beatCount) : '—';
    ui.metrics.beatNote.textContent = `${metrics.deliveredStimuli} stimuli delivered; ${metrics.beatCount} complete APDs`;
    ui.metrics.preset.textContent = formatPresetName(metrics.preset);

    if (Number.isFinite(metrics.minX) && Number.isFinite(metrics.maxX)) {
      ui.metrics.presetNote.textContent = `Initial X range: ${metrics.minX.toFixed(4)} to ${metrics.maxX.toFixed(4)}`;
    } else {
      ui.metrics.presetNote.textContent = 'No complete APDs detected in this run';
    }
  }

  function renderApdTable(result) {
    ui.apdTableBody.replaceChildren();
    result.apds.forEach((apd, index) => {
      const row = document.createElement('tr');
      [index + 1, apd.toFixed(2), result.xsInit[index].toFixed(4)].forEach((value) => {
        const cell = document.createElement('td');
        cell.textContent = String(value);
        row.appendChild(cell);
      });
      ui.apdTableBody.appendChild(row);
    });
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

  function getInputValues(validate = false) {
    const values = {};
    INPUT_IDS.forEach((id) => {
      const input = ui.inputs[id];
      const value = input.valueAsNumber;
      if (validate && !input.checkValidity()) {
        const issue = input.validity;
        let message = `${FIELD_LABELS[id]} must be a valid number.`;
        if (issue.valueMissing) message = `Enter ${FIELD_LABELS[id]}.`;
        else if (issue.rangeUnderflow || issue.rangeOverflow) {
          message = `${FIELD_LABELS[id]} must be between ${input.min} and ${input.max}.`;
        } else if (issue.stepMismatch) {
          message = `${FIELD_LABELS[id]} must use increments of ${input.step}.`;
        }
        const error = new Error(message);
        error.input = input;
        throw error;
      }
      values[id] = value;
    });
    return values;
  }

  function clearInputError() {
    ui.inputError.hidden = true;
    ui.inputError.textContent = '';
    INPUT_IDS.forEach((id) => {
      ui.inputs[id].removeAttribute('aria-invalid');
      ui.inputs[id].removeAttribute('aria-describedby');
    });
  }

  function showInputError(error) {
    const input = error.input;
    if (input) {
      input.setAttribute('aria-invalid', 'true');
      input.setAttribute('aria-describedby', 'inputError');
      input.focus();
    }
    ui.inputError.textContent = error.message;
    ui.inputError.hidden = false;
    setStatus('error', error.message);
  }

  function setResultsStale(stale) {
    resultsStale = stale;
    ui.resultState.textContent = lastResult
      ? (stale ? 'Results show the previous run. Run again to refresh them.' : 'Results match the current parameters.')
      : 'No completed simulation yet.';
    ui.resultState.classList.toggle('stale', stale);
    ui.exportBtn.disabled = stale || !lastResult || !!activeWorker;
  }

  function applyPreset(presetName) {
    const preset = PRESETS[presetName];
    if (!preset) return;

    setFormValues(preset);
    setActivePreset(presetName);
    clearInputError();
    setResultsStale(true);
    queueRun(true);
  }

  function setActivePreset(presetName) {
    currentPreset = presetName;

    Object.entries(ui.presetButtons).forEach(([name, button]) => {
      button.classList.remove('btn-pill-active', 'btn-secondary', 'btn-ghost');
      button.setAttribute('aria-pressed', String(name === presetName));

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
      button.setAttribute('aria-pressed', 'false');
    });
  }

  function exportPlotImage() {
    if (!ui.plot || !ui.plot.data || typeof Plotly === 'undefined' || resultsStale || activeWorker) {
      setStatus('error', 'Run the simulation with the current parameters before exporting.');
      return;
    }

    try {
      Promise.resolve(Plotly.downloadImage(ui.plot, {
        format: 'png',
        filename: 'lr1-simulation',
        width: 1600,
        height: 1400,
        scale: 2,
      })).catch((error) => setStatus('error', `Plot download failed: ${error.message}`));
    } catch (error) {
      setStatus('error', `Plot download failed: ${error.message}`);
    }
  }

  function updateRunningState(running) {
    ui.runBtn.textContent = running ? 'Run again' : 'Run simulation';
    ui.exportBtn.disabled = running || !lastResult || resultsStale;
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
