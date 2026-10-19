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
  let sweepWorker = null;
  let sweepRunId = 0;
  let lastSweep = null;

  function init() {
    ui = {
      form: document.getElementById('controlsForm'),
      plot: document.getElementById('myDiv1'),
      memoryPlot: document.getElementById('memoryPlot'),
      detailPlot: document.getElementById('detailPlot'),
      apdSelect: document.getElementById('apdSelect'),
      selectedEadCount: document.getElementById('selectedEadCount'),
      sweepPlot: document.getElementById('sweepPlot'),
      sweepTableBody: document.getElementById('sweepTableBody'),
      sweepBtn: document.getElementById('sweepBtn'),
      sweepCancelBtn: document.getElementById('sweepCancelBtn'),
      sweepStatus: document.getElementById('sweepStatus'),
      sweepInputs: Object.fromEntries(['Start', 'End', 'Step', 'Beats', 'Tail'].map((name) =>
        [name.toLowerCase(), document.getElementById(`sweep${name}`)])),
      status: document.getElementById('statusIndicator'),
      statusText: document.getElementById('statusText'),
      presetDetail: document.getElementById('presetDetail'),
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
        eadCount: document.getElementById('metricEadCount'),
        eadNote: document.getElementById('metricEadNote'),
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
      setSweepDefaults('default');
      markSweepStale();
      setActivePreset('default');
      setResultsStale(true);
      queueRun(true);
    });

    ui.exportBtn.addEventListener('click', exportPlotImage);
    ui.apdSelect.addEventListener('change', () => {
      if (lastResult) renderDetailPlot(lastResult, Number(ui.apdSelect.value)).catch(console.error);
    });
    ui.sweepBtn.addEventListener('click', startSweep);
    ui.sweepCancelBtn.addEventListener('click', () => cancelSweep('Sweep canceled.'));
    Object.values(ui.sweepInputs).forEach((input) => {
      input.addEventListener('input', () => {
        if (lastSweep && !sweepWorker) ui.sweepStatus.textContent = 'Sweep shows the previous range. Run it again to refresh.';
      });
    });
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
      if (ui.memoryPlot.data) Plotly.relayout(ui.memoryPlot, { height: window.innerWidth <= 720 ? 790 : 530 }).catch(console.error);
      if (ui.detailPlot.data) Plotly.relayout(ui.detailPlot, { height: window.innerWidth <= 720 ? 670 : 530 }).catch(console.error);
    }, 150));
  }

  function handleInputChange() {
    clearInputError();
    detectPresetFromCurrentValues();
    setResultsStale(true);
    markSweepStale();

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
        plotQueue = plotQueue.catch(() => {}).then(async () => {
          await renderPlots(result);
          await renderMemoryPlot(result);
          populateApdSelect(result);
          await renderDetailPlot(result, Number(ui.apdSelect.value));
        });
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

  function isPublishedChaosSetting(inputs) {
    return Math.abs(inputs.itofac - 1.05) < 1e-9 &&
      Math.abs(inputs.icalfac - 1.15) < 1e-9 &&
      Math.abs(inputs.ikifac - 2.2) < 1e-9 &&
      Math.abs(inputs.tauXfac - 5) < 1e-9 &&
      Math.abs(inputs.yshift - 8) < 1e-9 &&
      Math.abs(inputs.inafac - 1) < 1e-9 && Math.abs(inputs.ikfac - 1) < 1e-9;
  }

  function publishedMemoryMap(apd, pcl) {
    // Illustrative map parameters from Landaw and Qu, PRE 97, 042414 (2018), Fig. 10.
    const x = 0.07 * Math.pow((350 - apd) / (apd - 50), 1 / 25);
    const nextX = (0.6 - (0.6 - x) * Math.exp(-apd / 3000)) * Math.exp(-(pcl - apd) / 1000);
    return 50 + 300 / (1 + Math.pow(nextX / 0.07, 25));
  }

  function plotInto(element, traces, layout, filename) {
    const config = { responsive: true, displaylogo: false,
      toImageButtonOptions: { format: 'png', filename, scale: 2 } };
    return element.data ? Plotly.react(element, traces, layout, config)
      : Plotly.newPlot(element, traces, layout, config);
  }

  function renderMemoryPlot(result) {
    const mobile = window.innerWidth <= 720;
    const apds = result.apds;
    const bound = Math.max(400, ...apds, 0);
    const traces = [
      { x: result.xsInit, y: apds, type: 'scatter', mode: 'markers', name: 'Measured APD vs X',
        marker: { color: '#059669', size: 7 },
        hovertemplate: 'Initial X: %{x:.4f}<br>APD: %{y:.2f} ms<extra></extra>' },
      { x: [0, bound], y: [0, bound], type: 'scatter', mode: 'lines', name: 'Unchanged APD',
        xaxis: 'x2', yaxis: 'y2', line: { color: '#94a3b8', dash: 'dash', width: 1.5 }, hoverinfo: 'skip' },
      { x: apds.slice(0, -1), y: apds.slice(1), type: 'scatter', mode: 'markers',
        name: 'Measured next APD', xaxis: 'x2', yaxis: 'y2',
        marker: { color: '#dc2626', size: 7 },
        hovertemplate: 'APD n: %{x:.2f} ms<br>APD n+1: %{y:.2f} ms<extra></extra>' },
    ];
    if (isPublishedChaosSetting(result.inputs) &&
        result.metrics.beatCount === result.metrics.deliveredStimuli) {
      const xs = Array.from({ length: 301 }, (_, i) => 50.1 + i * 0.999);
      traces.push({ x: xs, y: xs.map((apd) => publishedMemoryMap(apd, result.inputs.pcl)),
        type: 'scatter', mode: 'lines', name: 'Published X-memory approximation',
        xaxis: 'x2', yaxis: 'y2', line: { color: '#7c3aed', width: 2.3 },
        hovertemplate: 'APD n: %{x:.1f} ms<br>Map APD n+1: %{y:.1f} ms<extra></extra>' });
    }
    const layout = {
      paper_bgcolor: '#fff', plot_bgcolor: '#fbfdff',
      height: mobile ? 790 : 530,
      margin: { l: 64, r: 25, t: 60, b: 95 },
      showlegend: true, legend: { orientation: 'h', y: -0.22, x: 0 },
      grid: { rows: mobile ? 2 : 1, columns: mobile ? 1 : 2,
        pattern: 'independent', roworder: 'top to bottom', xgap: 0.16, ygap: 0.25 },
      xaxis: { title: 'Initial X', gridcolor: '#e2e8f0', range: computeXScatterRange(result.metrics.minX, result.metrics.maxX) },
      yaxis: { title: 'APD (ms)', gridcolor: '#e2e8f0', range: computeApdRange(apds) },
      xaxis2: { title: 'APD n (ms)', gridcolor: '#e2e8f0', range: [0, bound] },
      yaxis2: { title: 'APD n+1 (ms)', gridcolor: '#e2e8f0', range: [0, bound] },
    };
    return plotInto(ui.memoryPlot, traces, layout, 'lr1-x-memory');
  }

  function populateApdSelect(result) {
    ui.apdSelect.replaceChildren();
    if (!result.apds.length) {
      const option = new Option('No complete APDs', '');
      ui.apdSelect.add(option);
      ui.apdSelect.disabled = true;
      return;
    }
    ui.apdSelect.disabled = false;
    result.apds.forEach((apd, index) => {
      ui.apdSelect.add(new Option(`#${index + 1} · ${apd.toFixed(1)} ms · ${result.eadCounts[index]} EADs`, String(index)));
    });
    const firstEad = result.eadCounts.findIndex((count) => count > 0);
    ui.apdSelect.value = String(firstEad >= 0 ? firstEad : result.apds.length - 1);
  }

  function renderDetailPlot(result, index) {
    if (!Number.isInteger(index) || index < 0 || index >= result.apds.length) {
      ui.selectedEadCount.textContent = 'No complete APD to inspect.';
      if (ui.detailPlot.data) Plotly.purge(ui.detailPlot);
      return Promise.resolve();
    }
    const start = result.apdStartTimes[index];
    const end = result.apdEndTimes[index];
    const relativeMs = [];
    const voltage = [];
    const ical = [];
    const ikr = [];
    const net = [];
    for (let i = 0; i < result.ts.length; i += 1) {
      const time = result.ts[i] * 1000;
      if (time < start - 50 || time > end + 100) continue;
      relativeMs.push(time - start);
      voltage.push(result.vs[i]);
      ical.push(result.icalTrace[i]);
      ikr.push(result.ikrTrace[i]);
      net.push(result.netTrace[i]);
    }
    const count = result.eadCounts[index];
    ui.selectedEadCount.textContent = `${count} detected EAD${count === 1 ? '' : 's'} · APD ${result.apds[index].toFixed(1)} ms`;
    const currentStart = Math.min(100, Math.max(10, result.apds[index] * 0.2));
    const currentTimes = [];
    const plateauIcal = [];
    const plateauIkr = [];
    const plateauNet = [];
    relativeMs.forEach((time, i) => {
      if (time >= currentStart) {
        currentTimes.push(time);
        plateauIcal.push(ical[i]);
        plateauIkr.push(ikr[i]);
        plateauNet.push(net[i]);
      }
    });
    const traces = [
      { x: relativeMs, y: voltage, type: 'scatter', mode: 'lines', name: 'Voltage',
        line: { color: '#2563eb', width: 2.3 },
        hovertemplate: 'Time from upstroke: %{x:.1f} ms<br>Voltage: %{y:.2f} mV<extra></extra>' },
      { x: currentTimes, y: plateauIcal, type: 'scatter', mode: 'lines', name: 'I Ca,L', xaxis: 'x2', yaxis: 'y2',
        line: { color: '#dc2626', width: 2 }, hovertemplate: 'Time: %{x:.1f} ms<br>I Ca,L: %{y:.3f}<extra></extra>' },
      { x: currentTimes, y: plateauIkr, type: 'scatter', mode: 'lines', name: 'I Kr', xaxis: 'x2', yaxis: 'y2',
        line: { color: '#7c3aed', width: 2 }, hovertemplate: 'Time: %{x:.1f} ms<br>I Kr: %{y:.3f}<extra></extra>' },
      { x: currentTimes, y: plateauNet, type: 'scatter', mode: 'lines', name: 'Net ionic current', xaxis: 'x2', yaxis: 'y2',
        line: { color: '#0f766e', width: 2 }, hovertemplate: 'Time: %{x:.1f} ms<br>Net current: %{y:.3f}<extra></extra>' },
    ];
    const layout = {
      paper_bgcolor: '#fff', plot_bgcolor: '#fbfdff', height: window.innerWidth <= 720 ? 670 : 530,
      margin: { l: 65, r: 25, t: 45, b: 90 }, showlegend: true,
      legend: { orientation: 'h', y: -0.22, x: 0 },
      grid: { rows: 2, columns: 1, pattern: 'independent', roworder: 'top to bottom', ygap: 0.18 },
      xaxis: { title: 'Time from upstroke (ms)', gridcolor: '#e2e8f0' },
      yaxis: { title: 'Voltage (mV)', gridcolor: '#e2e8f0' },
      xaxis2: { title: 'Time from upstroke (ms)', gridcolor: '#e2e8f0', matches: 'x' },
      yaxis2: { title: 'Current (model units)', gridcolor: '#e2e8f0' },
    };
    return plotInto(ui.detailPlot, traces, layout, 'lr1-action-potential-detail');
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
    ui.metrics.eadCount.textContent = String(metrics.totalEads);
    ui.metrics.eadNote.textContent = `${metrics.beatsWithEads} complete APDs with EADs`;
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
      [index + 1, apd.toFixed(2), result.xsInit[index].toFixed(4), result.eadCounts[index]].forEach((value) => {
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
    setSweepDefaults(presetName);
    markSweepStale();
    setActivePreset(presetName);
    clearInputError();
    setResultsStale(true);
    queueRun(true);
  }

  function setSweepDefaults(presetName) {
    const presets = {
      default: [350, 650, 30, 30, 8],
      chaos: [330, 530, 20, 40, 10],
      eads: [1500, 1650, 15, 24, 6],
    };
    const values = presets[presetName];
    if (!values) return;
    ['start', 'end', 'step', 'beats', 'tail'].forEach((key, index) => {
      ui.sweepInputs[key].value = String(values[index]);
    });
  }

  function sweepSettingsChanged(previous) {
    const current = getInputValues(false);
    return ['inafac', 'itofac', 'icalfac', 'ikfac', 'ikifac', 'tauXfac', 'yshift']
      .some((key) => !Number.isFinite(current[key]) || Math.abs(current[key] - previous[key]) > 1e-9);
  }

  function markSweepStale() {
    if (lastSweep && !sweepWorker && sweepSettingsChanged(lastSweep.inputs)) {
      ui.sweepStatus.textContent = 'Sweep shows earlier model settings. Run it again to refresh.';
    }
  }

  function cancelSweep(message) {
    sweepRunId += 1;
    if (sweepWorker) sweepWorker.terminate();
    sweepWorker = null;
    ui.sweepBtn.disabled = false;
    ui.sweepCancelBtn.hidden = true;
    ui.sweepStatus.textContent = message;
    ui.sweepStatus.classList.remove('error');
  }

  function startSweep() {
    let inputs;
    try {
      inputs = getInputValues(true);
      clearInputError();
    } catch (error) {
      showInputError(error);
      return;
    }
    const values = Object.fromEntries(Object.entries(ui.sweepInputs).map(([key, input]) => [key, input.valueAsNumber]));
    const invalid = Object.entries(ui.sweepInputs).find(([, input]) => !input.checkValidity());
    if (invalid || values.end < values.start || values.tail > values.beats ||
        Math.floor((values.end - values.start) / values.step) + 1 > 25 ||
        Math.floor((values.end - values.start) / values.step) + 1 < 2) {
      ui.sweepStatus.textContent = 'Enter 2–25 PCL points and a late-cycle count no greater than beats per PCL.';
      ui.sweepStatus.classList.add('error');
      if (invalid) invalid[1].focus();
      return;
    }
    if (sweepWorker) cancelSweep('Starting a new sweep…');
    const id = ++sweepRunId;
    let worker;
    try {
      worker = new Worker('assets/js/simulation-worker.js');
    } catch (error) {
      ui.sweepStatus.textContent = `Could not start sweep worker: ${error.message}`;
      ui.sweepStatus.classList.add('error');
      return;
    }
    sweepWorker = worker;
    ui.sweepBtn.disabled = true;
    ui.sweepCancelBtn.hidden = false;
    ui.sweepStatus.classList.remove('error');
    ui.sweepStatus.textContent = 'Running independent simulations…';
    worker.onmessage = async ({ data }) => {
      if (!data || data.id !== sweepRunId || worker !== sweepWorker) return;
      if (data.kind === 'sweep-progress') {
        ui.sweepStatus.textContent = `Completed ${data.completed} of ${data.total} PCL points…`;
        return;
      }
      worker.terminate();
      sweepWorker = null;
      ui.sweepBtn.disabled = false;
      ui.sweepCancelBtn.hidden = true;
      if (data.kind === 'error') {
        ui.sweepStatus.textContent = data.message;
        ui.sweepStatus.classList.add('error');
        return;
      }
      if (data.kind !== 'sweep-result') return;
      try {
        await renderSweepPlot(data.result);
        lastSweep = data.result;
        renderSweepTable(data.result);
        const stale = sweepSettingsChanged(data.result.inputs);
        ui.sweepStatus.textContent = `Sweep complete in ${(data.result.runtimeMs / 1000).toFixed(1)} s.${stale ? ' Settings changed during the sweep; these results use the earlier settings.' : ''}`;
      } catch (error) {
        ui.sweepStatus.textContent = `Could not draw sweep: ${error.message}`;
        ui.sweepStatus.classList.add('error');
      }
    };
    worker.onerror = () => {
      if (worker !== sweepWorker) return;
      worker.terminate();
      sweepWorker = null;
      ui.sweepBtn.disabled = false;
      ui.sweepCancelBtn.hidden = true;
      ui.sweepStatus.textContent = 'The sweep worker failed to load or run.';
      ui.sweepStatus.classList.add('error');
    };
    worker.postMessage({ kind: 'sweep', id, inputs, start: values.start, end: values.end,
      interval: values.step, beats: values.beats, tail: values.tail });
  }

  function renderSweepPlot(result) {
    ui.sweepPlot.hidden = false;
    const xApd = [];
    const yApd = [];
    result.points.forEach((point) => point.apds.forEach((apd) => {
      xApd.push(point.pcl);
      yApd.push(apd);
    }));
    const traces = [
      { x: xApd, y: yApd, type: 'scatter', mode: 'markers', name: 'Late APDs',
        marker: { color: '#dc2626', size: 6, opacity: 0.7 },
        hovertemplate: 'PCL: %{x} ms<br>APD: %{y:.2f} ms<extra></extra>' },
      { x: result.points.map((point) => point.pcl),
        y: result.points.map((point) => point.captureRatio),
        xaxis: 'x2', yaxis: 'y2', type: 'scatter', mode: 'lines+markers', name: 'Complete APDs / stimuli',
        line: { color: '#059669', width: 2 }, marker: { size: 7 },
        hovertemplate: 'PCL: %{x} ms<br>Complete APDs / stimuli: %{y:.2f}<extra></extra>' },
    ];
    const layout = { paper_bgcolor: '#fff', plot_bgcolor: '#fbfdff',
      height: window.innerWidth <= 720 ? 620 : 540,
      margin: { l: 65, r: 25, t: 35, b: 65 }, showlegend: false,
      grid: { rows: 2, columns: 1, pattern: 'independent', roworder: 'top to bottom', ygap: 0.18 },
      xaxis: { title: 'PCL (ms)', gridcolor: '#e2e8f0', range: [result.start - result.interval, result.end + result.interval] },
      yaxis: { title: 'Late APDs (ms)', gridcolor: '#e2e8f0' },
      xaxis2: { title: 'PCL (ms)', gridcolor: '#e2e8f0', matches: 'x' },
      yaxis2: { title: 'Complete APDs / stimuli', gridcolor: '#e2e8f0', range: [0, 1.1] },
    };
    return plotInto(ui.sweepPlot, traces, layout, 'lr1-pcl-sweep');
  }

  function renderSweepTable(result) {
    ui.sweepTableBody.replaceChildren();
    result.points.forEach((point) => {
      const row = document.createElement('tr');
      [point.pcl, point.apds.map((apd) => apd.toFixed(1)).join(', ') || '—',
        `${point.completeApds}/${point.deliveredStimuli} (${(100 * point.captureRatio).toFixed(0)}%)`]
        .forEach((value) => {
          const cell = document.createElement('td');
          cell.textContent = String(value);
          row.appendChild(cell);
        });
      ui.sweepTableBody.appendChild(row);
    });
  }

  function setActivePreset(presetName) {
    currentPreset = presetName;
    const descriptions = {
      default: 'Default: Ito off; baseline Gsi = 0.09 and GK1 = 0.6047 mS/cm²; X time constant ×1.',
      chaos: 'Chaos: gto,f = 0.21, Gsi = 0.1035, GK1 = 1.33034 mS/cm²; X time constant ×5; Ito inactivation shifted +8 mV.',
      eads: 'EADs: Ito off; baseline Gsi and GK1; X time constant ×10; preset PCL 1575 ms.',
    };
    ui.presetDetail.textContent = descriptions[presetName] || 'Custom: values are set by the controls below.';

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
    ui.presetDetail.textContent = 'Custom: values are set by the controls below.';
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
