/* global LR1CellIto */
'use strict';

importScripts('LR1CellIto.js');

const DT_MS = 0.01;
const WARMUP_MS = 10000;
const CONDITIONING_BEATS = 10;
const MAX_FIRST_STIM_MS = 100;
const SAMPLE_INTERVAL_MS = 5;
const MAX_PLOT_SAMPLES = 25000;
const LIMITS = {
  inafac: [0, 5], itofac: [0, 5], icalfac: [0, 5], ikfac: [0, 5],
  ikifac: [0, 5], tauXfac: [0.1, 25], yshift: [-40, 40],
  beats: [1, 250], pcl: [50, 5000],
};

self.onmessage = ({ data }) => {
  if (!data || data.kind !== 'run') return;

  try {
    self.postMessage({ kind: 'result', id: data.id, result: runSimulation(data.inputs, data.id) });
  } catch (error) {
    self.postMessage({
      kind: 'error',
      id: data.id,
      message: error instanceof Error ? error.message : 'Simulation failed.',
    });
  }
};

function validateInputs(inputs) {
  for (const [key, [min, max]] of Object.entries(LIMITS)) {
    const value = inputs[key];
    if (!Number.isFinite(value) || value < min || value > max) {
      throw new Error(`${key} must be between ${min} and ${max}.`);
    }
  }
  if (!Number.isInteger(inputs.beats) || !Number.isInteger(inputs.pcl)) {
    throw new Error('Beats and PCL must be whole numbers.');
  }
}

function checkCell(cell) {
  const state = [cell.v, cell.cai, cell.m, cell.h, cell.j, cell.xr,
    cell.d, cell.f, cell.xtos, cell.ytos, cell.xtof, cell.ytof];
  if (cell.cai <= 0 || !state.every(Number.isFinite)) {
    throw new Error('The model became numerically unstable for these parameters.');
  }
}

function runSimulation(inputs, id) {
  validateInputs(inputs);
  const startedAt = performance.now();
  const cell = new LR1CellIto();
  for (const key of ['inafac', 'itofac', 'icalfac', 'ikfac', 'ikifac', 'tauXfac', 'yshift']) {
    cell[key] = inputs[key];
  }

  const dt = cell.dt;
  if (dt !== DT_MS) throw new Error('The model time step has changed; update the pacing schedule.');
  const warmupSteps = Math.ceil(WARMUP_MS / dt);
  // Keep the first stimulus inside the first cycle, including short PCLs.
  const firstStimStep = Math.round(Math.min(MAX_FIRST_STIM_MS, inputs.pcl / 2) / dt);
  const pclSteps = Math.round(inputs.pcl / dt);
  const stimSteps = Math.round(cell.stimduration / dt);
  const conditioningSteps = CONDITIONING_BEATS * pclSteps;
  const recordingSteps = firstStimStep + inputs.beats * pclSteps;
  // Keep long traces small enough for Plotly; APD detection still uses every step.
  const sampleSteps = Math.max(
    Math.round(SAMPLE_INTERVAL_MS / dt),
    Math.ceil(recordingSteps / MAX_PLOT_SAMPLES),
  );
  const totalSteps = warmupSteps + conditioningSteps + recordingSteps;
  const progressStride = Math.max(100000, Math.ceil(totalSteps / 100));
  let stepsDone = 0;
  let nextProgressStep = progressStride;
  let nextCheckStep = 1000;

  function step(stimulus) {
    cell.stepdt(stimulus);
    stepsDone += 1;
    if (stepsDone >= nextCheckStep) {
      checkCell(cell);
      nextCheckStep += 1000;
      if (stepsDone >= nextProgressStep) {
        self.postMessage({ kind: 'progress', id, percent: Math.floor(100 * stepsDone / totalSteps) });
        nextProgressStep += progressStride;
      }
    }
  }

  function stimulusAt(stepNumber) {
    return stepNumber >= firstStimStep &&
      (stepNumber - firstStimStep) % pclSteps < stimSteps ? cell.stimulus : 0;
  }

  for (let i = 0; i < warmupSteps; i += 1) step(0);
  for (let i = 1; i < conditioningSteps; i += 1) step(stimulusAt(i));

  const ts = [0];
  const vs = [cell.v];
  const xs = [cell.xr];
  const apds = [];
  const xsInit = [];
  const apdCounter = [];
  const stimulusTimes = [];
  let startApdTime = null;
  let initialX = null;
  let minApd = Infinity;
  let maxApd = -Infinity;
  let minX = Infinity;
  let maxX = -Infinity;

  for (let i = 1; i < recordingSteps; i += 1) {
    const oldV = cell.v;
    const stimulus = stimulusAt(i);
    if (stimulus && (i - firstStimStep) % pclSteps === 0) {
      stimulusTimes.push(i * dt);
    }
    step(stimulus);
    const newV = cell.v;
    const timeMs = i * dt;

    if (oldV <= cell.threshold && newV > cell.threshold) {
      startApdTime = timeMs;
      initialX = cell.xr;
    } else if (oldV >= cell.threshold && newV < cell.threshold && startApdTime !== null) {
      const apd = timeMs - startApdTime;
      apds.push(apd);
      xsInit.push(initialX);
      apdCounter.push(apds.length);
      minApd = Math.min(minApd, apd);
      maxApd = Math.max(maxApd, apd);
      minX = Math.min(minX, initialX);
      maxX = Math.max(maxX, initialX);
      startApdTime = null;
      initialX = null;
    }

    if (i % sampleSteps === 0) {
      ts.push(timeMs / 1000);
      vs.push(newV);
      xs.push(cell.xr);
    }
  }

  checkCell(cell);
  if ((recordingSteps - 1) % sampleSteps !== 0) {
    ts.push((recordingSteps - 1) * dt / 1000);
    vs.push(cell.v);
    xs.push(cell.xr);
  }
  const hasApd = apds.length > 0;
  return {
    ts, vs, xs, apds, xsInit, apdCounter, stimulusTimes,
    runtimeMs: performance.now() - startedAt,
    sampleIntervalMs: sampleSteps * dt,
    inputs,
    metrics: {
      minApd: hasApd ? minApd : NaN,
      maxApd: hasApd ? maxApd : NaN,
      apdRange: hasApd ? maxApd - minApd : NaN,
      beatCount: apds.length,
      deliveredStimuli: stimulusTimes.length,
      minX: hasApd ? minX : NaN,
      maxX: hasApd ? maxX : NaN,
    },
  };
}
