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
  if (!data || !['run', 'sweep'].includes(data.kind)) return;

  try {
    if (data.kind === 'sweep') {
      self.postMessage({ kind: 'sweep-result', id: data.id, result: runSweep(data) });
    } else {
      self.postMessage({ kind: 'result', id: data.id, result: runSimulation(data.inputs, data.id) });
    }
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

function runSimulation(inputs, id, summaryOnly = false) {
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
      if (!summaryOnly && stepsDone >= nextProgressStep) {
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

  const ts = summaryOnly ? [] : [0];
  const vs = summaryOnly ? [] : [cell.v];
  const xs = summaryOnly ? [] : [cell.xr];
  const icalTrace = summaryOnly ? [] : [cell.ical];
  const ikrTrace = summaryOnly ? [] : [cell.ikr];
  const netTrace = summaryOnly ? [] : [cell.ina + cell.ical + cell.ito + cell.ikr + cell.ik1 + cell.ipk + cell.ib];
  const apds = [];
  const xsInit = [];
  const apdCounter = [];
  const apdStartTimes = [];
  const apdEndTimes = [];
  const eadCounts = [];
  const stimulusTimes = [];
  let startApdTime = null;
  let initialX = null;
  let eadCount = 0;
  let eadTrough = Infinity;
  let eadPeak = -Infinity;
  let eadRising = false;
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
      eadCount = 0;
      eadTrough = Infinity;
      eadPeak = -Infinity;
      eadRising = false;
    } else if (oldV >= cell.threshold && newV < cell.threshold && startApdTime !== null) {
      const apd = timeMs - startApdTime;
      apds.push(apd);
      xsInit.push(initialX);
      apdCounter.push(apds.length);
      apdStartTimes.push(startApdTime);
      apdEndTimes.push(timeMs);
      eadCounts.push(eadCount);
      minApd = Math.min(minApd, apd);
      maxApd = Math.max(maxApd, apd);
      minX = Math.min(minX, initialX);
      maxX = Math.max(maxX, initialX);
      startApdTime = null;
      initialX = null;
    } else if (startApdTime !== null && timeMs - startApdTime >= 100) {
      // A secondary rise of at least 5 mV from a trough below -15 mV.
      // Wait for a 5 mV fall before allowing another EAD on the same AP.
      if (!eadRising) {
        eadTrough = Math.min(eadTrough, newV);
        if (eadTrough <= -15 && newV - eadTrough >= 5) {
          eadCount += 1;
          eadRising = true;
          eadPeak = newV;
        }
      } else {
        eadPeak = Math.max(eadPeak, newV);
        if (eadPeak - newV >= 5) {
          eadRising = false;
          eadTrough = newV;
        }
      }
    }

    if (!summaryOnly && i % sampleSteps === 0) {
      ts.push(timeMs / 1000);
      vs.push(newV);
      xs.push(cell.xr);
      icalTrace.push(cell.ical);
      ikrTrace.push(cell.ikr);
      netTrace.push(cell.ina + cell.ical + cell.ito + cell.ikr + cell.ik1 + cell.ipk + cell.ib);
    }
  }

  checkCell(cell);
  if (!summaryOnly && (recordingSteps - 1) % sampleSteps !== 0) {
    ts.push((recordingSteps - 1) * dt / 1000);
    vs.push(cell.v);
    xs.push(cell.xr);
    icalTrace.push(cell.ical);
    ikrTrace.push(cell.ikr);
    netTrace.push(cell.ina + cell.ical + cell.ito + cell.ikr + cell.ik1 + cell.ipk + cell.ib);
  }
  const hasApd = apds.length > 0;
  return {
    ts, vs, xs, icalTrace, ikrTrace, netTrace,
    apds, xsInit, apdCounter, apdStartTimes, apdEndTimes, eadCounts, stimulusTimes,
    runtimeMs: performance.now() - startedAt,
    sampleIntervalMs: sampleSteps * dt,
    inputs,
    metrics: {
      minApd: hasApd ? minApd : NaN,
      maxApd: hasApd ? maxApd : NaN,
      apdRange: hasApd ? maxApd - minApd : NaN,
      beatCount: apds.length,
      deliveredStimuli: stimulusTimes.length,
      totalEads: eadCounts.reduce((sum, count) => sum + count, 0),
      beatsWithEads: eadCounts.filter((count) => count > 0).length,
      minX: hasApd ? minX : NaN,
      maxX: hasApd ? maxX : NaN,
    },
  };
}

function runSweep({ inputs, start, end, interval, beats, tail, id }) {
  if (![start, end, interval, beats, tail].every(Number.isInteger) ||
      start < LIMITS.pcl[0] || end > LIMITS.pcl[1] || end < start ||
      interval < 1 || beats < 1 || beats > LIMITS.beats[1] ||
      tail < 1 || tail > beats) {
    throw new Error('Enter valid whole-number PCL bounds, step, beats, and late-beat count.');
  }
  const pclValues = [];
  for (let pcl = start; pcl <= end; pcl += interval) pclValues.push(pcl);
  if (pclValues.length < 2 || pclValues.length > 25) {
    throw new Error('The sweep needs 2 to 25 PCL points.');
  }
  const estimatedSteps = pclValues.reduce((sum, pcl) =>
    sum + Math.ceil((WARMUP_MS + (CONDITIONING_BEATS + beats + 1) * pcl) / DT_MS), 0);
  if (estimatedSteps > 100000000) {
    throw new Error('This sweep is too large. Reduce the points or beats per PCL.');
  }
  const points = [];
  const startedAt = performance.now();
  pclValues.forEach((pcl, index) => {
    const result = runSimulation({ ...inputs, pcl, beats }, id, true);
    const cutoffMs = result.apdEndTimes.length
      ? result.stimulusTimes[result.stimulusTimes.length - tail]
      : Infinity;
    const lateApds = result.apds.filter((_, i) => result.apdEndTimes[i] >= cutoffMs);
    points.push({ pcl, apds: lateApds, completeApds: lateApds.length,
      deliveredStimuli: tail, captureRatio: lateApds.length / tail });
    self.postMessage({ kind: 'sweep-progress', id, completed: index + 1, total: pclValues.length });
  });
  return { points, inputs, start, end, interval, beats, tail, runtimeMs: performance.now() - startedAt };
}
