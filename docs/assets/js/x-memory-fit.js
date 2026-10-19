(function () {
  'use strict';

  // The paper's Hill parameters are illustrative (PRE 97, 042414, Fig. 10).
  // Its Fig. 9 comparison instead interpolates a measured APD-versus-X curve.
  function publishedReturnMap(apd, pcl) {
    const bounded = Math.min(349.999, Math.max(50.001, apd));
    const x = 0.07 * Math.pow((350 - bounded) / (bounded - 50), 1 / 25);
    const nextX = (0.6 - (0.6 - x) * Math.exp(-apd / 3000)) * Math.exp(-(pcl - apd) / 1000);
    return 50 + 300 / (1 + Math.pow(nextX / 0.07, 25));
  }

  function plateauFeatures(x, apd) {
    const t = apd / 350;
    return [x, t, t * t, t * t * t, x * t, x * t * t, x * t * t * t];
  }

  function dot(a, b) {
    let sum = 0;
    for (let i = 0; i < a.length; i += 1) sum += a[i] * b[i];
    return sum;
  }

  // Reorthogonalized QR avoids squaring the condition number of the regression matrix.
  function leastSquares(rows, targets) {
    const size = rows[0].length;
    const columns = Array.from({ length: size }, (_, j) => rows.map((row) => row[j]));
    const scales = columns.map((column) => Math.sqrt(dot(column, column)));
    const q = [];
    const r = Array.from({ length: size }, () => Array(size).fill(0));
    for (let j = 0; j < size; j += 1) {
      if (!(scales[j] > 0)) return null;
      const v = columns[j].map((value) => value / scales[j]);
      for (let pass = 0; pass < 2; pass += 1) {
        for (let i = 0; i < j; i += 1) {
          const projection = dot(q[i], v);
          r[i][j] += projection;
          for (let k = 0; k < v.length; k += 1) v[k] -= projection * q[i][k];
        }
      }
      const length = Math.sqrt(dot(v, v));
      if (length < 1e-10) return null;
      r[j][j] = length;
      q.push(v.map((value) => value / length));
    }
    const coefficients = Array(size).fill(0);
    for (let j = size - 1; j >= 0; j -= 1) {
      let value = dot(q[j], targets);
      for (let k = j + 1; k < size; k += 1) value -= r[j][k] * coefficients[k];
      coefficients[j] = value / r[j][j];
    }
    return coefficients.map((value, j) => value / scales[j]);
  }

  function responseEstimate(training, x, di, bandwidthX, bandwidthDi, omit = -1) {
    let minDistance = Infinity;
    const distances = training.map((point, index) => {
      if (index === omit) return Infinity;
      const dx = (x - point.x) / bandwidthX;
      const dd = (di - point.di) / bandwidthDi;
      const distance = dx * dx + dd * dd;
      minDistance = Math.min(minDistance, distance);
      return distance;
    });
    if (!Number.isFinite(minDistance)) return NaN;
    let numerator = 0;
    let denominator = 0;
    distances.forEach((distance, index) => {
      if (!Number.isFinite(distance)) return;
      const weight = Math.exp(-0.5 * (distance - minDistance));
      numerator += weight * training[index].apd;
      denominator += weight;
    });
    return numerator / denominator;
  }

  function chooseBandwidth(training) {
    const xChoices = [0.0002, 0.0003, 0.0005, 0.00075, 0.001, 0.0015, 0.002, 0.003, 0.005];
    const diChoices = [10, 20, 30, 50, 80, 120, 200, 500];
    let best = { error: Infinity, x: null, di: null };
    for (const bandwidthX of xChoices) {
      for (const bandwidthDi of diChoices) {
        let error = 0;
        training.forEach((point, index) => {
          const predicted = responseEstimate(training, point.x, point.di, bandwidthX, bandwidthDi, index);
          error += (predicted - point.apd) ** 2;
        });
        if (error < best.error) best = { error, x: bandwidthX, di: bandwidthDi };
      }
    }
    return best;
  }

  function calibratedMap(result) {
    const { apds, xsInit, xsRepolarized, apdStartTimes, apdEndTimes } = result;
    const n = apds.length;
    if (n < 40 || !xsRepolarized || xsRepolarized.length !== n ||
        apdStartTimes.length !== n || apdEndTimes.length !== n) return null;
    const trainingBeats = Math.floor(0.6 * n);
    const rows = [];
    const targets = [];
    for (let i = 0; i < trainingBeats; i += 1) {
      rows.push(plateauFeatures(xsInit[i], apds[i]));
      targets.push(xsRepolarized[i]);
    }
    const coefficients = leastSquares(rows, targets);
    if (!coefficients) return null;

    let diSquared = 0;
    let diLogRatio = 0;
    const training = [];
    for (let i = 0; i < trainingBeats - 1; i += 1) {
      const di = apdStartTimes[i + 1] - apdEndTimes[i];
      if (!(di > 0 && xsRepolarized[i] > 0 && xsInit[i + 1] > 0)) return null;
      diSquared += di * di;
      diLogRatio += di * Math.log(xsInit[i + 1] / xsRepolarized[i]);
      training.push({ x: xsInit[i + 1], di, apd: apds[i + 1] });
    }
    const tauRecovery = -diSquared / diLogRatio;
    if (!(tauRecovery > 0 && Number.isFinite(tauRecovery))) return null;
    const bandwidth = chooseBandwidth(training);
    const predictions = [];
    let errorSquared = 0;
    let publishedErrorSquared = 0;
    let plateauErrorSquared = 0;
    for (let i = trainingBeats; i < n - 1; i += 1) {
      const di = apdStartTimes[i + 1] - apdEndTimes[i];
      if (!(di > 0)) return null;
      const xAtRepolarization = dot(coefficients, plateauFeatures(xsInit[i], apds[i]));
      const predictedX = xAtRepolarization * Math.exp(-di / tauRecovery);
      const predictedApd = responseEstimate(training, predictedX, di, bandwidth.x, bandwidth.di);
      const publishedApd = publishedReturnMap(apds[i], result.inputs.pcl);
      if (![predictedX, predictedApd, publishedApd].every(Number.isFinite)) return null;
      const actualApd = apds[i + 1];
      errorSquared += (predictedApd - actualApd) ** 2;
      publishedErrorSquared += (publishedApd - actualApd) ** 2;
      plateauErrorSquared += (xAtRepolarization - xsRepolarized[i]) ** 2;
      predictions.push({ beat: i + 2, currentApd: apds[i], actualApd, predictedApd, di });
    }
    if (!predictions.length) return null;
    return {
      predictions,
      trainingBeats,
      tauRecovery,
      bandwidthX: bandwidth.x,
      bandwidthDi: bandwidth.di,
      rmse: Math.sqrt(errorSquared / predictions.length),
      publishedRmse: Math.sqrt(publishedErrorSquared / predictions.length),
      plateauXrmse: Math.sqrt(plateauErrorSquared / predictions.length),
    };
  }

  window.XMemoryMap = { calibratedMap, publishedReturnMap };
})();
