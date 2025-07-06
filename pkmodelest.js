/**
 * PKModelFit.js
 * A comprehensive JavaScript library for pharmacokinetic model fitting using ODEs.
 * Supports NONMEM-style datasets, multiple compartment models, various elimination types,
 * and routes of administration. Also supports simulation of multiple subjects and covariate testing.
 * Dependencies: numeric.js, ml-levenberg-marquardt (or a suitable optimizer)
 */

const pkModEstSim = (() => {

  // ========== UTILS ========== //
  function mean(arr) {
    return arr.reduce((a, b) => a + b, 0) / arr.length;
  }

  function stddev(arr) {
    const avg = mean(arr);
    const sqDiffs = arr.map(v => Math.pow(v - avg, 2));
    return Math.sqrt(mean(sqDiffs));
  }

  function seedRandom(seed) {
    let x = Math.sin(seed++) * 10000;
    return () => {
      x = Math.sin(seed++) * 10000;
      return x - Math.floor(x);
    };
  }

  function clone(obj) {
    return JSON.parse(JSON.stringify(obj));
  }

  function linearRegression(x, y) {
    const xMean = mean(x);
    const yMean = mean(y);
    const numerator = x.reduce((sum, xi, i) => sum + (xi - xMean) * (y[i] - yMean), 0);
    const denominator = x.reduce((sum, xi) => sum + Math.pow(xi - xMean, 2), 0);
    const slope = numerator / denominator;
    const intercept = yMean - slope * xMean;
    const residuals = y.map((yi, i) => yi - (slope * x[i] + intercept));
    const sse = residuals.reduce((sum, r) => sum + r * r, 0);
    const se = Math.sqrt(sse / (x.length - 2));
    const tStat = slope / (se / Math.sqrt(denominator));
    const pValue = 2 * (1 - normalCDF(Math.abs(tStat)));
    return { slope, intercept, tStat, pValue };
  }

  function normalCDF(z) {
    return 0.5 * (1 + Math.erf(z / Math.sqrt(2)));
  }

  // ========== ROUTES ========== //
  const Routes = {
    IV_BOLUS: 'iv_bolus',
    IV_INFUSION: 'iv_infusion',
    ORAL: 'oral'
  };

  // ========== ELIMINATION TYPES ========== //
  const Elimination = {
    LINEAR: 'linear',
    MM: 'michaelis_menten',
    COMBINED: 'combined'
  };

  // ========== COVARIATE TESTING FUNCTION ========== //
  function testCovariateEffect(dataset, parameterName, subjectData, covariateKey, transform = 'linear') {
    let values = [];
    let covariates = [];
    for (let id in subjectData) {
      const param = subjectData[id][parameterName];
      const cov = dataset.find(r => r.ID === id)[covariateKey];
      if (param !== undefined && cov !== undefined) {
        let transformedCov = cov;
        if (transform === 'log') transformedCov = Math.log(cov);
        if (transform === 'exp') transformedCov = Math.exp(cov);
        if (transform === 'power') transformedCov = Math.pow(cov, 0.75);
        values.push(param);
        covariates.push(transformedCov);
      }
    }
    return linearRegression(covariates, values);
  }

  // ========== MODEL ODE DEFINITIONS ========== //
  function odeSystem(t, y, p, model) {
    const A = y;
    const { CL, V1, V2, V3, Q2, Q3, KA, VMAX, KM } = p;

    let C1 = A[0] / V1;
    let dA = new Array(y.length).fill(0);

    if (model.route === Routes.ORAL) {
      dA[0] += KA * A[1];
      dA[1] = -KA * A[1];
    } else if (model.route === Routes.IV_INFUSION) {
      dA[0] += model.infusionRate;
    }

    if (model.elimination === Elimination.LINEAR) {
      dA[0] += -CL * C1;
    } else if (model.elimination === Elimination.MM) {
      dA[0] += -VMAX * C1 / (KM + C1);
    } else if (model.elimination === Elimination.COMBINED) {
      dA[0] += -CL * C1 - VMAX * C1 / (KM + C1);
    }

    if (model.compartments >= 2) {
      let C2 = A[2] / V2;
      dA[0] += -Q2 * (C1 - C2);
      dA[2] += Q2 * (C1 - C2);
    }
    if (model.compartments === 3) {
      let C3 = A[3] / V3;
      dA[0] += -Q3 * (C1 - C3);
      dA[3] += Q3 * (C1 - C3);
    }

    dA[y.length - 1] = -dA[0];
    return dA;
  }

  // ========== NUMERICAL INTEGRATION ========== //
  function simulateODE(times, y0, p, model) {
    const dt = 0.1;
    let result = [];
    let y = [...y0];
    let t = 0;
    for (let i = 0; i < times.length; i++) {
      while (t < times[i]) {
        const dy = odeSystem(t, y, p, model);
        for (let j = 0; j < y.length; j++) {
          y[j] += dy[j] * dt;
        }
        t += dt;
      }
      result.push([t, ...y]);
    }
    return result;
  }

  // ========== FITTING FUNCTION ========== //
  function fitModel(dataset, modelSpec, initParams, bounds, method = 'levenberg-marquardt', seed = 42) {
    let grouped = {};
    for (let row of dataset) {
      if (!grouped[row.ID]) grouped[row.ID] = [];
      grouped[row.ID].push(row);
    }

    function residuals(paramsArr) {
      let params = Object.fromEntries(Object.keys(initParams).map((k, i) => [k, paramsArr[i]]));
      let error = [];
      for (let id in grouped) {
        const patient = grouped[id];
        const times = patient.map(r => r.TIME);
        const y0 = [patient[0].AMT || 0, 0, 0, 0, 0];
        const sim = simulateODE(times, y0, params, modelSpec);
        for (let i = 0; i < patient.length; i++) {
          if (patient[i].EVID === 0) {
            let Cp = sim[i][1] / params.V1;
            error.push(patient[i].DV - Cp);
          }
        }
      }
      return error;
    }

    let optimizer;
    const paramKeys = Object.keys(initParams);
    let current = Object.values(initParams);

    if (method === 'gradient-descent') {
      let lr = 0.001;
      for (let epoch = 0; epoch < 200; epoch++) {
        let grad = new Array(current.length).fill(0);
        let baseLoss = residuals(current).reduce((s, r) => s + r * r, 0);
        for (let i = 0; i < current.length; i++) {
          let tmp = [...current];
          tmp[i] += 1e-4;
          let loss = residuals(tmp).reduce((s, r) => s + r * r, 0);
          grad[i] = (loss - baseLoss) / 1e-4;
        }
        for (let i = 0; i < current.length; i++) {
          current[i] -= lr * grad[i];
        }
      }
      optimizer = { solution: current };

    } else if (method === 'simulated-annealing') {
      let rand = seedRandom(seed);
      let T = 1.0, coolingRate = 0.95;
      let best = [...current];
      let bestLoss = residuals(best).reduce((s, r) => s + r * r, 0);
      for (let i = 0; i < 500; i++) {
        let candidate = best.map(v => v + (rand() - 0.5) * 0.1);
        let loss = residuals(candidate).reduce((s, r) => s + r * r, 0);
        if (loss < bestLoss || Math.exp((bestLoss - loss) / T) > rand()) {
          best = candidate;
          bestLoss = loss;
        }
        T *= coolingRate;
      }
      optimizer = { solution: best };

    } else if (method === 'levenberg-marquardt') {
      optimizer = numeric.uncmin(
        (p) => residuals(p).reduce((s, r) => s + r * r, 0),
        current
      );
    } else if (method === 'SAEM' || method === 'FOCEI') {
      throw new Error(`${method} is a population method requiring random effects modeling and is not implemented yet.`);
    }

    const estimates = Object.fromEntries(paramKeys.map((k, i) => [k, optimizer.solution[i]]));
    const residualVec = residuals(optimizer.solution);
    const mse = residualVec.reduce((s, v) => s + v * v, 0) / residualVec.length;
    const se = Math.sqrt(mse);
    const ci = Object.fromEntries(paramKeys.map(k => [k, [estimates[k] - 1.96 * se, estimates[k] + 1.96 * se]]));

    return { estimates, ci, residuals: residualVec, mse };
  }

  // ========== SIMULATION FUNCTION ========== //
  function simulateSubjects(nSubjects, times, modelSpec, trueParams, seed = 123) {
    const rand = seedRandom(seed);
    let results = [];
    for (let i = 0; i < nSubjects; i++) {
      let y0 = [trueParams.AMT || 0, 0, 0, 0, 0];
      let sim = simulateODE(times, y0, trueParams, modelSpec);
      results.push(sim.map(row => ({
        subject: i + 1,
        time: row[0],
        A1: row[1],
        A2: row[2],
        A3: row[3],
        A_abs: row[1],
        A_peri: row[2],
        A_elim: row[4],
        Cp: row[1] / trueParams.V1
      })));
    }
    return results;
  }

  return {
    fitModel,
    simulateSubjects,
    testCovariateEffect,
    Routes,
    Elimination
  };
})();