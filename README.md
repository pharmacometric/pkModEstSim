# 📘 pkModEstSim JavaScript Model Estimation and Simulation Library 

`pkModEstSim` is a comprehensive JavaScript library for simulating and fitting pharmacokinetic (PK) models using differential equations. It enables researchers, students, and analysts to perform modeling directly in the browser without requiring backend computation or installation of complex software.

---

## 🚀 Introduction

Pharmacokinetic models describe how a drug is absorbed, distributed, metabolized, and eliminated by the body. Common models include one- and multi-compartment systems, linear and non-linear elimination, and various routes of administration. This library simulates these systems and fits observed data using multiple estimation strategies.

---

## 📦 Features

- Supports 1-, 2-, and 3-compartment models
- Absorption methods: Oral (first-order), IV bolus, IV infusion
- Elimination: Linear, Michaelis-Menten, or combined
- Model fitting via:
  - Levenberg–Marquardt
  - Gradient Descent
  - Simulated Annealing
- Confidence interval calculation
- Covariate testing with linear, exponential, or power functions
- Simulation of subject profiles with reproducible random seeds
- Output from all compartments and amount eliminated

---

## 📂 Dataset Preparation

### 🧬 How to Create a NONMEM-style Dataset in R

```r
set.seed(123)
n <- 10
subjects <- data.frame(
  ID = 1:n,
  AGE = sample(25:60, n, replace = TRUE),
  WT = round(runif(n, 60, 90), 1),
  SEX = sample(c("M", "F"), n, replace = TRUE),
  RACE = sample(c("White", "Black", "Asian", "Hispanic"), n, replace = TRUE),
  CRCL = round(runif(n, 75, 110), 1)
)

time_points <- c(0, 1, 2, 4, 6, 8, 12)
pkdata <- do.call(rbind, lapply(1:n, function(i) {
  subject <- subjects[i, ]
  data.frame(
    ID = subject$ID,
    TIME = time_points,
    DV = round(10 / (1 + time_points) + rnorm(length(time_points), 0, 0.5), 2),
    AMT = ifelse(time_points == 0, 100, 0),
    EVID = ifelse(time_points == 0, 1, 0),
    AGE = subject$AGE,
    WT = subject$WT,
    SEX = subject$SEX,
    RACE = subject$RACE,
    CRCL = subject$CRCL
  )
}))
write.csv(pkdata, "nonmem_example_data.csv", row.names = FALSE)
```

---

## 💡 Using Dataset in JavaScript

```js
// Example NONMEM-style dataset (10 subjects, 7 time points)
      const dataset = [];
      const subjects = 10;
      const times = [0, 1, 2, 4, 6, 8, 12];
      for (let id = 1; id <= subjects; id++) {
        for (let t of times) {
          dataset.push({
            ID: id,
            TIME: t,
            DV: +(50 * Math.exp(-0.2 * t) + (Math.random() - 0.5)).toFixed(2),
            AMT: t === 0 ? 100 : 0,
            EVID: t === 0 ? 1 : 0
          });
        }
      }

      // Define model spec for 1-compartment with IV bolus and linear elimination
      const modelSpec = {
        compartments: 1,
        route: pkModEstSim.Routes.IV_BOLUS,
        elimination: pkModEstSim.Elimination.LINEAR
      };

      const initParams = {
        CL: 5,
        V1: 20
      };

      const result = pkModEstSim.fitModel(dataset, modelSpec, initParams);
```

---

## 🔬 Function Reference

### `fitModel(dataset, modelSpec, initParams, bounds, method, seed)`
Fits a PK model using the specified algorithm.

### `simulateSubjects(nSubjects, times, modelSpec, trueParams, seed)`
Simulates concentration-time data for a population.

### `testCovariateEffect(dataset, parameterName, subjectData, covariateKey, transform)`
Evaluates how covariates influence model parameters.

---

## 📊 Example Outputs

- Compartment concentration: A1/V1
- Peripheral compartments: A2, A3
- Amount eliminated: A_elim
- PK Metrics: Cmax, AUC, T1/2, trough

---

## 📘 Example Usage

### Model Fitting (selected examples)

```js
pkModEstSim.fitModel(exampleDataset, {
  compartments: 2,
  route: pkModEstSim.Routes.IV_BOLUS,
  elimination: pkModEstSim.Elimination.LINEAR
}, {
  CL: 1, V1: 15, Q2: 0.5, V2: 20
}, null, "levenberg-marquardt", 42);
```

### Simulation

```js
pkModEstSim.simulateSubjects(20, [0,1,2,4,6,8], {
  compartments: 1,
  route: pkModEstSim.Routes.ORAL,
  elimination: pkModEstSim.Elimination.LINEAR
}, {
  CL: 1, V1: 10, KA: 1
}, 123);
```

#### simulation example

```js
const modelSpec = {
        compartments: 1,
        route: pkModEstSim.Routes.IV_BOLUS,
        elimination: pkModEstSim.Elimination.LINEAR
      };

const simResult = pkModEstSim.simulateSubjects({
        n: 100,
        modelSpec,
        parameters: {
  CL: 1, V1: 10, KA: 1
},
        variability: {
          CL: 0.25,
          V1: 0.25
        },
        dosing: [{ time: 0, amt: 100 }],
        times
      });
	  
	  

      const simByTime = {};
      times.forEach(t => simByTime[t] = []);
      simResult.forEach(row => {
        simByTime[row.TIME].push(row.CONC);
      });

      const medianLine = [];
      const lower90Line = [];
      const upper90Line = [];

      for (let t of times) {
        const values = simByTime[t].sort((a, b) => a - b);
        const mid = Math.floor(values.length / 2);
        const median = values.length % 2 === 0
          ? (values[mid - 1] + values[mid]) / 2
          : values[mid];
        medianLine.push({ x: t, y: median });
        lower90Line.push({ x: t, y: values[Math.floor(0.05 * values.length)] });
        upper90Line.push({ x: t, y: values[Math.floor(0.95 * values.length)] });
      }

```

---

## 📚 Future Enhancements

- Random effect modeling
- PK-PD extensions
- Output to CSV, PDF
- Embedded plotting tools

MIT License © 2025 — pkModEstSim Contributors
