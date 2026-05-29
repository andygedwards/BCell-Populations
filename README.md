# Heterogeneous Beta Cell Population Simulation Pipeline

**Author:** Andy Edwards  
**Model Base:** Modified version of Riz et al. (2014) pancreatic beta cell electrophysiology model with low voltage INa variant

This repository provides a three-script pipeline for constructing a heterogeneous population of pancreatic beta cells, simulating their glucose responses, and post-processing the results to classify electrophysiological phenotypes.

---

## Repository Structure

```
├── Baseline models/          # ODE model definition files (parameters, states, RHS)
│   ├── Riz_2014_original/
│   └── Riz_2014_INa_low/
├── Create Population/        # Output of Step 1 (populated on execution)
├── Glucose screen/           # Output of Steps 2 & 3 (populated on execution)
│   ├── data/                 # Raw simulation output (Step 2)
│   └── results/              # Post-processed results (Step 3)
├── Optimization/             # Pre-computed optimization results (.mat files)
│   ├── Cost/                 # Cost function definitions
│   ├── GA/                   # Genetic algorithm runners
│   ├── DE/                   # Differential evolution runners & results
│   └── PSO/                  # Particle swarm optimization runner
├── Shared/                   # Utility functions used across scripts
├── create_het_pop.m          # Step 1: Build heterogeneous cell population
├── sim_gluc_responses.m      # Step 2: Simulate glucose responses
└── run_population_postprocess.m  # Step 3: Classify and post-process results
```

> **Note:** The `Create Population/`, `Glucose screen/data/`, and `Glucose screen/results/` directories are empty at distribution and will be populated when the scripts are executed in order.

---

## Prerequisites

- MATLAB with the Parallel Computing Toolbox (recommended for `parfor` loops)
- All three scripts must be run from the **repository root directory**

---

## Step 1 — `create_het_pop.m`

Constructs a heterogeneous virtual cell population by sampling model parameters from optimized distributions, then computes voltage-clamp summary metrics for each cell.

### User-Configurable Inputs

| Variable | Default | Description |
|---|---|---|
| `nTrials` | `3000` | Number of cells in the population |
| `base_paramFile` | `@Riz2014_init_parameters_INa_low` | Baseline parameter function handle |
| `Experiments` | `{'Combined'}` | Optimization type(s) to use. Options: `'Combined'` (published), `'Na'`, `'Tail'` |
| `frac_high` | `[0.370]` | Fraction of cells assigned to the high-INa subpopulation. The published optimal value is `0.370`. Leave empty to use values from the optimization result. |

### Required Input Files (pre-existing in `Optimization/`)

- `Combined_refined.mat` or `Na_refined.mat` — refined optimization results containing fitted parameter structs with fields including `Total_Error`, `V_hNa`, `g_Na`, `V_mNa`, `V_hNa_low`, `g_Na_low`, `V_mNa_low`, `g_CaL`, `g_CaPQ`, `g_CaT`, `Frachigh`

### Outputs

Saved to `Create Population/<Experiment>/<nTrials>/<frac_high>/` (subfolders `A`/`B` created automatically when duplicate `frac_high` values exist):

| File | Contents |
|---|---|
| `modParam_scaling_<nTrials>_<frac_high>.mat` | Parameter scaling factor matrix (`nTrials × nParams`) |
| `modParam_stdev_<nTrials>_<frac_high>.mat` | Standard deviations used for each parameter |
| `modParam_names_<nTrials>_<frac_high>.mat` | Parameter name list |
| `peak_INa_<nTrials>_<frac_high>.mat` | Peak sodium current per cell |
| `v_half_act_<nTrials>_<frac_high>.mat` | Voltage of half-activation per cell |
| `v_half_<nTrials>_<frac_high>.mat` | Voltage of half-inactivation per cell |
| `peak_ICa_<nTrials>_<frac_high>.mat` | Peak calcium current per cell |
| `late_ICa_<nTrials>_<frac_high>.mat` | Late calcium current per cell |
| `EE_<nTrials>_<frac_high>.mat` | Exocytosis efficiency per cell |
| `TE_<nTrials>_<frac_high>.mat` | Total exocytosis per cell |
| `frac_high_<nTrials>_<frac_high>.mat` | `frac_high` value used |
| `modParam_scaling_*.png` (figures) | Distribution plots for each varied parameter |

### Variable Parameters Sampled

Non-INa channel parameters: `V_GK_max`, `K_GK`, `g_Kv`, `g_BK`, `g_CaL`, `g_CaPQ`, `g_CaT`, `g_KATP_hat`  
High-INa subpopulation: `V_hNa`, `n_hNa`, `g_Na`, `V_mNa`  
Low-INa subpopulation: `V_hNa_low`, `n_hNa_low`, `g_Na_low`, `V_mNa_low`

---

## Step 2 — `sim_gluc_responses.m`

Simulates the electrophysiological response of each cell in the population to low (2 mM) and high (20 mM) glucose conditions using the stiff ODE solver `ode15s`.

### User-Configurable Inputs

| Variable | Default | Description |
|---|---|---|
| `Optimizations` | `{'Combined'}` | Must match the experiment used in Step 1 |
| `nTrials` | `3000` | Must match the population size created in Step 1 |
| `frac_high` | `[0.37]` | Must match the population created in Step 1 |
| `ka` | `4e-05` | Electro-metabolic coupling parameter (`k_A`). Published value: `4e-05` |
| `ka_default` | `1e-04` | Default model value of `k_A`; used to scale `ka` relative to baseline |

### Required Input Files

- Population files from Step 1 (`modParam_scaling_*.mat`, `modParam_names_*.mat`) located in `Create Population/`

### Simulation Protocol (per cell)

1. **Initialization:** 60,000 ms run at 2 mM glucose to reach steady state
2. **Low glucose:** 600,000 ms ODE integration from the initialized state at 2 mM glucose
3. **High glucose:** 600,000 ms ODE integration from the same initialized state at 20 mM glucose

### Outputs

Saved to `Glucose screen/data/<Experiment>/<nTrials>/<frac_high>/ka_<ka>/`:

| Location | File | Contents |
|---|---|---|
| `cells/` | `<i>_low.mat` | Time vector `T`, state matrix `Y`, and parameter vector `param` for cell `i` at low glucose |
| `cells/` | `<i>_high.mat` | Same as above for high glucose condition |
| `metadata/` | `metadata.mat` | `modParam_names`, `nModParams`, `nTrials`, `model`, `modParam_inds` |

State matrix `Y` columns follow the Riz2014 state variable order; membrane voltage is column 15 and the metabolic variable `a` is column 12.

---

## Step 3 — `run_population_postprocess.m`

Loads the simulated glucose responses, detects action potential peaks and troughs, classifies each cell's electrophysiological phenotype, and saves summary results.

### User-Configurable Inputs

| Variable | Default | Description |
|---|---|---|
| `nTrials` | `3000` | Must match Steps 1 & 2 |
| `init_frac_high` | `[0.37]` | Population(s) to analyze. Leave empty (`[]`) to process all available subfolders |
| `ka` | `0.00004` | Selects which `ka` subfolder to load from `Glucose screen/data/` |
| `Experiment` | `'Combined'` | Must match Steps 1 & 2 |
| `plot_peaks_troughs` | `0` | Set to `1` to overlay detected peaks/troughs on time-series plots |
| `plotting.number` | `[]` | Vector of specific cell indices to plot (range `1:nTrials`) |
| `plotting.pheno` | `{}` | Plot cells matching given phenotype(s): `'Silent'`, `'Bursting'`, `'Spiking'`, `'Depolarized'`, `'Other'` |
| `plotting.transitions` | `{}` | Plot cells matching low→high glucose transitions: `'SiSi'`, `'SiB'`, `'SiSp'`, `'SiD'`, `'SiO'` |
| `get_user_input` | `1` | If `1`, presents a dialog to manually correct automated phenotype classifications |

### Phenotype Classification

Each cell is classified independently under low and high glucose using the `Andrean` function based on detected peaks and troughs in membrane voltage (threshold: −40 mV). Possible phenotypes:

- `Silent` — no action potentials
- `Bursting` — periodic bursts of action potentials
- `Spiking` — continuous spiking
- `Depolarized` — sustained depolarization
- `Other` — does not meet criteria for the above

Metabolic state (active/inactive) is assessed separately based on whether the metabolic variable `a` exceeds 2.0 mM during the second half of the simulation.

### Outputs

Saved to `Glucose screen/results/<Experiment>/<nTrials>/<frac_high>/ka_<ka>/`:

| File | Contents |
|---|---|
| `low_glucose.mat` | `class_low`, `outs_low`, `param_low`, `metab_low` — phenotype, output metrics, parameters, and metabolic state for each cell at low glucose |
| `high_glucose.mat` | `class_high`, `outs_high`, `param_high`, `metab_high` — same for high glucose |
| `metadata.mat` | `modParam_names`, `nModParams`, `nTrials`, `model`, `modParam_inds` |
| `Summary/` | Directory reserved for summary figures (populated by downstream analysis) |

---

## Running the Pipeline

```matlab
% From the repository root in MATLAB:

% Step 1 — build the population (~minutes with parfor)
run('create_het_pop.m')

% Step 2 — simulate glucose responses (~hours with parfor)
run('sim_gluc_responses.m')

% Step 3 — post-process and classify
run('run_population_postprocess.m')
```

Ensure `nTrials`, `frac_high`, `ka`, and `Experiment` are set consistently across all three scripts before running.

---

## Citation / Reference

This pipeline implements the beta cell population model described in:

> Riz, M., Pedersen, M.G. et al. (2014). *Missing inactivation of Na⁺ current...* (see manuscript for full citation)

Please cite the associated publication if you use this code.
