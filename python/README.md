# python

The models and the pipeline that turns their solutions into the paper's exhibits. `InformalSavings/`, `US/` and
`paper/` are folders of scripts. Each script puts its own folder on `sys.path` and imports its siblings by bare
name, so run them from the repository root as `python python/<folder>/<script>.py`. `gridsearch/` is a package,
installed by `pip install -e python` (`pyproject.toml`). Docstrings cite equations of the technical documentation
by label (e.g. `eq:governmentBudget`); look the label up in `docs/equationLabels.csv` for its number and page in
`docs/TechnicalDocumentation.pdf`.

## `paper/`: the pipeline (it owns no economics)

| Stage | Files |
|---|---|
| everything | `replicate.py`: every step, in order of dependence, at the paper's settings (`--dry`, `--list`, `--quick`, `--only`, `--from`) |
| (0) data | `dataTargets.py` (Argentina's capital-output target, Penn World Table via FRED), `oecdFigure1.py` (figure 1's cross-section), `oecd29Data.py` (the 29 countries' inputs), `oecd26Sample.py` (the sample rule), `oecd26PensionTiers.py` (replacement rates and spending by tier; Windows only, it calls PowerShell), `oecd26Topups.py` (top-up spending and shares), `oecd26Theta.py` (the measure of the design, `data/oecd26_theta.csv`). These reach the network and write `data/` only. |
| (i) calibration | `runCalibration.py` (Argentina), `runCalibrationUS.py` (the U.S., France and the UK, and the cost of redistribution) |
| (ii) experiments | `runShocks.py` (Argentina's reform and comparative statics), `runShocksUS.py` (the counterfactuals, the endogenous design, the cross-country test) |
| (iii) build | `build.py` (the registry of every output; `--list`, `--only`, `--no-copy`), `config.py` (paths, and every setting behind a published number), `datasets.py` (the only module that knows the layout of `results/`), the builders `tables.py`, `figures.py` (Argentina), `tablesUS.py`, `figuresUS.py`, `tablesOECD26.py`, `figuresOECD26.py`, `tablesHorizon.py`, `tablesOA.py`, `figuresOA.py` (the supplemental appendix's own exhibits), `siteData.py` (the datasets of its interactive figures) |
| tools | `compareResults.py`: every changed csv under `results/` against a git ref, rows aligned on the file's keys, per-column differences above `--tol` |

Stages (i) and (ii) launch the model scripts below as child processes, skip work whose output exists (`--force`
redoes it) and print what they would run with `--dry` (`--list` where they have it). Stage (iii) imports no model
code and unpickles nothing: it reads `results/` and writes `results/paper/` in seconds. Always pass `--no-copy` in
this repository. Without it, `build.py` also copies its output into `writing/Paper`, the paper's source, which is
not part of this repository.

## `InformalSavings/`: Argentina (section 5, appendices B and E)

The overlapping-generations economy with an informal type that saves, calibrated to Argentina's 2010 targets.

| | |
|---|---|
| `base.py`, `model.py`, `policy.py` | the economic equilibrium (`base`), the model class `ModelInformalSavings` with its steady state, initial state and calibration (`model`), the politico-economic equilibrium under log and CRRA preferences, a backward recursion over the state (`policy`) |
| `test.py` | the loader of `data/ArgentinaTest.xlsx` (not a test) |
| `calibrateRhoGrid.py` | the calibration over the grid of ρ, one pickled instance per point in `results/calibration/instances*/` (`--epsPre0`: the ε = 0 variant) |
| `shockUniversal.py`, `shockEEOnly.py` | the reform (the universal pension), with taxes re-optimised and with taxes held at the baseline path |
| `sweepEpsThetaGrid.py` | comparative statics over the (ε, θ) plane |
| `reformDiagnostics.py`, `stationaryApprox.py`, `runHorizon.py` | the reform's first order conditions and savings by group; stationary against date-specific policy functions; appendix D's horizon and residual checks |

## `US/`: the U.S., France, the UK and the OECD cross-section (sections 6 and 7, appendices C, F and G)

The same economy without the informal type, with the endogenous design of section 7 and the cross-country test.

| | |
|---|---|
| `base.py`, `model.py`, `policy.py` | the economic equilibrium, `ModelUS` with its steady state and calibration, the tax policy under log and CRRA preferences |
| `modelFR.py` | `ModelFR(ModelUS)`: the calibration protocol of France and the UK |
| `modelESC.py`, `policyESC.py` | the endogenous design θ (chosen one period ahead, or once and for all) and the cost of redistribution |
| `shocks.py` | the counterfactuals, one function per scenario |
| `oecd26.py` | `ModelOECD`: each OECD country from `data/oecd29*.csv` and `data/oecd26_theta.csv` |
| `test.py`, `testEU.py` | the loaders of `data/USMain_test.xlsx` and of `data/FRMain.xlsx`, `data/UKMain.xlsx` (not tests) |
| drivers | `calibrateRhoGrid.py`, `calibrateRhoGridEU.py` (calibrations), `runShocksUS.py`, `runShocksXi.py` (counterfactuals), `runESC.py`, `runESCcrra.py`, `runESCxi.py`, `collectESCexperiments.py` (the endogenous design), `runTauTheta.py` (the 2020 tax at a pinned design), `runOECD26.py` (the cross-country test), `stationaryApprox.py`, `stationaryApproxESC.py` (stationary policies), `runHorizon.py` (appendix D), `c5Numbers.py` (the smaller numbers the text quotes) |

The U.S. commands carry `--theta C --main`: the paper's measure of the design, `config.US['theta']`.

## `gridsearch/`

Root finding and selection on grids, shared by both models: `robustRoot.py` (a bounded root problem as an
unconstrained one), `roots1d.py` (every root or maximum along a grid, and the selection among candidate
equilibria), `cartesian.py` (product grids), `interp.py` (gridded interpolants and smoothing), `continuation.py`
(a parameter march with extrapolated starting values), `testing.py` (the PASS/FAIL harness some drivers report
with).
