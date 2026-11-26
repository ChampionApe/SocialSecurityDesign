# The Political Economy of Pension Design: replication package

Rasmus K. S. Berg (University of Copenhagen) and Martín Gonzalez-Eiras (University of Bologna)

[![replicate](https://github.com/ChampionApe/SocialSecurityDesign/actions/workflows/replicate.yml/badge.svg)](https://github.com/ChampionApe/SocialSecurityDesign/actions/workflows/replicate.yml)

**Supplemental appendix:** https://championape.github.io/SocialSecurityDesign ·
**The paper's tables and figures, as built here:** [`docs/PaperExhibits.pdf`](docs/PaperExhibits.pdf) ·
**The supplemental appendix's:** [`docs/SupplementalExhibits.pdf`](docs/SupplementalExhibits.pdf) ·
**Technical documentation** (models, equilibria, algorithms, data): [`docs/TechnicalDocumentation.pdf`](docs/TechnicalDocumentation.pdf)

## Overview

The code solves an overlapping-generations model with heterogeneous households, in which a government without
commitment chooses the pension tax each period by probabilistic voting and, in section 7, the electorate also
chooses the pension design. It produces every table and figure of the paper and of its supplemental appendix.
There are two arms:

- **Argentina** (section 5, appendix E), `python/InformalSavings/`: the 2005–10 extension of pensions to the
  informal sector.
- **The U.S., France, the UK and 21 OECD countries** (sections 6 and 7, appendices C, F and G), `python/US/`.

`python/paper/` runs them in stages: (0) data, (i) calibration, (ii) experiments, (iii) build. Stage (iii) reads
only `results/` and writes every exhibit in seconds. The committed `results/` hold the output of stages (i) and (ii),
so the exhibits can be rebuilt without solving anything. Re-solving everything takes about a day on a laptop. One
command runs all of it: `python python/paper/replicate.py`. Working with an AI assistant such as Claude Code:
`CLAUDE.md` orients it.

## Data availability and provenance

### Statement about rights

- [x] The authors have legitimate access to and permission to use the data used in this manuscript.
- [x] The authors have permission to redistribute the data contained in this replication package. Every file in
  `data/` holds aggregated statistics (group means and shares, country-level indicators) computed from the public
  sources below. None holds microdata. Each provider's terms (table below) allow the statistics to be published
  with attribution. The WID and WIID values are licensed for non-commercial use only and keep their providers'
  share-alike licences (`LICENSE-data`).

### Summary of availability

- [x] All data **are** publicly available. The household surveys behind the group statistics (LIS, the BHPS, the
  CPS, the EPH) are available to registered users from their providers. This package contains only the statistics
  tabulated from them.

### Sources

`data/README.md` describes every file: content, source, vintage, the code that reads it and where the paper uses
it. The `*_sources.csv` files and `oecd26_pensionTiers_harmonised.csv` give each value's URL and retrieval date.
The authors' own compilations are licensed CC BY 4.0 (`LICENSE-data`). The third-party data in them remain under
their providers' terms. The Eurostat, OECD and other statistics are modified here: they are selected, converted,
combined and, where the files say so, fitted. Their providers are not responsible for these uses.

| Source | Files | Terms | Cite as |
|---|---|---|---|
| OECD: Social Expenditure Database (SOCX), *Pensions at a Glance* 2013/2021/2025, earnings database, hours worked, national accounts (table 29), Global Pension Statistics, Income Distribution Database, *Education at a Glance*, historical population | workbooks, `oecdFigure1.csv`, `oecd29*.csv`, `oecd26_*.csv` | content published from 1 July 2024: CC BY 4.0; earlier content: the OECD Terms and Conditions (https://www.oecd.org/en/about/oecd-open-by-default-policy.html) | OECD, the dataset and its year |
| UN, *World Population Prospects* 2019, 2022, 2024 (medium variant) | workbooks, `oecd29.csv`, `oecd29_population.csv` | CC BY 3.0 IGO (https://creativecommons.org/licenses/by/3.0/igo/) | United Nations, DESA, Population Division (2024), World Population Prospects 2024 |
| Penn World Table 11.0 (through FRED) | `argentina_*.csv`, `oecd29.csv` | CC BY 4.0 (https://www.rug.nl/ggdc/productivity/pwt/) | Feenstra, Inklaar and Timmer (2015), *American Economic Review* 105(10) |
| World Bank, World Development Indicators | `oecdFigure1.csv` (`pop*_wb`) | CC BY 4.0, the World Bank Terms of Use for Datasets (https://www.worldbank.org/en/about/legal/terms-of-use-for-datasets) | World Bank, WDI |
| UNU-WIDER, World Income Inequality Database, version of 8 September 2026 | `oecdFigure1.csv` (`gini`, `gini_wiidMarket` and their year, id and source columns) | free download (https://www.wider.unu.edu/database/world-income-inequality-database-wiid). The database names no licence of its own, so these columns are taken to carry UNU-WIDER's default, CC BY-NC-SA 3.0 IGO (https://creativecommons.org/licenses/by-nc-sa/3.0/igo/): non-commercial use, share-alike | UNU-WIDER, World Income Inequality Database (WIID), version 8 September 2026 |
| World Inequality Database | `oecdFigure1.csv` (`gini_widPretax`, `gini_widPretax_year`), `oecd29_distribution.csv` (sources `wid_ptinc`, `wid_diinc`, validation only) | CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/): non-commercial use, share-alike | WID.world |
| Comparative Study of Electoral Systems, modules 5 and 6, Integrated Module Dataset | workbooks, `oecd29.csv` | free public datasets (registration at cses.org is voluntary). The package holds turnout rates computed by the authors, not survey records | CSES, the module and its release, with the release's DOI |
| ILOSTAT (SDG 8.3.1, 1.3.1); ILO modelled estimates | `oecd29.csv`, `oecd26_sample.csv`, `oecd29_distribution.csv` | CC BY 4.0 for ILO databases from 3 May 2023 (https://www.ilo.org/node/9321) | ILO, ILOSTAT |
| Eurostat (ESSPROS, EU-SILC, `nama_10_gdp`) | `oecd26_pensionTiers_harmonised.csv`, `oecd26_topups.csv`, `oecd29_distribution.csv` | reuse authorised with acknowledgement of the source; modifications stated (above) (https://ec.europa.eu/eurostat/help/copyright-notice) | Eurostat, the dataset code |
| Luxembourg Income Study, 2018 wave (France) | `FRMain.xlsx` | the microdata are reached only through LIS's remote system, which returns aggregates; the package holds three group means and shares tabulated there, published under the LIS citation rules | LIS Cross-National Data Center in Luxembourg |
| Harmonised British Household Panel Survey, wave 17 (UK Data Service) | `UKMain.xlsx` | End User Licence (UK Data Service CD137): outputs may be published if they meet its disclosure-control standards (clause 8), with the acknowledgement and citation of clause 11. The package holds three group means and shares, not data | University of Essex, ISER, the BHPS release |
| CPS March Supplement 2019; U.S. Census, *Voting and Registration in the Election of November 2020* | `USMain_test.xlsx` | works of the U.S. Government, in the public domain (17 U.S.C. § 105) | U.S. Census Bureau |
| INDEC, Encuesta Permanente de Hogares, second half of 2006 | `ArgentinaTest.xlsx` | public microdata; INDEC asks to be cited as the primary source (https://www.indec.gob.ar/Institucional/Indec/BasesDeDatos). The package holds quartile statistics | Fuente: Instituto Nacional de Estadística y Censos (INDEC) |
| National administrative sources for the minimum-pension top-ups (INPS via Itinerari Previdenziali, DREES, the Austrian federal budget, Bundestag and Destatis, the Belgian Court of Audit, CFP/IGFSS, NAV, ZUS) | `oecd26_topups.csv` | published statistics, cited row by row in the file | each row's source |
| Cetrángolo and Grushka (2020); Rofman and Oliveri (2012); MedPAC, *Data Book* 2020 | `ArgentinaTest.xlsx`, `USMain_test.xlsx` | published figures | as in the paper's references |

### Not included

- The raw files that stage (0) downloads (about 730 MB; `logs/`, ignored by git). The scripts download them again
  (`--refetch` where offered), and every derived value records its URL and retrieval date.
- The household microdata and the authors' tabulations of them: the group statistics of `FRMain.xlsx` (LIS) and
  `UKMain.xlsx` (BHPS), and the per-cell sources of `ArgentinaTest.xlsx` (EPH). See *Known gaps*.

## Computational requirements

### Software

- Python 3.14 (3.14.6 used), with the pinned packages of `requirements.txt` (numpy, pandas, scipy, matplotlib,
  openpyxl, and `symMaps`/`pyDbs` for the model's symbol tables) and this repository's `gridsearch` package
  (`pip install -e python`).
- Only for stage (0)'s `python/paper/oecd26PensionTiers.py`: Windows PowerShell.
- Only to typeset `docs/*.pdf`, which are built by the authors and committed: TeX Live 2026.

### Randomness

None: every solver is deterministic. Results are bitwise reproducible within a process but not across processes
or platforms, where the order of floating-point operations can differ. A re-solve therefore matches the committed
csvs to round-off, and every printed digit of the exhibits.

### Memory, run time, storage

The results were computed on a laptop with an Intel Core i7-1185G7 (4 cores, 8 threads), 32 GB of RAM and
Windows 11. The repository takes about 30 MB, and stage (0)'s raw downloads another 730 MB. Run times, cold, of
the steps of `python/paper/replicate.py`:

| Step | What | Time |
|---|---|---|
| `build` | stage (iii): every table and figure from `results/` | seconds |
| `calibARG` | stage (i), Argentina, the 16 points of the ρ grid | ~3 h |
| `calibUS` | stage (i), the U.S., France and the UK at ρ = 0.5, 1, 2, with the cost of redistribution | ~40 min |
| `shocksARG` | stage (ii), Argentina's reform and the (ε, θ) plane | ~1 h |
| `shocksUS` | stage (ii), the counterfactuals, the endogenous design and the cross-country test | ~3 h, plus 1.5–6 h for the supplemental appendix's data variants |
| `argEps0`, `argDiag` | the ε = 0 Argentine calibration and its reform; the reform's diagnostics | not timed (the sweep and reform of `calibARG` and `shocksARG` at ε = 0) |
| `c5` | the smaller U.S. numbers the text quotes | ~1 h for each ρ of the CRRA design decomposition |
| `horizonUS`, `horizonARG` | appendix D's checks | ~40 min plus ~4 h for the design; ~45 min |
| `--prepub` | the pre-publication checks (exact CRRA costs, stationary policies), read by no exhibit of the paper | ~9 h |
| `--data` | stage (0), with the network | minutes |

## Description of the code

`python/README.md` describes every file. In short:

- `python/paper/`: the pipeline: `replicate.py`, the stage drivers, `build.py` and its builders, `config.py` (every
  setting behind a published number) and `datasets.py` (the layout of `results/`).
- `python/InformalSavings/` and `python/US/`: the models and their drivers.
- `python/gridsearch/`: root finding and equilibrium selection on grids.

Docstrings cite the technical documentation's equations by label; `docs/equationLabels.csv` maps each label to
its number and page.

## Instructions to replicators

Set up, from the repository root (on Windows use `.venv\Scripts\python` and set `PYTHONUTF8=1`):

```
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip install -e python
```

**Rebuild every exhibit from the committed results (seconds):**

```
python python/paper/build.py --list        # every output buildable
python python/paper/build.py --no-copy     # rewrites results/paper; git diff shows no change but the figures' dates
```

**Re-solve everything (about a day):**

```
python python/paper/replicate.py --dry     # the commands, in order
python python/paper/replicate.py --list    # what exists
python python/paper/replicate.py           # all steps but stage (0); --prepub adds the checks, --data stage (0)
```

The stage drivers skip work whose output exists. A full re-solve therefore needs `--force`, or a fresh copy of
`results/` without the files to recompute. `python python/paper/compareResults.py` then reports every changed csv
against the committed version, aligned on the files' keys. Every csv that stages (i) and (ii) write counts the
equilibria and candidates at each state of its solve (`nEqMax`, `nCandMax`, `nFallback`). In the committed results
every count reads one equilibrium.

The [replicate](https://github.com/ChampionApe/SocialSecurityDesign/actions/workflows/replicate.yml) workflow
checks each commit on a clean Linux machine: it installs `requirements.txt`, rebuilds every exhibit and re-solves
one experiment from the committed calibration.

## List of tables and figures

Every exhibit of the paper and the supplemental appendix, its output in `python/paper/build.py`, its builder and
the files it reads. The rebuild command is `python python/paper/build.py --no-copy --only <output>`.

<!-- exhibits:begin -->
<!-- generated by public/export.py from docs/ (commit 0a3c2fa); do not edit -->
| Exhibit | Output | Built by | Reads |
|---|---|---|---|
| Figure 1 | `OECDdata` | `python/paper/oecdFigure1.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv` |
| Figure 2 | `coverageArg` | Not produced by this repository’s pipeline: drawn from Cetr'angolo and Grushka (2020), table 6. | -- |
| Table 1 | `ArgentinaCalibration` | `python/paper/tables.py` | `data/ArgentinaTest.xlsx`, `results/paper/calibrationSummary.csv` |
| Table 2 | `ArgentinaUniversal` | `python/paper/tables.py` | `data/ArgentinaTest.xlsx`, `results/calibration/informalSavings_rhoGridCommonX.csv`, `results/paper/calibrationSummary.csv`, `results/shocks/eeOnly_match_rho*_commonX.csv`, `results/shocks/universal_match_rho*_commonX.csv` |
| Figure 3 | `ARG_CRRA_LOG` | `python/paper/figures.py` | `data/ArgentinaTest.xlsx`, `results/calibration/informalSavings_rhoGridCommonX.csv`, `results/paper/calibrationSummary.csv`, `results/shocks/eeOnly_match_rho*_commonX.csv (16 files)`, `results/shocks/universal_match_rho*_commonX.csv (16 files)` |
| Table 3 | `USFRCalibration` | `python/paper/tablesUS.py` | `data/USMain_test.xlsx`, `results/paper/usCalibrationSummary.csv` |
| Figure 4 | `US_overview` | `python/paper/figuresUS.py` | `results/shocks/US_shocksCommonX.csv` |
| Figure 5 | `OECD26_Spending` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table 4 | `US_ESC_Calibration` | `python/paper/tablesUS.py` | `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv` |
| Figure 6 | `US_ESC_overview` | `python/paper/figuresUS.py` | `data/USMain_test.xlsx`, `results/esc_tau1/escExperiments.csv`, `results/paper/usCalibrationSummary.csv` |
| Figure 7 | `OECD26_Layer1` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26_tau1/anchor.csv`, `results/oecd26_tau1/layer1.csv` |
| Table 5 | `OECD26_Slopes` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26_tau1/anchor.csv`, `results/oecd26_tau1/layer1.csv` |
| Figure 8 | `OECD26_ThetaFit` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26_tau1/anchor.csv`, `results/oecd26_tau1/layer1.csv` |
| Table 6 | `NUM_Residuals` | `python/paper/tablesHorizon.py` | `results/numerical/ARG_residuals_commonX.csv`, `results/numerical/US_residuals.csv` |
| Table 7 | `US_householdheterogeneity` | `python/paper/tablesUS.py` | `results/paper/usCalibrationSummary.csv` |
| Table 8 | `FR_householdheterogeneity` | `python/paper/tablesUS.py` | `results/paper/usCalibrationSummary.csv` |
| Table 9 | `OECD26_Omega` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer2.csv` |
| Figure A1.1 | `ARG_LOG_FourInOne` | `python/paper/figures.py` | `data/ArgentinaTest.xlsx`, `results/paper/calibrationSummary.csv`, `results/sweeps/epsThetaGrid_rho*_commonX.csv` |
| Table A1.1 | `ARG_SpendingBaseline` | `python/paper/tablesOA.py` | `data/ArgentinaTest.xlsx`, `results/arg_spending052/calibration/informalSavings_rhoGridCommonX.csv`, `results/arg_spending052/shocks/universal_match_rho*_commonX.csv`, `results/calibration/informalSavings_rhoGridCommonX.csv`, `results/calibration/informalSavings_rhoGridCommonXEpsPre0.csv`, `results/paper/calibrationSummary.csv`, `results/shocks/universal_match_rho*_commonX.csv`, `results/shocks/universal_match_rho*_commonX_epsPre0.csv` |
| Figure A2.1 | `ARG_RhoGrid` | `python/paper/figuresOA.py` | `data/ArgentinaTest.xlsx`, `results/calibration/informalSavings_rhoGridCommonX.csv` |
| Figure A2.2 | `ARG_ReformPaths` | `python/paper/siteData.py` | `data/ArgentinaTest.xlsx`, `results/calibration/informalSavings_rhoGridCommonX.csv`, `results/paper/calibrationSummary.csv`, `results/shocks/eeOnly_match_rho*_commonX.csv (16 files)`, `results/shocks/universal_match_rho*_commonX.csv (16 files)` |
| Table A3.1 | `ARG_LeisureCalibration` | `python/paper/tablesOA.py` | `data/ArgentinaTest.xlsx`, `results/paper/calibrationSummary.csv` |
| Table A3.2 | `ARG_LeisureRho` | `python/paper/tablesOA.py` | `data/ArgentinaTest.xlsx`, `results/calibration/informalSavings_rhoGrid.csv`, `results/calibration/informalSavings_rhoGridCommonX.csv`, `results/paper/calibrationSummary.csv` |
| Figure B1.1 | `OECD_RhoGrid` | `python/paper/figuresOA.py` | `data/USMain_test.xlsx`, `results/calibration/FR_rhoGridCommonX.csv`, `results/calibration/UK_rhoGridCommonX.csv`, `results/calibration/US_rhoGridCommonX.csv` |
| Figure B2.1 | `OECD_Counterfactuals` | `python/paper/siteData.py` | `results/paper/usCalibrationSummary.csv`, `results/shocks/US_shocks.csv`, `results/shocks/US_shocksCommonX.csv`, `results/shocks/US_shocksXiCommonX.csv` |
| Figure B2.2 | `US_TaxState` | `python/paper/figuresOA.py` | `data/USMain_test.xlsx`, `results/numerical/C5_taxPolicyState.csv` |
| Table B2.1 | `US_Strategic` | `python/paper/tablesOA.py` | `data/USMain_test.xlsx`, `results/numerical/C5_strategicTerm.csv` |
| Table B3.1 | `OECD_LeisureCalibration` | `python/paper/tablesOA.py` | `results/paper/usCalibrationSummary.csv` |
| Table B3.2 | `OECD_LeisureGroups` | `python/paper/tablesOA.py` | `results/paper/usCalibrationSummary.csv` |
| Figure B3.1 | `OECD_Counterfactuals` | `python/paper/siteData.py` | `results/paper/usCalibrationSummary.csv`, `results/shocks/US_shocks.csv`, `results/shocks/US_shocksCommonX.csv`, `results/shocks/US_shocksXiCommonX.csv` |
| Table B4.1 | `OECD_Calibration` | `python/paper/tablesUS.py` | `data/USMain_test.xlsx`, `results/paper/usCalibrationSummary.csv` |
| Table B4.2 | `UK_householdheterogeneity` | `python/paper/tablesUS.py` | `results/paper/usCalibrationSummary.csv` |
| Figure C1.1 | `ESC_Paths` | `python/paper/siteData.py` | `data/USMain_test.xlsx`, `results/esc/escCalibration.csv`, `results/esc/escCalibrationCRRA.csv`, `results/esc/escPath.csv`, `results/esc/escPathCRRA.csv`, `results/esc/escPathHost.csv`, `results/esc/escXiRobustness.csv`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv`, `results/esc_tau1/escPath.csv`, `results/esc_tau1/escPathCRRA.csv`, `results/esc_tau1/escPathHost.csv`, `results/esc_tau1/escXiRobustness.csv`, `results/paper/usCalibrationSummary.csv` |
| Figure C1.2 | `ESC_Paths` | `python/paper/siteData.py` | `data/USMain_test.xlsx`, `results/esc/escCalibration.csv`, `results/esc/escCalibrationCRRA.csv`, `results/esc/escPath.csv`, `results/esc/escPathCRRA.csv`, `results/esc/escPathHost.csv`, `results/esc/escXiRobustness.csv`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv`, `results/esc_tau1/escPath.csv`, `results/esc_tau1/escPathCRRA.csv`, `results/esc_tau1/escPathHost.csv`, `results/esc_tau1/escXiRobustness.csv`, `results/paper/usCalibrationSummary.csv` |
| Figure C1.3 | `US_ESC_Path_linear` | `python/paper/figuresOA.py` | `data/USMain_test.xlsx`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv`, `results/esc_tau1/escPath.csv` |
| Figure C1.4 | `US_ESC_Path` | `python/paper/figuresOA.py` | `data/USMain_test.xlsx`, `results/esc/escCalibration.csv`, `results/esc/escCalibrationCRRA.csv`, `results/esc/escPath.csv` |
| Table C2.1 | `ESC_Timing` | `python/paper/tablesOA.py` | `data/USMain_test.xlsx`, `results/esc/escCalibration.csv`, `results/numerical/C5_criticalTurnout.csv` |
| Table C2.2 | `ESC_CostCalibration_linear` | `python/paper/tablesOA.py` | `data/USMain_test.xlsx`, `results/esc_tau1/escCalibration.csv` |
| Table C2.3 | `ESC_CostCalibration` | `python/paper/tablesOA.py` | `data/USMain_test.xlsx`, `results/esc/escCalibration.csv` |
| Figure C2.1 | `ESC_Counterfactuals` | `python/paper/siteData.py` | `data/USMain_test.xlsx`, `results/esc/escExperiments.csv`, `results/esc/escXiRobustness.csv`, `results/esc_tau1/escExperiments.csv`, `results/esc_tau1/escXiRobustness.csv`, `results/paper/usCalibrationSummary.csv` |
| Figure C2.2 | `US_ESC_TauTheta_linear` | `python/paper/figuresOA.py` | `data/USMain_test.xlsx`, `data/oecd26_theta.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escTauTheta.csv` |
| Figure C2.3 | `US_ESC_TauTheta` | `python/paper/figuresOA.py` | `data/USMain_test.xlsx`, `data/oecd26_theta.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/esc/escCalibration.csv`, `results/esc/escTauTheta.csv` |
| Table C2.4 | `ESC_Decompositions_linear` | `python/paper/tablesOA.py` | `data/USMain_test.xlsx`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv`, `results/esc_tau1/escExperiments.csv`, `results/numerical_tau1/C5_designDecomp.csv`, `results/numerical_tau1/C5_designDecompCRRA_rho*.csv (2 files)`, `results/numerical_tau1/C5_lambdaElasticity.csv` |
| Table C2.5 | `ESC_Decompositions` | `python/paper/tablesOA.py` | `data/USMain_test.xlsx`, `results/esc/escCalibration.csv`, `results/esc/escCalibrationCRRA.csv`, `results/esc/escExperiments.csv`, `results/numerical/C5_designDecomp.csv`, `results/numerical/C5_designDecompCRRA_rho*.csv (2 files)`, `results/numerical/C5_lambdaElasticity.csv` |
| Table C2.6 | `ESC_TauPower` | `python/paper/tablesOA.py` | `data/USMain_test.xlsx`, `results/esc/escCalibration.csv`, `results/esc/escCalibrationCRRA.csv`, `results/esc/escShocks.csv`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv`, `results/esc_tau1/escShocks.csv` |
| Figure C3.1 | `ESC_Counterfactuals` | `python/paper/siteData.py` | `data/USMain_test.xlsx`, `results/esc/escExperiments.csv`, `results/esc/escXiRobustness.csv`, `results/esc_tau1/escExperiments.csv`, `results/esc_tau1/escXiRobustness.csv`, `results/paper/usCalibrationSummary.csv` |
| Figure C3.2 | `ESC_Paths` | `python/paper/siteData.py` | `data/USMain_test.xlsx`, `results/esc/escCalibration.csv`, `results/esc/escCalibrationCRRA.csv`, `results/esc/escPath.csv`, `results/esc/escPathCRRA.csv`, `results/esc/escPathHost.csv`, `results/esc/escXiRobustness.csv`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv`, `results/esc_tau1/escPath.csv`, `results/esc_tau1/escPathCRRA.csv`, `results/esc_tau1/escPathHost.csv`, `results/esc_tau1/escXiRobustness.csv`, `results/paper/usCalibrationSummary.csv` |
| Table C4.1 | `OECD26_Variants` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26/layer2.csv`, `results/oecd26_theta_A/anchor.csv`, `results/oecd26_theta_A/calibration/US_rhoGridCommonX.csv`, `results/oecd26_theta_A/layer1.csv`, `results/oecd26_theta_A/layer2.csv`, `results/oecd26_theta_B/anchor.csv`, `results/oecd26_theta_B/calibration/US_rhoGridCommonX.csv`, `results/oecd26_theta_B/layer1.csv`, `results/oecd26_theta_B/layer2.csv` |
| Figure C4.1 | `OECD26_ThetaFitPopulation` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26/layer2.csv` |
| Figure C4.2 | `OECD26_Layer1_nuNone` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table C4.2 | `OECD26_Omega_nuNone` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer2.csv` |
| Figure C4.3 | `OECD26_Layer1_nuExit` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table C4.3 | `OECD26_Omega_nuExit` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer2.csv` |
| Figure C4.4 | `OECD26_ThetaFitMeasure` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26/layer2.csv`, `results/oecd26_theta_A/anchor.csv`, `results/oecd26_theta_A/calibration/US_rhoGridCommonX.csv`, `results/oecd26_theta_A/layer1.csv`, `results/oecd26_theta_A/layer2.csv`, `results/oecd26_theta_B/anchor.csv`, `results/oecd26_theta_B/calibration/US_rhoGridCommonX.csv`, `results/oecd26_theta_B/layer1.csv`, `results/oecd26_theta_B/layer2.csv` |
| Figure C4.5 | `OECD26_Layer1_thetaA` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/oecd26_theta_A/anchor.csv`, `results/oecd26_theta_A/calibration/US_rhoGridCommonX.csv`, `results/oecd26_theta_A/layer1.csv` |
| Table C4.4 | `OECD26_Omega_thetaA` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/oecd26_theta_A/anchor.csv`, `results/oecd26_theta_A/calibration/US_rhoGridCommonX.csv`, `results/oecd26_theta_A/layer2.csv` |
| Figure C4.6 | `OECD26_Layer1_thetaB` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/oecd26_theta_B/anchor.csv`, `results/oecd26_theta_B/calibration/US_rhoGridCommonX.csv`, `results/oecd26_theta_B/layer1.csv` |
| Table C4.5 | `OECD26_Omega_thetaB` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/oecd26_theta_B/anchor.csv`, `results/oecd26_theta_B/calibration/US_rhoGridCommonX.csv`, `results/oecd26_theta_B/layer2.csv` |
| Figure C4.7 | `OECD26_Layer1_cses` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table C4.6 | `OECD26_Slopes_cses` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26/layer2.csv` |
| Figure C4.8 | `OECD26_Layer1_omegaUS` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table C4.7 | `OECD26_Slopes_omegaUS` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26/layer2.csv` |
| Table C4.8 | `OECD26_GivenDesign` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table C4.9 | `OECD26_USInputs_linear` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26_tau1/anchor.csv`, `results/oecd26_tau1/layer1.csv` |
| Table C4.10 | `OECD26_USInputs` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table C4.11 | `OECD26_SlopesNoUS` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table C4.12 | `OECD26_SlopesTau1` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26_tau1/anchor.csv`, `results/oecd26_tau1/layer1.csv` |
| Table C5.1 | `US_ESC_Country_linear` | `python/paper/tablesUS.py` | `results/esc_tau1/escCountry.csv` |
| Table C5.2 | `US_ESC_Country` | `python/paper/tablesUS.py` | `results/esc/escCountry.csv` |
| Table C5.3 | `UK_ESC_Calibration_linear` | `python/paper/tablesUS.py` | `results/esc_tau1/escCalibrationCRRAUK.csv`, `results/esc_tau1/escCountry.csv` |
| Table C5.4 | `UK_ESC_Calibration` | `python/paper/tablesUS.py` | `results/esc/escCalibrationCRRAUK.csv`, `results/esc/escCountry.csv` |
| Table D1.1 | `OECD_Countries` | `python/paper/tablesOA.py` | `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv` |
| Table D1.2 | `OECD_Sources` | `python/paper/tablesOA.py` | `data/oecd26_pensionTiers_harmonised.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecd29_sources.csv`, `data/oecdFigure1.csv`, `data/oecdFigure1_sources.csv` |
| Table D1.3 | `OECD_Correlations` | `python/paper/tablesOA.py` | `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/paper/oecdCorrelations.csv` |
| Figure D1.2 | `OECD_InequalityMeasures` | `python/paper/figuresOA.py` | `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv` |
| Table D1.4 | `OECD_InequalityCorrelations` | `python/paper/tablesOA.py` | `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/paper/oecdInequalityCorrelations.csv` |
| Figure D1.3 | `OECD26_Layer1Inequality` | `python/paper/figuresOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv` |
| Table D1.5 | `OECD26_Layer1InequalitySlopes` | `python/paper/tablesOECD26.py` | `data/oecd26_sample.csv`, `data/oecd26_theta.csv`, `data/oecd29_distribution.csv`, `data/oecdFigure1.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/paper/oecd26Layer1InequalitySlopes.csv` |
| Table D2.1 | `NUM_Stationary` | `python/paper/tablesOA.py` | `results/numerical/ARG_stationaryApprox_commonX.csv`, `results/numerical/US_ESC_stationaryApprox.csv`, `results/numerical/US_stationaryApprox.csv` |
| Table D2.2 | `NUM_Selection` | `python/paper/tablesOA.py` | `results/arg_spending052/calibration/informalSavings_rhoGridCommonX.csv`, `results/arg_spending052/shocks/universal_match_rho*_commonX.csv`, `results/calibration/FR_rhoGrid.csv`, `results/calibration/FR_rhoGridCommonX.csv`, `results/calibration/UK_rhoGrid.csv`, `results/calibration/UK_rhoGridCommonX.csv`, `results/calibration/US_rhoGrid.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/calibration/informalSavings_rhoGrid.csv`, `results/calibration/informalSavings_rhoGridCommonX.csv`, `results/esc/escCalibration.csv`, `results/esc/escCalibrationCRRA.csv`, `results/esc/escCalibrationCRRAUK.csv`, `results/esc/escCalibrationCosts.csv`, `results/esc/escCountry.csv`, `results/esc/escPath.csv`, `results/esc/escPathCRRA.csv`, `results/esc/escPathHost.csv`, `results/esc/escShocks.csv`, `results/esc/escShocksCRRA.csv`, `results/esc/escShocksCosts.csv`, `results/esc/escTauTheta.csv`, `results/esc/escXiRobustness.csv`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv`, `results/esc_tau1/escCalibrationCRRAUK.csv`, `results/esc_tau1/escCalibrationCosts.csv`, `results/esc_tau1/escCountry.csv`, `results/esc_tau1/escPath.csv`, `results/esc_tau1/escPathCRRA.csv`, `results/esc_tau1/escPathHost.csv`, `results/esc_tau1/escShocks.csv`, `results/esc_tau1/escShocksCRRA.csv`, `results/esc_tau1/escShocksCosts.csv`, `results/esc_tau1/escTauTheta.csv`, `results/esc_tau1/escXiRobustness.csv`, `results/numerical/ARG_horizon_commonX.csv`, `results/numerical/ARG_stationaryApprox_commonX.csv`, `results/numerical/US_ESC_stationaryApprox.csv`, `results/numerical/US_derivativeScheme.csv`, `results/numerical/US_horizon.csv`, `results/numerical/US_horizonDesign.csv`, `results/numerical/US_stationaryApprox.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26/layer2.csv`, `results/oecd26_tau1/anchor.csv`, `results/oecd26_tau1/layer1.csv`, `results/oecd26_tau1/layer2.csv`, `results/oecd26_theta_A/anchor.csv`, `results/oecd26_theta_A/layer1.csv`, `results/oecd26_theta_A/layer2.csv`, `results/oecd26_theta_B/anchor.csv`, `results/oecd26_theta_B/layer1.csv`, `results/oecd26_theta_B/layer2.csv`, `results/shocks/US_shocks.csv`, `results/shocks/US_shocksCommonX.csv`, `results/shocks/US_shocksXiCommonX.csv`, `results/shocks/eeOnly_match_rho*.csv (16 files)`, `results/shocks/eeOnly_match_rho*_commonX.csv (16 files)`, `results/shocks/universal_match_rho*.csv (16 files)`, `results/shocks/universal_match_rho*_commonX.csv (16 files)`, `results/sweeps/epsThetaGrid_rho*.csv`, `results/sweeps/epsThetaGrid_rho*_commonX.csv` |
| Table D2.3 | `NUM_SelectionFiles` | `python/paper/tablesOA.py` | `results/arg_spending052/calibration/informalSavings_rhoGridCommonX.csv`, `results/arg_spending052/shocks/universal_match_rho*_commonX.csv`, `results/calibration/FR_rhoGrid.csv`, `results/calibration/FR_rhoGridCommonX.csv`, `results/calibration/UK_rhoGrid.csv`, `results/calibration/UK_rhoGridCommonX.csv`, `results/calibration/US_rhoGrid.csv`, `results/calibration/US_rhoGridCommonX.csv`, `results/calibration/informalSavings_rhoGrid.csv`, `results/calibration/informalSavings_rhoGridCommonX.csv`, `results/esc/escCalibration.csv`, `results/esc/escCalibrationCRRA.csv`, `results/esc/escCalibrationCRRAUK.csv`, `results/esc/escCalibrationCosts.csv`, `results/esc/escCountry.csv`, `results/esc/escPath.csv`, `results/esc/escPathCRRA.csv`, `results/esc/escPathHost.csv`, `results/esc/escShocks.csv`, `results/esc/escShocksCRRA.csv`, `results/esc/escShocksCosts.csv`, `results/esc/escTauTheta.csv`, `results/esc/escXiRobustness.csv`, `results/esc_tau1/escCalibration.csv`, `results/esc_tau1/escCalibrationCRRA.csv`, `results/esc_tau1/escCalibrationCRRAUK.csv`, `results/esc_tau1/escCalibrationCosts.csv`, `results/esc_tau1/escCountry.csv`, `results/esc_tau1/escPath.csv`, `results/esc_tau1/escPathCRRA.csv`, `results/esc_tau1/escPathHost.csv`, `results/esc_tau1/escShocks.csv`, `results/esc_tau1/escShocksCRRA.csv`, `results/esc_tau1/escShocksCosts.csv`, `results/esc_tau1/escTauTheta.csv`, `results/esc_tau1/escXiRobustness.csv`, `results/numerical/ARG_horizon_commonX.csv`, `results/numerical/ARG_stationaryApprox_commonX.csv`, `results/numerical/US_ESC_stationaryApprox.csv`, `results/numerical/US_derivativeScheme.csv`, `results/numerical/US_horizon.csv`, `results/numerical/US_horizonDesign.csv`, `results/numerical/US_stationaryApprox.csv`, `results/oecd26/anchor.csv`, `results/oecd26/layer1.csv`, `results/oecd26/layer2.csv`, `results/oecd26_tau1/anchor.csv`, `results/oecd26_tau1/layer1.csv`, `results/oecd26_tau1/layer2.csv`, `results/oecd26_theta_A/anchor.csv`, `results/oecd26_theta_A/layer1.csv`, `results/oecd26_theta_A/layer2.csv`, `results/oecd26_theta_B/anchor.csv`, `results/oecd26_theta_B/layer1.csv`, `results/oecd26_theta_B/layer2.csv`, `results/shocks/US_shocks.csv`, `results/shocks/US_shocksCommonX.csv`, `results/shocks/US_shocksXiCommonX.csv`, `results/shocks/eeOnly_match_rho*.csv (16 files)`, `results/shocks/eeOnly_match_rho*_commonX.csv (16 files)`, `results/shocks/universal_match_rho*.csv (16 files)`, `results/shocks/universal_match_rho*_commonX.csv (16 files)`, `results/sweeps/epsThetaGrid_rho*.csv`, `results/sweeps/epsThetaGrid_rho*_commonX.csv` |
| Table D2.4 | `NUM_HorizonPaths` | `python/paper/tablesHorizon.py` | `data/ArgentinaTest.xlsx`, `data/USMain_test.xlsx`, `results/numerical/ARG_horizon_commonX.csv`, `results/numerical/US_horizon.csv`, `results/numerical/US_horizonDesign.csv` |
<!-- exhibits:end -->

## Known gaps

- Figure 5.1 (beneficiaries of Argentina's pension system) is drawn from Cetrángolo and Grushka (2020), table 6.
  Its data and script are not in the repository.
- The household statistics of `ArgentinaTest.xlsx` (EPH), `FRMain.xlsx` (LIS) and `UKMain.xlsx` (BHPS) were
  tabulated outside this repository, and the tabulation code is not included.

## License and citation

Code: MIT (`LICENSE`). The authors' data compilations and documentation: CC BY 4.0 (`LICENSE-data`). To cite the
paper or this package, see `CITATION.cff`.
