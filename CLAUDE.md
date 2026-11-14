# Social Security Design and Its Political Support: replication package

This repository reproduces every table and figure of Berg and Gonzalez-Eiras (2026), *Social Security Design and
Its Political Support*, and of its supplemental appendix (https://championape.github.io/SocialSecurityDesign).
The paper studies pension systems in an overlapping-generations economy with heterogeneous households. Each
period, a government that cannot commit chooses the pension tax by probabilistic voting. In section 7 the
electorate also chooses the design θ, the earnings link of benefits: θ = 1 is fully Bismarckian, θ = 0 fully
Beveridgean (flat). There are two quantitative arms. Argentina's 2005–10 extension of pensions to the informal
sector uses `python/InformalSavings` (section 5). The U.S., France, the UK and a cross-section of 21 OECD
countries use `python/US` (sections 6 and 7).

The authors develop the code in a private repository. This one is a snapshot of it, and each commit names the
private revision it came from. Do not expect a history of the code here.

## Map

- `README.md`: the replication instructions (AEA template): data sources and terms, requirements, run times, and
  the table from every exhibit to its builder and inputs. Read it first.
- `python/README.md`: what each code file does. `python/paper/` is the pipeline; `python/InformalSavings/`,
  `python/US/` hold the models; `python/gridsearch/` is a numerical package.
- `data/README.md`: every input file, its source and what reads it.
- `results/`: committed outputs of every step, the pickled calibrated instances among them. `results/paper/`
  holds the built tables (`.tex`), figures (`.pdf`, `.svg`) and the site's datasets.
- `docs/`: `PaperExhibits.pdf` and `SupplementalExhibits.pdf` (every exhibit as built, each with its output id,
  builder, inputs and rebuild command), `TechnicalDocumentation.pdf` (the models, equilibria and algorithms;
  the paper cites it), and `equationLabels.csv`.

## Setup

Python 3.14 (3.14.6 was used): `python -m venv .venv`, then `pip install -r requirements.txt` and
`pip install -e python`. On Windows set `PYTHONUTF8=1` for every run, because the scripts print Greek letters.
Run every script from the repository root.

## Ground rules

1. `results/` is output and `data/` is input. Never edit a file under `results/` by hand. A generated `.tex`
   carries a `%% GENERATED` banner; change its source and rebuild instead. The sources are
   `python/paper/config.py`, which records every setting behind a published number, and the experiment that
   writes the csv.
2. Stage (iii) is cheap and safe: `python python/paper/build.py --no-copy` (or `--only <output>`) rebuilds from
   `results/` in seconds. Always pass `--no-copy`, since this repository has no `writing/`. `build.py --list`
   reports what is buildable.
3. Stages (i) and (ii) and the other model runs take hours (`README.md`, run times). Ask the user before
   launching one. Show the plan first with `python python/paper/replicate.py --dry` or `--list`, or with the
   driver's own `--dry`/`--list`. Drivers skip work whose output exists and `--force` redoes it, so a forced run
   overwrites committed results.
4. For an experiment of your own (another parameter, another country), write to new files. Never overwrite the
   committed csvs. `python python/paper/compareResults.py` lists what changed under `results/` against `HEAD`.
5. Never write the `data/*.xlsx` workbooks with openpyxl: it drops the cached formula values and every steady
   state fails. Stage (0) (`--data`) reaches the network; run it only when asked.

## Answering questions

- Where a number in the paper comes from: the exhibit table in `README.md` (or the provenance block under each
  exhibit in `docs/PaperExhibits.pdf`) gives the output id and the csvs. The generated `.tex` names its source in
  its `%% Source:` banner, and each table row ends in a `% row: <key>` comment.
- What equation a function implements: docstrings cite labels such as `eq:governmentBudget` or
  `alg:LOG:gridsearch`. Look them up in `docs/equationLabels.csv` (number, page, part) and read that page of
  `docs/TechnicalDocumentation.pdf`. Its four parts are the informal analytical model, the informal savings model
  (Argentina), the U.S. model, and the OECD data. Rows of part `paper` are the paper's own labels (`sec:esc`,
  numbered as in the paper), which the code also cites. A label found in several parts belongs to the code's own
  folder: `python/US` → `US`, `python/InformalSavings` → `informalSavings`.
- How a solved model looks: load a pickled instance from `results/calibration/instances*/` under the matching
  model folder's `sys.path`, as `python/US/runShocksUS.py` does. Pickles need the pinned package versions.
- Why something is done a certain way: comments cite numbered findings ("#13", "finding #7") and notes, logs or
  archives of the authors' private repository (`notes/`, `archive/`, `logs/`, `RESEARCH_LOG.md`, `writing/…`,
  labels like A14, C5, W25). Those files are not part of this repository, so do not search for them. The code
  and the technical documentation are self-contained without them.

## Glossary

θ: the design, the earnings-related share of benefits. τ: the pension tax rate. ε: Argentina's universal pension
(the informal type's benefit). ρ: the inverse of the elasticity of intertemporal substitution (LOG is ρ = 1, CRRA
the rest). X: the disutility of hours; `commonX` calibrates one X for all types, the headline. EE: the economic
equilibrium at given policy; PEE: the politico-economic equilibrium. ESC: the endogenous choice of the design
(section 7); its cost of redistribution has slope κ (`sizePiR`). The `oecd26` files belong to the cross-country
test, with layer 1 'the model' and layer 2 'own ω'. Variant C of `data/oecd26_theta.csv` is the paper's measure
of θ (`--theta C --main`). `c5Numbers.py` and the `C5_*.csv` files are named after an internal work item: they hold
the smaller numbers the text quotes.
