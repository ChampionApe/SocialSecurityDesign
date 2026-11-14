import os, pandas as pd, numpy as np
from model import ModelInformalSavings

# Read in data from excel: data/ArgentinaTest.xlsx
PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'ArgentinaTest.xlsx')
wb = pd.read_excel(PATH, sheet_name = None, header = None)
t_ss  = 3 # number of periods in steady state

# Read population sheet:
dft = pd.DataFrame(wb['population'].values)
dft = pd.DataFrame(dft.iloc[1:,].values, columns = dft.iloc[0,:]).set_index('t')
datesLog = dft.index
νLog = dft['Worker-to-retiree'].values.astype(float)
TLog = len(datesLog)
dates = datesLog.union(pd.Index([datesLog[-1]+30*i for i in range(1, t_ss+1)]))

# Read heterogeneity sheet:
dfj = pd.DataFrame(wb['heterogeneity'].values)
dfj = pd.DataFrame(dfj.iloc[1:,].values, columns = dfj.iloc[0,:])

# Read calibration sheet:
dfc = pd.DataFrame(wb['calibration'].values)
dfc = pd.Series(dfc.iloc[1,:].values, index = dfc.iloc[0,:].values)
workweek = dfc['Average workweek']

# Income and hours relative to the average FORMAL household, eq (calibration:z): the population-weighted
# mean over j>0, with the informal row excluded from the denominator. The formal z_i only need the right
# ratios (the eigenvector problem is scale-free), but z_0 is read as a level by eq (calibration:eta0/X0).
formalMean = lambda x: np.average(x[1:].astype(float), weights = dfj['Population shares'].values[1:].astype(float))

# Create objects to parse to model initialization:
kwargs = {'T': TLog+t_ss, 'nj':dfj.shape[0]}
pars = {'α': dfc['Capital income share'], 'ξ': dfc['Labor supply elasticity'], 'ν': np.hstack([νLog, np.full(kwargs['T']-TLog, νLog[-1])]),
        'RR0': dfc['Replacement rate'], 'RRGroups': tuple(eval(dfc['Replacement rate groups'])), 't0': dates.get_loc(dfc['Calibration year']),
        'γj': dfj['Population shares'].values.astype(float), 
        'μj': dfj['Voting shares'].values.astype(float), 
        'Xj': 1,
        'zxj': (dfj['Hours'].values/formalMean(dfj['Hours'].values)).astype(float),
        'zηj': (dfj['Income'].values/formalMean(dfj['Income'].values)).astype(float)}
pars['h0'] = workweek / (7*12) # estimate of share of time spend on labor
pars['s0'] = dfc['Savings rate']   # reported only; the target that identifies beta is KY0 below

# The tax target is DERIVED from pension spending over GDP (the workbook's datum), with GDP formal plus
# informal output and the informal share y the model's own at the calibration (ModelInformalSavings.
# taxTarget): it follows alpha, gamma_0 and z_0^eta. A literal tau0 would silently go stale with any of them.
pars['pensionSpending0'] = float(dfc['Pension spending'])
pars['τ0'], pars['yInformal0'] = ModelInformalSavings.taxTarget(
    pars['pensionSpending0'], float(dfc['Capital income share']), float(pars['γj'][0]), float(pars['zηj'][0]))

# The capital-output target does NOT live in the workbook: it is a window mean over an external series,
# so python/paper/dataTargets.py derives it and this csv is the record (target, window, source, date).
# The workbook's 'Savings rate' is what it superseded -- archive/notes/argentina_calibrationTarget.md.
TARGETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data',
                       'argentina_calibrationTargets.csv')
pars['KY0'] = float(pd.read_csv(TARGETS).set_index('target').loc['capitalOutputRatio', 'value'])



# Initialize LOG model with some parameter values:
pars.update({'ρ': 1, 'ω': 2, 'β': .6})
mLOG = ModelInformalSavings(pars = pars, **kwargs)
mLOG.db['dates'] = dates
mLOG.db['workweek'] = workweek

