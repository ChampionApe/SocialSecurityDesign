r""" Loader of the U.S. workbook data/USMain_test.xlsx: `pars` and `kwargs` for ModelUS, the dates, the workweek.

Importing this module has side effects the suites and drivers rely on: `pars` is updated with the LOG anchor point
(ρ = 1, ω = 2, β = 0.6), so every `testmod.pars | {...}` carries them unless it overrides, and `mLOG`, a ModelUS on
`pars` with db['dates'] and db['workweek'] installed, is built at import (the structure suites read it).
"""
import os, pandas as pd, numpy as np
from model import ModelUS

# Read in data from excel: data/USMain_test.xlsx
PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'USMain_test.xlsx')
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

# Create objects to parse to model initialization:
kwargs = {'T': TLog+t_ss, 'nj':dfj.shape[0]}
pars = {'α': dfc['Capital income share'], 'ξ': dfc['Labor supply elasticity'], 'ν': np.hstack([νLog, np.full(kwargs['T']-TLog, νLog[-1])]),
        'τ0': dfc['Pension tax'], 'RR0': dfc['Replacement rate'], 'RRGroups': tuple(eval(dfc['Replacement rate groups'])), 't0': dates.get_loc(dfc['Calibration year']),
        'γj': dfj['Population shares'].values.astype(float),
        'μj': dfj['Voting shares'].values.astype(float),
        'Xj': 1,
        'zxj': (dfj['Hours'].values/(dfj['Hours'].values.mean())).astype(float),
        'zηj': (dfj['Income'].values/(dfj['Income'].values.mean())).astype(float)}
pars['h0'] = workweek / (7*12) # estimate of share of time spend on labor
pars['R0'] = dfc['30y interest']
pars['Ushare0'] = dfc['Universal share']


# Initialize LOG model with some parameter values:
pars.update({'ρ': 1, 'ω': 2, 'β': .6})
mLOG = ModelUS(pars = pars, **kwargs)
mLOG.db['dates'] = dates
mLOG.db['workweek'] = workweek
