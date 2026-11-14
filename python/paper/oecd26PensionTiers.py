r""" Stage (0), DRAFT: harmonised cross-country data on pension tiers, first-tier benefits and the financing
of each tier, for the 26 OECD countries of the cross-section. Not wired into build.py, REPLICATION.md, the
tests or any calibration; notes/data_thetaTau.md section A reads it.

Run:  .venv\Scripts\python.exe -I python\paper\oecd26PensionTiers.py [--force] [--refetch]
      -> data/oecd26_pensionTiers_harmonised.csv (long: iso3, item, tier, earnings_multiple, sex, value, unit,
         year, source, url, retrieved, note)

`main` is the only code here that touches the network: `fetch` downloads every raw file missing from
logs/thetaTau_raw/harmonised/<source>/ (gitignored), then `collect` reads them. main skips when the output
exists unless --force; raw files are reused unless --refetch. Raw files are untrusted data: never executed,
nothing written into their directories. A raw file's mtime is its retrieval date.

Sources: Pensions at a Glance (PaG) 2021 and 2025 StatLinks (publication tables, as published), the OECD
Data Explorer (PaG 2025 replacement rates DF_PRR; SOCX aggregates; asset-backed pensions DF_BC/DF_PA;
national-accounts table 29), the SOCX programme-level workbook (webfs.oecd.org), Eurostat ESSPROS
(spr_exp_pens, spr_exp_fol, spr_exp_fsu, spr_rec_sumt, spr_pns_ben) and the PaG 2013 progressivity table.
Values are copied as published: no interpolation, no filling, a missing cell is no row.
"""
import os, re, sys, argparse, datetime, urllib.request
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

C26 = ['AUS', 'AUT', 'BEL', 'CAN', 'CZE', 'DNK', 'FIN', 'FRA', 'DEU', 'GRC', 'HUN', 'ISL', 'IRL', 'ITA',
       'JPN', 'LUX', 'NLD', 'NZL', 'NOR', 'POL', 'PRT', 'ESP', 'SWE', 'CHE', 'GBR', 'USA']
NAME = {'Australia': 'AUS', 'Austria': 'AUT', 'Belgium': 'BEL', 'Canada': 'CAN', 'Czech Republic': 'CZE',
        'Czechia': 'CZE', 'Denmark': 'DNK', 'Finland': 'FIN', 'France': 'FRA', 'Germany': 'DEU',
        'Greece': 'GRC', 'Hungary': 'HUN', 'Iceland': 'ISL', 'Ireland': 'IRL', 'Italy': 'ITA', 'Japan': 'JPN',
        'Luxembourg': 'LUX', 'Netherlands': 'NLD', 'New Zealand': 'NZL', 'Norway': 'NOR', 'Poland': 'POL',
        'Portugal': 'PRT', 'Spain': 'ESP', 'Sweden': 'SWE', 'Switzerland': 'CHE', 'United Kingdom': 'GBR',
        'United States': 'USA'}
# Eurostat geo codes of the EU/EFTA members of the 26 (and the UK, in ESSPROS until it left)
ESTAT = {'AT': 'AUT', 'BE': 'BEL', 'CZ': 'CZE', 'DK': 'DNK', 'FI': 'FIN', 'FR': 'FRA', 'DE': 'DEU', 'EL': 'GRC',
         'HU': 'HUN', 'IS': 'ISL', 'IE': 'IRL', 'IT': 'ITA', 'LU': 'LUX', 'NL': 'NLD', 'NO': 'NOR', 'PL': 'POL',
         'PT': 'PRT', 'ES': 'ESP', 'SE': 'SWE', 'CH': 'CHE', 'UK': 'GBR'}
OUT = os.path.join(C.DATA, 'oecd26_pensionTiers_harmonised.csv')
RAWDIR = os.path.join(C.REPO, 'logs', 'thetaTau_raw', 'harmonised')
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) '
                    'Chrome/124.0 Safari/537.36'}
SDMX = 'https://sdmx.oecd.org/public/rest/data/'
ESDMX = 'https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/'
K26, KEU = '+'.join(C26), '+'.join(ESTAT)
Y0, Y1 = 2015, 2022


def _key(*parts):
    return '.'.join(parts)


# PaG StatLinks: code -> (edition, table as published). Edition document ids: 2021 ca401ebd-en, 2025 e40274c1-en.
STATLINK = {
    # PaG 2021 (pension rules in 2020; "Version 1 - Last updated: 08-Dec-2021")
    'b2f0ws': (2021, 'Table 4.1 Gross pension replacement rates by earnings, mandatory schemes'),
    'uwc8nt': (2021, 'Table 4.2 Gross pension replacement rates from mandatory public, mandatory private and voluntary private schemes'),
    'cofr31': (2021, 'Table 4.4 Net pension replacement rates by earnings'),
    'scxdj3': (2021, 'Table 4.5 Gross and net pension replacement rates from mandatory and voluntary schemes'),
    '2xd6uo': (2021, 'Table 3.1 Structure of retirement-income provision through mandatory schemes'),
    'uxgr8v': (2021, 'Table 3.2 Current level and recipients of first-tier benefits'),
    'kpoilu': (2021, 'Figure 3.2 Non-contributory first-tier benefits'),
    'k23yri': (2021, 'Figure 3.3 Contributory first-tier benefits'),
    'h4uwj7': (2021, 'Table 8.1 Mandatory pension contribution rates in 2020'),
    '9oqgj4': (2021, 'Table 8.2 Social insurance contribution rates in 2020'),
    '92exj3': (2021, 'Table 8.3 Cash expenditure on old-age and survivors benefits: public schemes'),
    'kp4d5j': (2021, 'Table 8.4 Cash expenditures on old-age and survivor benefits: private schemes'),
    '8gtspo': (2021, 'Figure 9.1 Minimum or mandatory contribution rates in mandatory and auto-enrolment plans'),
    'gkvjah': (2021, 'Table 9.1 Coverage of retirement savings plans'),
    'qi2ev9': (2021, 'Table 9.2 Assets in retirement savings plans and public pension reserve funds'),
    'lqhcun': (2021, 'Figure 7.1 Income sources of older people'),
    # PaG 2025 (pension rules in 2024), codes read off the publication pdf
    'oeavq4': (2025, 'Table 4.1 Gross pension replacement rates by earnings, mandatory schemes'),
    'pru4oj': (2025, 'Table 4.2 Gross pension replacement rates from mandatory public, mandatory private and voluntary private schemes'),
    'ic8ung': (2025, 'Table 4.4 Net pension replacement rates by earnings'),
    'f2v8nc': (2025, 'Table 4.5 Gross and net pension replacement rates from mandatory and voluntary schemes'),
    'mh9ruy': (2025, 'Table 3.1 Structure of retirement-income provision through mandatory schemes'),
    'i28xu0': (2025, 'Table 3.2 Current level and recipients of first-tier benefits'),
    's5w1y7': (2025, 'Figure 3.2 Non-contributory first-tier benefits'),
    '7ohtnu': (2025, 'Figure 3.3 Contributory first-tier benefits'),
    'orhufb': (2025, 'Table 8.1 Mandatory contribution rates in 2024'),
    '92ur17': (2025, 'Table 8.2 Public expenditure on cash old-age and survivor benefits'),
    'madn2r': (2025, 'Table 8.3 Private pension-benefit expenditures'),
    'vwcalj': (2025, 'Figure 9.1 Minimum or mandatory contribution rates in mandatory and auto-enrolment plans'),
    'f31vde': (2025, 'Table 9.1 Participation rate in pension plans'),
    '9ku103': (2025, 'Table 9.2 Assets earmarked for retirement'),
    'zx5sen': (2025, 'Figure 7.1 Income sources of older people'),
}

# key: (raw subdirectory, file name, url)
RAW = {
    **{'sl' + k: ('pag%d' % v[0], k + '.xlsx', 'https://stat.link/' + k) for k, v in STATLINK.items()},
    'pag2025pdf': ('pag2025', 'e40274c1-en.pdf', 'https://www.oecd.org/content/dam/oecd/en/publications/reports/'
                   '2025/11/pensions-at-a-glance-2025_76510fe4/e40274c1-en.pdf'),
    'pag2013pdf': ('pag2013', 'pension_glance-2013-en.pdf', 'https://www.oecd.org/content/dam/oecd/en/publications/'
                   'reports/2013/11/pensions-at-a-glance-2013_g1g3404e/pension_glance-2013-en.pdf'),
    # PaG 2013 table 4.19, Gini coefficients on pension entitlements and earnings (progressivity index);
    # the StatLink http://dx.doi.org/10.1787/888932907566 redirects here. Legacy .xls.
    'pag2013prog': ('pag2013', '812013221P1T014.xls', 'https://statlinks.oecdcode.org/812013221P1T014.xls'),
    'pag2015pdf': ('pag2015', 'pension_glance-2015-en.pdf','https://www.oecd.org/content/dam/oecd/en/publications/'
                   'reports/2015/12/pensions-at-a-glance-2015_g1g5a8b2/pension_glance-2015-en.pdf'),
    'prr2025':  ('pag2025', 'DF_PRR.csv', SDMX + 'OECD.ELS.SPD,DSD_PAG@DF_PRR,1.0/'
                 + _key(K26, 'A', '', '', '', '', '') + '?format=csvfilewithlabels'),
    'socx':     ('socx', 'DF_SOCX_AGG_oldAgeSurvivors.csv', SDMX + 'OECD.ELS.SPD,DSD_SOCX_AGG@DF_SOCX_AGG,1.0/'
                 + _key(K26, 'A', 'SOCX', 'PT_B1GQ+XDC', 'ES10+ES20+ES30+ES10_20+ES20_30', 'C+_T',
                        'TP01+TP11+TP111+TP112+TP113+TP21+TP211+TP212', '')
                 + '?startPeriod=%d&endPeriod=%d&format=csvfilewithlabels' % (Y0, Y1)),
    'socxDet':  ('socx', 'SOCX2025-All-nat-cur.xlsx',
                 'https://webfs.oecd.org/els-com/Social-Expenditure-Database/SOCX2025-All-nat-cur.xlsx'),
    'gpsBC':    ('gps', 'DF_BC.csv', SDMX + 'OECD.DAF.CM,DSD_FP@DF_BC,1.0/'
                 + _key(K26, 'A', '3100+3200+3300+4000', 'PT_B1GQ', '', '', '')
                 + '?startPeriod=%d&endPeriod=%d&format=csvfilewithlabels' % (Y0, Y1)),
    'gpsPA':    ('gps', 'DF_PA.csv', SDMX + 'OECD.DAF.CM,DSD_FP@DF_PA,1.0/'
                 + _key(K26, 'A', '1000', 'PT_B1GQ', '', '', '')
                 + '?startPeriod=%d&endPeriod=%d&format=csvfilewithlabels' % (Y0, Y1)),
    'na29':     ('na29', 'DF_TABLE29.csv', SDMX + 'OECD.SDD.NAD,DSD_NAPENS@DF_TABLE29,1.0/'
                 + _key('A', K26, 'S14', 'D61+D6111+D6121+D6131+D6141+D619+D62+D61SC_P+LE', '', '', '', '', '', '')
                 + '?startPeriod=%d&endPeriod=%d&format=csvfilewithlabels' % (Y0, Y1)),
    'esPens':   ('essp', 'spr_exp_pens.csv', ESDMX + 'spr_exp_pens/' + _key('A', '', 'TOTAL+MT+NMT', 'PC_GDP', KEU)
                 + '?format=SDMX-CSV&startPeriod=%d&endPeriod=%d' % (Y0, Y1)),
    'esFol':    ('essp', 'spr_exp_fol.csv', ESDMX + 'spr_exp_fol/' + _key('A', 'TOTAL', '', 'TOTAL+MT+NMT', 'PC_GDP', KEU)
                 + '?format=SDMX-CSV&startPeriod=%d&endPeriod=%d' % (Y0, Y1)),
    'esFsu':    ('essp', 'spr_exp_fsu.csv', ESDMX + 'spr_exp_fsu/' + _key('A', 'TOTAL', '', 'TOTAL+MT+NMT', 'PC_GDP', KEU)
                 + '?format=SDMX-CSV&startPeriod=%d&endPeriod=%d' % (Y0, Y1)),
    'esRec':    ('essp', 'spr_rec_sumt.csv', ESDMX + 'spr_rec_sumt/' + _key('A', '', 'PC_GDP+PC_REC', KEU)
                 + '?format=SDMX-CSV&startPeriod=%d&endPeriod=%d' % (Y0, Y1)),
    'esBen':    ('essp', 'spr_pns_ben.csv', ESDMX + 'spr_pns_ben/' + _key('A', '', 'TOTAL', '', 'T', 'PER', KEU)
                 + '?format=SDMX-CSV&startPeriod=%d&endPeriod=%d' % (Y0, Y1)),
}


def _path(key):
    d, f = RAW[key][:2]
    return os.path.join(RAWDIR, d, f)


def _retrieved(key):
    p = _path(key)
    return datetime.date.fromtimestamp(os.path.getmtime(p)).isoformat() if os.path.exists(p) else ''


# --- fetching (network; called from main only) -----------------------------------------------------------
def _get(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers = UA), timeout = 300).read()


def fetch(refetch = False, only = None):
    """ Download every raw file that is missing (all of them with refetch). A failed download is reported
    and skipped, so that one dead link does not stop the rest. """
    for key, (d, f, url) in RAW.items():
        if only and key not in only:
            continue
        p = _path(key)
        if os.path.exists(p) and not refetch:
            continue
        os.makedirs(os.path.dirname(p), exist_ok = True)
        print('  fetching', os.path.relpath(p, C.REPO), flush = True)
        try:
            data = _get(url)
        except Exception as e:                                   # noqa: BLE001 -- report, keep going
            print('    FAILED', type(e).__name__, e)
            continue
        with open(p, 'wb') as fh:
            fh.write(data)
    csvFromXls()


def csvFromXls():
    """ PaG 2013 table 4.19 is a legacy .xls (no xlrd in pyenv.md): Excel, through COM from PowerShell, saves
    its first sheet as csv into pag2013_csv/, outside the raw file's directory. Macros disabled (AutomationSecurity
    = 3), opened read-only. Skips when the csv exists. """
    import subprocess
    src = _path('pag2013prog')
    dst = os.path.join(RAWDIR, 'pag2013_csv', '812013221P1T014.csv')
    if not os.path.exists(src) or os.path.exists(dst):
        return
    os.makedirs(os.path.dirname(dst), exist_ok = True)
    ps = ("$xl = New-Object -ComObject Excel.Application; $xl.Visible = $false; $xl.DisplayAlerts = $false; "
          "$xl.AutomationSecurity = 3; try {{ $wb = $xl.Workbooks.Open('{}', 0, $true); "
          "$wb.Worksheets.Item(1).Activate(); $wb.SaveAs('{}', 62); $wb.Close($false) }} "
          "finally {{ $xl.Quit() }}").format(src, dst)
    subprocess.run(['powershell', '-NoProfile', '-Command', ps], check = True)


# --- reading ----------------------------------------------------------------------------------------------
COLS = ['iso3', 'item', 'tier', 'earnings_multiple', 'sex', 'value', 'unit', 'year', 'source', 'url',
        'retrieved', 'note']
EDITION_YEAR = {2021: 2020, 2025: 2024}        # PaG edition -> year of the pension rules it models
NZ_LOW = {2021: 'NZL low earner at 0.60 AW, not 0.5 (table footnote)',
          2025: 'NZL low earner at 0.63 AW, not 0.5 (StatLink footnote; the pdf footnote says 0.61)'}


def _grid(path, sheet = None):
    import openpyxl                              # read only; nothing in this module writes a workbook
    wb = openpyxl.load_workbook(path, read_only = True, data_only = True)
    ws = wb[sheet] if sheet else wb[wb.sheetnames[0]]
    rows = [list(r) for r in ws.iter_rows(values_only = True)]
    version = ''
    for sn in wb.sheetnames:
        for r in wb[sn].iter_rows(max_row = 12, values_only = True):
            for v in r:
                if isinstance(v, str) and v.strip().startswith('Version'):
                    version = v.strip()
    return rows, version


def _num(v):
    """ (float or None, flag): '(59.8)' -> 59.8; '7.0 [w]' -> 7.0, '[w]'; 'x', '..', '<5', ' ' -> None, text. """
    if v is None:
        return None, ''
    if isinstance(v, (int, float)):
        return float(v), ''
    s = str(v).strip()
    m = re.fullmatch(r'\(?\s*(-?\d+(?:\.\d+)?)\s*\)?\s*(\[[a-z, ]+\])?', s)
    if m:
        return float(m.group(1)), (m.group(2) or '')
    return None, s


def _iso(name):
    """ ISO3 of a PaG row label ('Austria*', 'Denmark (QMO)', 'Czechia'), with the stripped suffix. """
    if not isinstance(name, str):
        return None, ''
    s = name.strip()
    m = re.fullmatch(r'(.+?)\s*(\*+|\(.*\))?', s, re.S)
    if not m:
        return None, ''
    base, suffix = m.group(1).strip(), (m.group(2) or '')
    return NAME.get(base), suffix


class Rows(list):
    def add(self, iso3, item, tier, em, sex, value, unit, year, source, url, retrieved, note = ''):
        if iso3 in C26 and value is not None and not (isinstance(value, float) and np.isnan(value)):
            self.append(dict(zip(COLS, [iso3, item, tier, em, sex, value, unit, year, source, url, retrieved, note])))


def _sl(code):
    """ path, source label, url, retrieved, rules year, edition of a StatLink. """
    ed, title = STATLINK[code]
    key = 'sl' + code
    rows, version = _grid(_path(key))
    src = 'OECD, Pensions at a Glance {}, {} (StatLink {}{})'.format(ed, title, code,
                                                                     ', ' + version if version else '')
    return rows, src, RAW[key][2], _retrieved(key), EDITION_YEAR[ed], ed


def _panels(rows, n):
    """ Every (iso3, suffix, the n cells after the label) in a one- or two-panel PaG table. """
    out = []
    for r in rows:
        for j, v in enumerate(r):
            iso, suf = _iso(v)
            if iso:
                out.append((iso, suf, (list(r[j + 1:j + 1 + n]) + [None]*n)[:n]))
    return out


def pagReplacement(R):
    for code in ('b2f0ws', 'oeavq4', 'cofr31', 'ic8ung'):       # tables 4.1 (gross) and 4.4 (net), men/women
        rows, src, url, ret, yr, ed = _sl(code)
        item = 'grr_mandatory' if code in ('b2f0ws', 'oeavq4') else 'nrr_mandatory'
        for iso, suf, c in _panels(rows, 8):
            age, agew = _num(c[0])[0], _num(c[1])[0]
            if item == 'grr_mandatory':
                R.add(iso, 'pension_age', 'mandatory', '', 'M', age, 'years', yr, src, url, ret)
                R.add(iso, 'pension_age', 'mandatory', '', 'F', agew if agew is not None else age, 'years', yr,
                      src, url, ret, '' if agew is not None else
                      'same as men (the table shows women only where different)')
            for k, em in enumerate(('0.5', '1', '2')):
                m, w = _num(c[2 + 2*k])[0], _num(c[3 + 2*k])[0]
                nz = NZ_LOW[ed] if iso == 'NZL' and em == '0.5' else ''
                R.add(iso, item, 'mandatory', em, 'M', m, '% of individual gross earnings' if 'grr' in item
                      else '% of individual net earnings', yr, src, url, ret, nz)
                R.add(iso, item, 'mandatory', em, 'F', w if w is not None else m,
                      '% of individual gross earnings' if 'grr' in item else '% of individual net earnings', yr, src,
                      url, ret, '; '.join(x for x in (nz, '' if w is not None else
                                                       'same as men (the table shows women only where different)') if x))
    tiers = ('public_mandatory', 'private_mandatory', 'total_mandatory', 'voluntary', 'total_withVoluntary')
    for code in ('uwc8nt', 'pru4oj'):                          # table 4.2: gross by tier, men
        rows, src, url, ret, yr, ed = _sl(code)
        for iso, suf, c in _panels(rows, 15):
            for t, tier in enumerate(tiers):
                for k, em in enumerate(('0.5', '1', '2')):
                    v = _num(c[3*t + k])[0]
                    note = NZ_LOW[ed] if iso == 'NZL' and em == '0.5' else ''
                    if v is None and tier in ('private_mandatory', 'voluntary') and _num(c[k])[0] is not None:
                        note = '; '.join(x for x in (note, 'blank in the table: no such tier') if x)
                    R.add(iso, 'grr_tier', tier, em, 'M', v, '% of individual gross earnings', yr, src, url, ret, note)
    tiers = ('gross_mandatory', 'net_mandatory', 'gross_withVoluntary', 'net_withVoluntary')
    for code in ('scxdj3', 'f2v8nc'):                          # table 4.5
        rows, src, url, ret, yr, ed = _sl(code)
        for iso, suf, c in _panels(rows, 12):
            for t, tier in enumerate(tiers):
                for k, em in enumerate(('0.5', '1', '2')):
                    R.add(iso, 'rr_mandatoryVoluntary', tier, em, 'M', _num(c[3*t + k])[0],
                          '% of individual earnings', yr, src, url, ret,
                          NZ_LOW[ed] if iso == 'NZL' and em == '0.5' else '')
    # PaG 2025 replacement rates on the Data Explorer (the 2020 vintage is no longer served there)
    d = pd.read_csv(_path('prr2025'))
    src = 'OECD Data Explorer, OECD.ELS.SPD:DSD_PAG@DF_PRR(1.0) Pension replacement rates (Pensions at a Glance 2025)'
    for r in d.to_dict('records'):
        tier = ('gross_' if r['MEASURE'].startswith('G') else 'net_') + ('mandatory' if r['OPTIONALITY'] == 'M'
                                                                         else 'withVoluntary')
        R.add(r['REF_AREA'], 'prr_dataExplorer', tier, {'50': '0.5', '100': '1', '200': '2'}[r['MEASURE'][4:]],
              r['SEX'], float(r['OBS_VALUE']), r['Unit of measure'], int(r['TIME_PERIOD']), src, RAW['prr2025'][2],
              _retrieved('prr2025'), NZ_LOW[2025] if r['REF_AREA'] == 'NZL' and r['MEASURE'].endswith('R50') else '')


def pagFirstTier(R):
    kinds = ('residenceBasic', 'targeted', 'contributoryBasic', 'minimum')
    for code in ('uxgr8v', 'i28xu0'):                          # table 3.2
        rows, src, url, ret, yr, ed = _sl(code)
        for iso, suf, c in _panels(rows, 8):
            for k, kind in enumerate(kinds):
                R.add(iso, 'firstTier_value', kind, '', 'T', _num(c[k])[0], '% of gross average earnings', yr, src,
                      url, ret, 'benefit for new pensioners in %d' % yr)
                R.add(iso, 'firstTier_recipients', kind, '', 'T', _num(c[4 + k])[0], '% of population aged 65+',
                      yr, src, url, ret, 'can exceed 100: recipients include people under 65 or abroad')
    for code, tiers in (('kpoilu', ('residenceBasic', 'targeted_topUp', 'total_nonContributory')),
                        ('s5w1y7', ('residenceBasic', 'targeted_topUp', 'total_nonContributory')),
                        ('k23yri', ('contributoryBasic', 'minimum_topUp', 'total_contributory')),
                        ('7ohtnu', ('contributoryBasic', 'minimum_topUp', 'total_contributory'))):
        rows, src, url, ret, yr, ed = _sl(code)
        for iso, suf, c in _panels(rows, 3):
            for k, tier in enumerate(tiers):
                R.add(iso, 'firstTier_value_figure', tier, '', 'T', _num(c[k])[0], '% of gross average earnings',
                      yr, src, url, ret, 'figure: the top-up is what is paid on top of the basic amount'
                      + ('; subtitle says 2018 but values equal table 3.2 (2020)' if code == 'k23yri' else ''))
    cols = ('residenceBasic', 'targeted', 'contributoryBasic', 'minimum', 'secondTier_public', 'secondTier_private')
    for code in ('2xd6uo', 'mh9ruy'):                          # table 3.1, panel A (and B where given)
        rows, src, url, ret, yr, ed = _sl(code)
        panel = 'A'
        for r in rows:
            if any(isinstance(v, str) and v.startswith('Panel B') for v in r):
                panel = 'B'
            for j, v in enumerate(r):
                iso, suf = _iso(v)
                if not iso:
                    continue
                for k, col in enumerate(cols):
                    x = r[j + 1 + k] if j + 1 + k < len(r) else None
                    if isinstance(x, str) and x.strip():
                        tick = not re.search(r'[A-Za-z]', x)
                        R.add(iso, 'structure_panel' + panel, col, '', 'T', 1.0, 'indicator (1 = present)', yr,
                              src, url, ret, '' if tick else x.strip())


def pagContributions(R):
    tiers = ('employee_public', 'employer_public', 'employee_private', 'employer_private', 'total_nominal')
    for code in ('h4uwj7', 'orhufb', '9oqgj4'):                # tables 8.1 (and 8.2 of 2021: IRL, ESP, GBR)
        rows, src, url, ret, yr, ed = _sl(code)
        item = 'contribRate_socialInsurance' if code == '9oqgj4' else 'contribRate_mandatory'
        for iso, suf, c in _panels(rows, 7):
            note0 = {'*': 'rate also finances disability or invalidity benefits',
                     '**': 'see the table note (**)'}.get(suf, suf)
            for k, tier in enumerate(tiers):
                v, flag = _num(c[k])
                R.add(iso, item, tier, '', 'T', v, '% of gross earnings', yr, src, url, ret,
                      '; '.join(x for x in (note0, flag) if x))
            v, flag = _num(c[6])
            R.add(iso, item, 'effective_AW', '1', 'T', v, '% of gross earnings', yr, src, url, ret,
                  '; '.join(x for x in (note0, flag) if x))
            if c[5] is not None:
                v, flag = _num(c[5])
                unit = 'multiple of gross AW' if code != '9oqgj4' else '% of gross AW'
                if v is not None:
                    R.add(iso, 'contribCeiling', 'public/private', '', 'T', v, unit, yr, src, url, ret, flag)
                else:                                          # 'None', '1.08 / 8.62': kept as text
                    R.append(dict(zip(COLS, [iso, 'contribCeiling', 'public/private', '', 'T', np.nan, unit, yr,
                                             src, url, ret, 'published: ' + flag])))
    for code, sheet in (('8gtspo', None), ('vwcalj', 'Data')):  # figure 9.1: funded (quasi-)mandatory and AE plans
        rows, src, url, ret, yr, ed = _sl(code)
        rows = _grid(_path('sl' + code), sheet)[0] if sheet else rows
        for iso, suf, c in _panels(rows, 4):
            for k, tier in enumerate(('employees', 'employers', 'state', 'total_unsplit')):
                R.add(iso, 'contribRate_fundedPlans', tier, '1', 'T', _num(c[k])[0], '% of earnings', yr, src, url,
                      ret, suf)


def pagFunded(R):
    for code, item in (('gkvjah', 'coverage_pag'), ('f31vde', 'coverage_pag')):
        rows, src, url, ret, yr, ed = _sl(code)
        for iso, suf, c in _panels(rows, 5):
            for k, tier in enumerate(('mandatoryQuasiMandatory', 'autoEnrolment', 'voluntary_occupational',
                                      'voluntary_personal', 'voluntary_total')):
                v, flag = _num(c[k])
                if v is not None or flag not in ('', 'x', '..'):
                    R.append(dict(zip(COLS, [iso, item, tier, '', 'T', v if v is not None else np.nan,
                                             '% of population 15-64', 'latest (see table notes)', src, url, ret,
                                             'published: ' + flag if v is None else ''])))
    for code in ('qi2ev9', '9ku103'):
        rows, src, url, ret, yr, ed = _sl(code)
        for iso, suf, c in _panels(rows, 4):
            R.add(iso, 'assets_pag', 'retirementSavingsPlans' if ed == 2021 else 'pensionProviders', '', 'T',
                  _num(c[0])[0], '% of GDP', 'end-%d or latest' % (2020 if ed == 2021 else 2024), src, url, ret)
            R.add(iso, 'assets_pag', 'publicPensionReserveFunds', '', 'T', _num(c[2])[0], '% of GDP',
                  'end-%d or latest' % (2020 if ed == 2021 else 2024), src, url, ret)
    for code, latest in (('kp4d5j', 2017), ('madn2r', 'latest: 2021, 2022 (*) or 2023 (**)')):
        rows, src, url, ret, yr, ed = _sl(code)
        hdr = next(r for r in rows if any(v == 1990 for v in r))
        col = [i for i, v in enumerate(hdr) if v in (2017, 'Latest')][0]
        iso = None
        for r in rows:
            i, suf = _iso(r[0])
            if i:
                iso = i
            elif r[0] not in (None, '') or not isinstance(r[1], str):
                continue
            if isinstance(r[1], str) and r[1].strip() and iso:
                R.add(iso, 'privateSpending_pag', r[1].strip(), '', 'T', _num(r[col])[0], '% of GDP', latest, src,
                      RAW['sl' + code][2], ret, 'scheme type: m = mandatory, q = quasi-mandatory, v = voluntary'
                      + ('; ' + suf if suf else ''))
    for code in ('lqhcun', 'zx5sen'):                          # income sources of people over 65
        rows, src, url, ret, yr, ed = _sl(code)
        rows = _grid(_path('sl' + code), 'Data')[0]
        tiers = ('publicTransfers', 'occupationalTransfers', 'capital', 'work')
        if ed == 2021:
            for iso, suf, c in _panels(rows, 4):
                for k, tier in enumerate(tiers):
                    R.add(iso, 'incomeSources65', tier, '', 'T', _num(c[k])[0], '% of gross household income',
                          '2018 or latest', src, url, ret)
        else:                                                  # right block: name, iso3, year, %W, %K, %pub, %occ
            for r in rows:
                for j, v in enumerate(r[:-6]):
                    if isinstance(v, str) and v in C26 and isinstance(r[j + 1], (int, float)):
                        for k, tier in ((2, 'work'), (3, 'capital'), (4, 'publicTransfers'),
                                        (5, 'occupationalTransfers')):
                            R.add(v, 'incomeSources65', tier, '', 'T', _num(r[j + k])[0],
                                  '% of gross household income', int(r[j + 1]), src, url, ret)


def progressivity(R):
    p = os.path.join(RAWDIR, 'pag2013_csv', '812013221P1T014.csv')
    if not os.path.exists(p):
        print('  missing', os.path.relpath(p, C.REPO), '(convert the .xls with Excel, see csvFromXls)')
        return
    d = pd.read_csv(p, sep = ';', decimal = ',', header = None, skiprows = 11, encoding = 'utf-8-sig')
    src = ('OECD, Pensions at a Glance 2013, Table 4.19 Gini coefficients on pension entitlements and earnings '
           '(StatLink http://dx.doi.org/10.1787/888932907566, Version 1 16-Oct-2013)')
    for r in d.itertuples(index = False):
        iso = NAME.get(str(r[0]).strip())
        for k, tier in ((1, 'pensionGini_OECDdistribution'), (2, 'progressivityIndex_OECDdistribution'),
                        (3, 'pensionGini_nationalDistribution'), (4, 'progressivityIndex_nationalDistribution'),
                        (5, 'wageGini_national')):
            v = r[k] if pd.notna(r[k]) else None
            R.add(iso, 'progressivity2013', tier, '', 'T', float(v) if v is not None else None,
                  'index (100 = flat basic, 0 = proportional)' if 'Index' in tier else 'Gini x 100', 2012,
                  src, RAW['pag2013prog'][2], _retrieved('pag2013prog'),
                  '100 x (1 - pension Gini / earnings Gini), mandatory schemes only; pension rules of 2012')


SOCX_SRC = {'ES10': 'public', 'ES20': 'mandatoryPrivate', 'ES30': 'voluntaryPrivate'}
SOCX_YEARS = (2019, 2020)


def socx(R):
    d = pd.read_csv(_path('socx'))
    src = 'OECD Data Explorer, OECD.ELS.SPD:DSD_SOCX_AGG@DF_SOCX_AGG(1.0) Social expenditure aggregates'
    url, ret = RAW['socx'][2], _retrieved('socx')
    x = d[(d.UNIT_MEASURE == 'PT_B1GQ') & d.EXPEND_SOURCE.isin(SOCX_SRC) & d.TIME_PERIOD.isin(SOCX_YEARS)]
    for r in x.itertuples(index = False):
        if r.SPENDING_TYPE == '_T' and r.EXPEND_SOURCE != 'ES30':
            continue
        cash = r.SPENDING_TYPE == 'C'
        R.add(r.REF_AREA, 'socx_' + r.PROGRAMME_TYPE + ('_cash' if cash else '_cashInKind'),
              SOCX_SRC[r.EXPEND_SOURCE], '', 'T', round(float(r.OBS_VALUE), 4), '% of GDP', int(r.TIME_PERIOD), src,
              url, ret, '' if cash else 'voluntary private is published only as cash plus in kind')
    # programme-level public and mandatory-private old-age/survivors cash pensions, % of GDP with the GDP that
    # SOCX itself implies (XDC / PT_B1GQ of public old-age and survivors cash)
    agg = d[(d.EXPEND_SOURCE == 'ES10') & (d.SPENDING_TYPE == 'C') & (d.PROGRAMME_TYPE == 'TP01')]
    gdp = (agg[agg.UNIT_MEASURE == 'XDC'].query("PRICE_BASE == 'V'").set_index(['REF_AREA', 'TIME_PERIOD']).OBS_VALUE
           / agg[agg.UNIT_MEASURE == 'PT_B1GQ'].set_index(['REF_AREA', 'TIME_PERIOD']).OBS_VALUE * 100)
    w = pd.read_excel(_path('socxDet'), sheet_name = 'SOCX_NatCur', engine = 'openpyxl')
    cou = pd.read_excel(_path('socxDet'), sheet_name = 'SOCX_COU', header = None, engine = 'openpyxl')
    code2iso = {}
    for r in cou.iloc[4:, :3].dropna().itertuples(index = False):
        try:
            code2iso[int(r[2])] = r[1]
        except (TypeError, ValueError):
            pass
    w['iso3'] = w.COUNTRY.map(code2iso)
    w = w[w.iso3.isin(C26) & (w.TYPEXP == 1) & (w.NPROG > 0)
          & (((w.BRANCH == 1) & w.TYPROG.isin([1, 3])) | ((w.BRANCH == 2) & (w.TYPROG == 1)))]
    srcw = ('OECD SOCX, programme-level data, SOCX2025-All-nat-cur.xlsx (millions of national currency; % of GDP '
            'computed here with the GDP implied by DF_SOCX_AGG)')
    for r in w.to_dict('records'):
        for y in SOCX_YEARS:
            v = r.get(y)
            try:
                v = float(v)
            except (TypeError, ValueError):
                continue
            g = gdp.get((r['iso3'], y))
            if g is None or np.isnan(g) or v == 0:
                continue
            R.add(r['iso3'], 'socx_programme_' + ('survivorsPension' if r['BRANCH'] == 2 else
                                                  'oldAgePension' if r['TYPROG'] == 1 else 'oldAgeOtherCash'),
                  {10: 'public', 20: 'mandatoryPrivate', 30: 'voluntaryPrivate'}[r['SOURCE']] + ':'
                  + programmeClass(r['iso3'], str(r['Title_Emm'])), '', 'T', round(100*v/g, 4), '% of GDP', y, srcw,
                  RAW['socxDet'][2], _retrieved('socxDet'), '{} [{}]'.format(str(r['Title_Emm']).strip(), r['Code']))


def programmeClass(iso, title):
    """ This collection's reading of a SOCX programme title (not a SOCX field): govEmployees, meansTested,
    targetedPensionTested (PaG's 'targeted' where SOCX says non-means-tested), basic, earningsRelated, other. """
    t = title.lower()
    if (iso == 'FIN' and '(kel)' in t) or (iso == 'SWE' and 'tax financed' in t):
        return 'targetedPensionTested'        # national pension / guarantee pension, reduced by own pension
    if re.search(r'civil serv|public sector occup|public employee|government employ|public service mutual|'
                 r'local government employees|military|soldiers|\babp\b|kvtel|state employee|other public sector|'
                 r'magistrates|public servants|central government \(non means|administrations publiques centrales|'
                 r'public sector employers', t):
        return 'govEmployees'
    if re.search(r'veteran|war |victims|reparation|injury|accident|resistance|lag\)|calimi|care|seamen|mining|'
                 r'night shift', t):
        return 'other'
    if ('non means' not in t and 'no means' not in t and 'sans conditions' not in t
            and re.search(r'means.tested|sous conditions de ressources|income support|supplemental security|'
                          r'guaranteed income|minimum guaranteed|low income|'
                          r'complementary benefits|additional income support|social welfare|'
                          r'social assistance|^allowance$|^age pension|dept of social protection|'
                          r'solidarity', t)):
        return 'meansTested'
    if re.search(r'\baow\b|new zealand superannuation|old age security|basic pension|flat rate|^national pension$|'
                 r'national basic', t):
        return 'basic'
    return 'earningsRelated'


def essp(R):
    url0 = 'https://ec.europa.eu/eurostat/databrowser/view/{}/default/table'
    for key, ds, dim, keep, item in (('esPens', 'spr_exp_pens', 'spdepb', None, 'essp_pensions_'),
                                     ('esFol', 'spr_exp_fol', 'spdep', ('CASH', 'CASH_P', 'CASH_P_OLD_PEN',
                                                                        'CASH_P_AOLD_PEN', 'CASH_P_PART_PEN'), 'essp_oldAge_'),
                                     ('esFsu', 'spr_exp_fsu', 'spdep', ('CASH', 'CASH_P_SRV_PEN'), 'essp_survivors_')):
        d = pd.read_csv(_path(key))
        d = d[d.TIME_PERIOD.isin(SOCX_YEARS)]
        src = 'Eurostat ESSPROS, {} (all schemes), last update {}'.format(ds, d['LAST UPDATE'].iloc[0])
        for r in d.itertuples(index = False):
            if keep and getattr(r, dim) not in keep:
                continue
            R.add(ESTAT.get(r.geo), item + getattr(r, dim), {'MT': 'meansTested', 'NMT': 'nonMeansTested',
                                                             'TOTAL': 'total'}[r.spdepm], '', 'T',
                  float(r.OBS_VALUE), '% of GDP', int(r.TIME_PERIOD), src, RAW[key][2], _retrieved(key),
                  'flag ' + r.OBS_FLAG if isinstance(r.OBS_FLAG, str) else '')
    d = pd.read_csv(_path('esRec'))
    d = d[d.TIME_PERIOD.isin(SOCX_YEARS)]
    src = 'Eurostat ESSPROS, spr_rec_sumt (all schemes, all functions), last update {}'.format(d['LAST UPDATE'].iloc[0])
    for r in d.itertuples(index = False):
        R.add(ESTAT.get(r.geo), 'essp_receipts', r.sptype, '', 'T', float(r.OBS_VALUE),
              {'PC_REC': '% of total receipts', 'PC_GDP': '% of GDP'}[r.unit], int(r.TIME_PERIOD), src,
              RAW['esRec'][2], _retrieved('esRec'), 'all social protection, not pensions alone')
    d = pd.read_csv(_path('esBen'))
    d = d[d.TIME_PERIOD.isin(SOCX_YEARS) & d.spdepb.isin(['OLD_TOT', 'OLD', 'SRV', 'OLD_SRV'])]
    src = 'Eurostat ESSPROS, spr_pns_ben (all schemes), last update {}'.format(d['LAST UPDATE'].iloc[0])
    for r in d.itertuples(index = False):
        R.add(ESTAT.get(r.geo), 'essp_beneficiaries_' + r.spdepb, {'MT': 'meansTested', 'NMT': 'nonMeansTested',
                                                                   'TOTAL': 'total'}[r.spdepm], '', 'T',
              float(r.OBS_VALUE), 'persons', int(r.TIME_PERIOD), src, RAW['esBen'][2], _retrieved('esBen'),
              'a person with both kinds counts in each')


def gps(R):
    d = pd.read_csv(_path('gpsBC'))
    d = d[d.TIME_PERIOD.isin(SOCX_YEARS)]
    src = 'OECD Data Explorer, OECD.DAF.CM:DSD_FP@DF_BC(1.0) Benefit payments and contributions to asset-backed pensions'
    for (iso, m, y), g in d.groupby(['REF_AREA', 'MEASURE', 'TIME_PERIOD']):
        cells = {(r.PLAN_TYPE, r.DEFINITION_TYPE, r.VEHICLE_TYPE): r.OBS_VALUE for r in g.itertuples(index = False)}
        pick = ([('_T', '_T', '_T')] if ('_T', '_T', '_T') in cells else
                [k for k in (('OCC', '_T', '_T'), ('PER', '_T', '_T')) if k in cells] or
                [k for k in (('_T', '_T', 'PF'),) if k in cells])
        for k in pick:
            R.add(iso, {3100: 'gps_contributions', 4000: 'gps_benefits'}[int(m)], '/'.join(k), '', 'T',
                  round(float(cells[k]), 4), '% of GDP', int(y), src, RAW['gpsBC'][2], _retrieved('gpsBC'),
                  'plan/definition/vehicle; the broadest cell published' + ('' if k == ('_T', '_T', '_T') else
                                                                            ' (no all-plan total)'))
    d = pd.read_csv(_path('gpsPA'))
    d = d[d.TIME_PERIOD.isin(SOCX_YEARS) & (d.VEHICLE_TYPE == '_T')]
    src = 'OECD Data Explorer, OECD.DAF.CM:DSD_FP@DF_PA(1.0) Assets earmarked for retirement'
    for r in d.itertuples(index = False):
        R.add(r.REF_AREA, 'gps_assets', 'allVehicles', '', 'T', round(float(r.OBS_VALUE), 3), '% of GDP',
              int(r.TIME_PERIOD), src, RAW['gpsPA'][2], _retrieved('gpsPA'))


def na29(R):
    d = pd.read_csv(_path('na29'))
    d = d[d.TIME_PERIOD.between(2018, 2021) & d.PENSION_FUNDTYPE.isin(['S13PS', 'S13PBI', 'S13PBX', 'S12PBI', 'S13PC',
                                                                      'S12P'])
          & d.TRANSACTION.isin(['D6111', 'D6121', 'D6131', 'D6141', 'D619', 'D62'])]
    src = 'OECD Data Explorer, OECD.SDD.NAD:DSD_NAPENS@DF_TABLE29(1.0) Social insurance pensions (national accounts table 29)'
    for r in d.itertuples(index = False):
        R.add(r.REF_AREA, 'na29_' + r.TRANSACTION, r.PENSION_FUNDTYPE, '', 'T', float(r.OBS_VALUE),
              'millions of national currency', int(r.TIME_PERIOD), src, RAW['na29'][2], _retrieved('na29'),
              'S13PS social security; S13PBI/S13PBX/S12PBI government-managed DB for its employees; S12P not managed '
              'by government; social assistance (tax-financed basic or means-tested pensions) is outside the table')


def collect():
    R = Rows()
    for f in (pagReplacement, pagFirstTier, pagContributions, pagFunded, progressivity, socx, essp, gps, na29):
        n = len(R)
        f(R)
        print('  {:18s} {:6d} rows'.format(f.__name__, len(R) - n), flush = True)
    df = pd.DataFrame(R, columns = COLS)
    df.insert(0, 'status', 'DRAFT: not wired into anything')
    return df


def main():
    p = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    p.add_argument('--force', action = 'store_true', help = 'rebuild the csv even if it exists')
    p.add_argument('--refetch', action = 'store_true', help = 'download the raw files again')
    p.add_argument('--fetch-only', action = 'store_true', help = 'download missing raw files and stop')
    args = p.parse_args()
    if os.path.exists(OUT) and not (args.force or args.refetch or args.fetch_only):
        print('up to date: {} (--force to rebuild)'.format(os.path.relpath(OUT, C.REPO)))
        return
    fetch(args.refetch)
    if args.fetch_only:
        return
    df = collect()
    df.to_csv(OUT, index = False)
    print('{}: {} rows, {} items'.format(os.path.relpath(OUT, C.REPO), len(df), df['item'].nunique()))


if __name__ == '__main__':
    main()
