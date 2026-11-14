r""" Stage (0), DRAFT: spending and recipients of the minimum-pension top-ups and tested old-age benefits in the
ten countries whose share s was missing (ITA FRA AUT DEU BEL PRT NOR POL LUX HUN), 2019 and 2020. Not wired
into build.py, REPLICATION.md, the tests or any calibration; notes/data_thetaTau.md section G reads it.

Run:  .venv\Scripts\python.exe -I python\paper\oecd26Topups.py [--force] [--refetch]
      -> data/oecd26_topups.csv (one row per scheme x year; first line a `#` comment, read with comment="#"),
         and prints s = top-up spending / SOCX public old-age and survivors cash (TP01, national currency), 2020

`main` is the only code here that touches the network: `fetch` downloads every raw file missing from
logs/thetaTau_raw/topups/<ISO3>/ (gitignored; --refetch downloads them all again) and writes pdftotext's
-layout and -raw text (UTF-8) beside each pdf, as <name>.txt and <name>_raw.txt. Everything else reads files:
the scheme values were read by hand off those raw files and are typed in ROWS with the file they come from and
a quote, and `check` raises unless every quote appears in that file's text (whitespace collapsed); SOCX programme values are read from the workbook in
logs/thetaTau_raw/harmonised/socx/ by programme code; GDP is Eurostat nama_10_gdp (CP_MNAC); PaG's 0.5 AW
pension and first-tier levels come from data/oecd26_pensionTiers_harmonised.csv. main skips writing when the
csv exists unless --force, and prints s either way. Raw files are untrusted data: never executed.

Units: spending_nc in millions of national currency; spending_pctGDP against Eurostat GDP of the same year.
reachedAtHalfAW compares PaG 2021's public pension at 0.5 AW (table 4.2, men, % of AW) with the benefit's
2020 level in % of AW (PaG 2021 table 3.2 unless the note says otherwise). No interpolation: a year without a
published value has no row, and a row's `year` is the year the value is for.
"""
import os, re, sys, json, html, shutil, argparse, subprocess, urllib.request
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

OUT = os.path.join(C.DATA, 'oecd26_topups.csv')
RAWDIR = os.path.join(C.REPO, 'logs', 'thetaTau_raw', 'topups')
SOCXDIR = os.path.join(C.REPO, 'logs', 'thetaTau_raw', 'harmonised', 'socx')
TIERS = os.path.join(C.DATA, 'oecd26_pensionTiers_harmonised.csv')
RET = '2026-10-30'
ISO = dict(ITA='IT', FRA='FR', AUT='AT', DEU='DE', BEL='BE', PRT='PT', NOR='NO', POL='PL', LUX='LU', HUN='HU')
CUR = dict(ITA='EUR', FRA='EUR', AUT='EUR', DEU='EUR', BEL='EUR', PRT='EUR', LUX='EUR', NOR='NOK', POL='PLN',
           HUN='HUF')
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) '
                    'Chrome/124.0 Safari/537.36'}
WB = 'https://web.archive.org/web/{}id_/{}'
ESGDP = ('https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/nama_10_gdp?format=JSON&unit=CP_MNAC'
         '&na_item=B1GQ&time=2019&time=2020&time=2021&time=2023' + ''.join('&geo=' + g for g in ISO.values()))

# key -> (file under RAWDIR, url it was downloaded from, short source name)
RAW = {
    'itaIX': ('ITA/IP_IXRapporto_2020data.pdf',
              'https://www.itinerariprevidenziali.it/wp-content/uploads/2025/11/32056669IXRapporto-IlBil.pdf',
              'Itinerari Previdenziali, Il Bilancio del Sistema Previdenziale italiano, Rapporto n. 9 (2022, data 2020)'),
    'itaVIII': ('ITA/IP_VIIIRapporto_2019data.pdf', 'https://www.money.it/IMG/pdf/32055154ottavorapportosu_2_.pdf',
                'Itinerari Previdenziali, Rapporto n. 8 (2021, data 2019)'),
    'itaX': ('ITA/IP_XRapporto_2021data.pdf',
             'https://www.bollettinoadapt.it/wp-content/uploads/2023/01/32058292XRapportosuIlBil_.pdf',
             'Itinerari Previdenziali, Rapporto n. 10 (2023, data 2021)'),
    'fraDD': ('FRA/DREES_DD_dispositifsSolidarite_fin2020.pdf',
              WB.format('20250410115046', 'https://drees.solidarites-sante.gouv.fr/sites/default/files/2025-04/'
                        'DDSolidarit%C3%A9_2020_MEL%20(1).pdf'),
              'DREES, Les Dossiers de la DREES 130 (April 2025), Fin 2020, les dispositifs de solidarite representent '
              '20 % des depenses de retraite (EIR 2020)'),
    'fraMS': ('FRA/DREES_Minima22_fiche07_depenses.pdf',
              WB.format('20230927110609', 'https://drees.solidarites-sante.gouv.fr/sites/default/files/2022-09/'
                        'Minima22%20Fiche%207%20-%20Les%20d%C3%A9penses%20d%E2%80%99allocations%20des%20minima'
                        '%20sociaux.pdf'),
              'DREES, Minima sociaux et prestations sociales, edition 2022, fiche 07, tableau 1'),
    'autBVA21': ('AUT/BMF_UG22_Teilheft_2021.pdf',
                 'https://service.bmf.gv.at/Budget/Budgets/2021/bfg/teilhefte/UG22/UG22_Teilheft_2021.pdf',
                 'BMF, Bundesvoranschlag 2021, Teilheft UG 22, Detailbudget 22.01.02 Ausgleichszulagen variabel, Erfolg 2019'),
    'autBVA22': ('AUT/BMF_UG22_Teilheft_2022.pdf',
                 'https://service.bmf.gv.at/Budget/Budgets/2022/bfg/teilhefte/UG22/UG22_Teilheft_2022.pdf',
                 'BMF, Bundesvoranschlag 2022, Teilheft UG 22, Detailbudget 22.01.02 Ausgleichszulagen variabel, Erfolg 2020'),
    'autORF': ('AUT/ORF_wien_3116351.html', 'https://wien.orf.at/radio/stories/3116351/',
               'ORF Wien, Alles zur Ausgleichszulage bei der Pension (press)'),
    'deuBT': ('DEU/BT-Drs_20-12534.pdf', 'https://dserver.bundestag.de/btd/20/125/2012534.pdf',
              'Deutscher Bundestag, Drucksache 20/12534 (Antwort der Bundesregierung, 09.08.2024), zu Frage 11'),
    'deuDS': ('DEU/destatis_PD21_188_228.html',
              'https://www.destatis.de/DE/Presse/Pressemitteilungen/2021/04/PD21_188_228.html',
              'Destatis, Pressemitteilung 188 of 15 April 2021'),
    'deuKOM': ('DEU/kommunal_grundsicherung-im-alter-2026.html',
               'https://kommunal.de/grundsicherung-im-alter-2026-rekord-auslaender-kosten',
               'kommunal.de (2026), citing a government answer reported by BILD (press)'),
    'belRH': ('BEL/Kamer_DOC55-0416-002_Rekenhof_IGO.pdf', 'https://www.dekamer.be/FLWB/PDF/55/0416/55K0416002.pdf',
              'Belgische Kamer, DOC 55 0416/002 (12 March 2020), advies van het Rekenhof'),
    'prtCFP': ('PRT/CFP_rel05-2021_SS_CGA_2020.pdf', 'https://www.cfp.pt/uploads/publicacoes_ficheiros/cfp-rel-05-2021.pdf',
               'CFP, Evolucao orcamental da Seguranca Social e da CGA em 2020, Relatorio 05/2021, Quadro 4 (IGFSS)'),
    'prtSS': ('PRT/SS_em_numeros_2023.pdf',
              'https://sisscontent.seg-social.pt/documents/10152/113014/Seguranca_Social_em_numeros_2023_PT.pdf/'
              '628ae617-32d8-426d-ac5c-a3d39be1932d',
              'Seguranca Social em Numeros, June 2023 (II, Estatisticas da Seguranca Social)'),
    'norNAV': ('NOR/NAV_statistikknotat_alderspensjon_des2020.pdf',
               'https://www.nav.no/_/attachment/download/ba80b626-55c1-4808-92fe-462141df3a8a:f0b2f182a09509a5b5b352c46'
               'a41e3edcb1b3ce2/Statistikknotat%20per%20desember%202020%20-%20alderspensjon%20v101.pdf',
               'NAV, Statistikknotat alderspensjon per desember 2020'),
    'polZUS': ('POL/ZUS_Struktura_marzec2020.pdf',
               'https://www.zus.pl/documents/10182/39637/Struktura+wysoko%C5%9Bci+%C5%9Bwiadcze%C5%84+wyp%C5%82acanych'
               '+przez+ZUS+po+waloryzacji+w+marcu+2020+r.pdf/3df99481-587f-8efe-8f1f-a77871e4f8c8',
               'ZUS, Struktura wysokosci swiadczen wyplacanych przez ZUS po waloryzacji w marcu 2020, tabela 2'),
    'polPRAWO': ('POL/prawo_518354.html',
                 'https://www.prawo.pl/kadry/emerytury-nizsze-niz-minimalna-emerytura-ile-i-komu-budzet,518354.html',
                 'prawo.pl, data provided by ZUS (press)'),
    'luxPJ': ('LUX/paperjam_complement-pension-minimum.html', 'https://paperjam.lu/article/complement-pension-minimum-bou',
              'Paperjam, Le complement pour pension minimum, citing the ministry (press)'),
    'hunPOP': ('HUN/eurostat_demo_pjangroup_65plus.json',
               'https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/demo_pjangroup?format=JSON&sex=T'
               '&age=Y65-69&age=Y70-74&age=Y75-79&age=Y80-84&age=Y_GE85&geo=HU&time=2020&time=2021',
               'Eurostat demo_pjangroup, population on 1 January, 65+'),
    'gdp': ('gdp/eurostat_nama_10_gdp_CPMNAC.json', ESGDP, 'Eurostat nama_10_gdp, B1GQ, CP_MNAC'),
}
SOCX_XLSX = os.path.join(SOCXDIR, 'SOCX2025-All-nat-cur.xlsx')
SOCX_AGG = os.path.join(SOCXDIR, 'DF_SOCX_AGG_oldAgeSurvivors.csv')
SOCX_URL = 'https://webfs.oecd.org/els-com/Social-Expenditure-Database/SOCX2025-All-nat-cur.xlsx'
PAG = 'PaG 2021 table 3.2 (StatLink uxgr8v)'

# s roles: 'main' enters the country's s, 'alt' an alternative numerator printed beside it, None neither
# (a component of another row, reached at 0.5 AW, or a check). thr: ('pag', tier) reads PaG table 3.2's level of
# that first tier (% AW, 2020), a float is a level in % AW given in the note, None means the note decides `reach`.
ROWS = [
    # ---------------- ITA (INPS archive and Casellario, via Itinerari Previdenziali)
    dict(iso3='ITA', scheme='Integrazione al trattamento minimo (IVS pensions)', year=2019, nc=7470, rec=2778509,
         base='integrated IVS pensions in payment, 31 December', raw='itaIX', conf='medium',
         q=['Integrazioni al minimo 2.778.509 2.648.653 7.470 7.024'], thr=('pag', 'minimum'), s=None,
         socx='part | inside the INPS fund pensions (old age [380.10.1.1.1.6], [380.10.1.1.1.11] and others; survivors); '
              'the part on invalidity pensions is SOCX incapacity',
         note='IX report table 5.6. The VIII report (table 6.6) gives 7,470.1 and 2,778,509; ACLI gives 8,830 and '
              '3,181,525 for 2016. IVS = old-age, survivors and invalidity pensions; no split by category found, so an '
              'upper bound of the old-age and survivors part. Not available to NDC-only careers.'),
    dict(iso3='ITA', scheme='Integrazione al trattamento minimo (IVS pensions)', year=2020, nc=7024, rec=2648653,
         base='integrated IVS pensions in payment, 31 December', raw='itaIX', conf='medium',
         q=['Integrazioni al minimo 2.778.509 2.648.653 7.470 7.024', '2.1) Integrazioni al minimo 7.024 0,42'],
         thr=('pag', 'minimum'), s='main',
         socx='part | inside the INPS fund pensions (old age [380.10.1.1.1.6], [380.10.1.1.1.11] and others; survivors); '
              'the part on invalidity pensions is SOCX incapacity',
         note='IX report table 5.6 and box 2 (0.42% of GDP there). The X report (table 5.6) revises 2020 to 7,030.1 and '
              '2,648,508. IVS includes invalidity pensions: an upper bound of the old-age and survivors part.'),
    dict(iso3='ITA', scheme='Maggiorazioni sociali (IVS pensions)', year=2019, nc=1409, rec=851317,
         base='pensions with a maggiorazione, 31 December', raw='itaIX', conf='medium',
         q=['Maggiorazioni sociali 851.317 1.143.670 1.409 2.633'], thr=None, reach='no', s=None,
         q2=('itaVIII', 'le maggiorazioni sociali (1,409 miliardi)'),
         socx='part | paid with the pension, financed by GIAS ([380.10.1.1.1.10] old age, GIAS)',
         note='IX report table 5.6 (the text under the VIII report\'s table 8.4 also gives 1,409). Income-tested supplements (L. 140/1985, '
              'L. 544/1988, art. 38 L. 448/2001): the guaranteed total lies below PaG\'s 0.5 AW pension (37.3% AW), '
              'so not reached; limits not collected here.'),
    dict(iso3='ITA', scheme='Maggiorazioni sociali on IVS pensions (L. 140/1985 and trattamenti minimi)', year=2020,
         nc=343 + 1162, rec=None, base='', raw='itaIX', conf='medium',
         q=['2.3) Maggiorazioni sociali L.n. 140/1985 343 0,02', 'Trattamenti minimi (Leggi varie) 1.162 0,07'],
         thr=None, reach='no', s='main',
         socx='part | paid with the pension, financed by GIAS ([380.10.1.1.1.10] old age, GIAS)',
         note='IX report box 2, items 2.3 (343) + 2.4 (1,162). Its table 5.6 gives 2,633 for 2020 (1,143,670 pensions), '
              'which also counts the maggiorazioni on civil-invalidity pensions (box 2 item 1.2: 968) and is not '
              'comparable with 2019; that count is the only 2020 recipient figure. Not reached (see 2019 row).'),
    dict(iso3='ITA', scheme='Pensioni e assegni sociali', year=2019, nc=4661, rec=796727,
         base='social pensions and allowances in payment, 31 December', raw='itaIX', conf='medium',
         q=['Pensioni e assegni sociali 796.727 803.441 4.661 4.803'], thr=('pag', 'targeted'), s=None,
         socx='yes | Old age pension (means-tested) Social Assistance Funds ... [380.10.1.1.1.15]',
         note='IX report table 5.6, without the maggiorazioni paid on them; with them, the GIAS outlay was 4,720 '
              '(VIII report, note 46).'),
    dict(iso3='ITA', scheme='Pensioni e assegni sociali, with their maggiorazioni', year=2020, nc=4849, rec=803441,
         base='social pensions and allowances in payment, 31 December (table 5.6)', raw='itaIX', conf='medium',
         q=['maggiorazioni 4.849 0,29', 'Pensioni e assegni sociali 796.727 803.441 4.661 4.803'],
         thr=('pag', 'targeted'), s='main',
         socx='yes | Old age pension (means-tested) Social Assistance Funds ... [380.10.1.1.1.15]',
         note='IX report box 2 item 1.4 (0.29% of GDP); table 5.6 gives 4,803 without the maggiorazioni. The SOCX '
              'line is 0.29% of GDP too, so the assegno sociale is the 0.30% found earlier, to rounding.'),
    dict(iso3='ITA', scheme='Quattordicesima (somma aggiuntiva, L. 127/2007)', year=2019, nc=1478, rec=3044619,
         base='pensioners receiving it', raw='itaIX', conf='medium',
         q=['Quattordicesima 3.044.619 2.893.837 1.478 1.403'], thr=2 * 22.7, s=None,
         socx='part | paid with the pension', note='IX report table 5.6. Income limit 2 x trattamento minimo '
         '(2 x 22.7% AW, PaG table 3.2): PaG\'s 0.5 AW pension is below it, so it is reached and stays out of s.'),
    dict(iso3='ITA', scheme='Quattordicesima (somma aggiuntiva, L. 127/2007)', year=2020, nc=1403, rec=2893837,
         base='pensioners receiving it', raw='itaIX', conf='medium',
         q=['Quattordicesima 3.044.619 2.893.837 1.478 1.403', '127/2007) 1.467 0,09'], thr=2 * 22.7, s=None,
         socx='part | paid with the pension', note='IX report table 5.6 (box 2 item 2.2 gives 1,467). Reached at '
         '0.5 AW (income limit 2 x trattamento minimo), so not in s.'),
    dict(iso3='ITA', scheme='Importo aggiuntivo (art. 70 L. 388/2000)', year=2019, nc=50, rec=330703,
         base='pensioners receiving it', raw='itaIX', conf='medium',
         q=['Importo aggiuntivo 330.703 541.535 50 82'], thr=1.5 * 22.7, s=None,
         socx='part | paid with the pension', note='IX report table 5.6. Income limit 1.5 x trattamento minimo.'),
    dict(iso3='ITA', scheme='Importo aggiuntivo (art. 70 L. 388/2000)', year=2020, nc=82, rec=541535,
         base='pensioners receiving it', raw='itaIX', conf='medium',
         q=['Importo aggiuntivo 330.703 541.535 50 82'], thr=1.5 * 22.7, s='main',
         socx='part | paid with the pension', note='IX report table 5.6. Income limit 1.5 x trattamento minimo.'),
    # ---------------- FRA (DREES)
    dict(iso3='FRA', scheme='Minima de pension, all schemes (MiCo, Miga, PMR, CD), direct-right pensions', year=2020,
         nc=8530, rec=6005400, base='beneficiaries at 31 December 2020', raw='fraDD', conf='high',
         q=['Minima de pension 8,53 Mᵈ€ 2,9 % 1 420 € 12,6 % 6 005 400'], thr=22.1, s='main',
         socx='yes | inside the old-age pension programmes (regime general [250.10.1.1.1.2], autres regimes '
              '[250.10.1.1.1.3], administrations publiques centrales [250.10.1.1.1.4]); no separate line',
         note='Dossier 130, tableau 3 (EIR 2020): 2.9% of direct-right pensions. Threshold: MiCo majore EUR 702.54 a '
              'month (tableau 1) x 12 / PaG AW 38,188 = 22.1% AW (PaG table 3.2: 20.2, the base MiCo). PaG\'s 0.5 AW '
              'pension (30.1% AW) is all tiers; MiCo tops up the CNAV pension alone and is capped at a total pension of '
              'EUR 1,191.56 a month, above the 0.5 AW total (EUR 958), so reach turns on the CNAV part, which PaG does '
              'not publish; PaG\'s replacement rates are equal at 0.5 and 1 AW, so MiCo does not bind in PaG.'),
    dict(iso3='FRA', scheme='of which: CNAV minimum contributif (component of the row above, not additive)', year=2020,
         nc=6100, rec=None, base='', raw='fraDD', conf='high', q=['Minima de pension 6,10 Mᵈ€ 5,5 %'], thr=22.1, s=None,
         socx='yes | Pension de vieillesse sans conditions de ressources - Regime general [250.10.1.1.1.2]',
         note='Dossier 130, tableau 7 (5.5% of CNAV direct-right pensions). Private-sector aligned base schemes '
              '(CNAV, MSA salaries, RSI): 5,567,100 beneficiaries, EUR 1,230 a year on average (tableau 5).'),
    dict(iso3='FRA', scheme='of which: civil-service minimum garanti (SRE civil, CNRACL, SRE militaire; component, not '
         'additive)', year=2020, nc=310 + 1230 + 50, rec=596100, base='fonction publique base-scheme beneficiaries',
         raw='fraDD', conf='high', q=['Minima de pension 6,10 Mᵈ€ 5,5 % - - 0,31 Mᵈ€ 0,7% 1,23 Mᵈ€ 6,0% 0,05 Mᵈ€ 0,6%',
                                       'Minima de pension 2,3% 2 650 € 596 100'], thr=None, reach='no', s=None,
         socx='yes | government-employee programmes ([250.10.1.1.1.4] and CNRACL inside [250.10.1.1.1.3])',
         note='Dossier 130, tableau 7 (0.31 + 1.23 + 0.05 bn) and tableau 5 (596,100). Civil servants: outside PaG\'s '
              'private-sector worker.'),
    dict(iso3='FRA', scheme='Minimum vieillesse (ASPA and ASV), allocations', year=2020, nc=3530, rec=635300,
         base='recipients at end 2020', raw='fraMS', conf='high',
         q=['(ASV et Aspa)3 Minimum vieillesse 635 300 14,2 3 530 +11,7 11,6 476'], thr=('pag', 'targeted'), s='main',
         socx='yes | Pension de vieillesse sous conditions de ressources - Regime general [250.10.1.1.3.2] and autres '
              'regimes [250.10.1.1.3.3]',
         note='Excludes the ASV first-floor allocations. EIR 2020 (Dossier 130, tableau 3) gives 3,810 and 717,800 '
              'beneficiaries. 2019: +11.7% real growth 2019-20 (no 2019 level collected).'),
    # ---------------- AUT
    dict(iso3='AUT', scheme='Ausgleichszulage, federal reimbursement to the pension insurers (all pension types)',
         year=2019, nc=979.82575110, rec=None, base='', raw='autBVA21', conf='high',
         q=['979,825.751,10'], thr=('pag', 'minimum'), s=None,
         socx='part | old age [40.10.1.1.1.2] and survivors [40.10.2.1.1.2]; the part on invalidity pensions is SOCX '
              'incapacity [40.10.3.1.1.2]',
         note='Erfolg 2019, Ergebnishaushalt (Finanzierungshaushalt 979,005,121.46). The federal budget reimburses the '
              'whole AZ; the three SOCX lines sum to 976.8. Recipients: 205,306 in December 2019 (BMSGPK evaluation, '
              'secondary host; data/oecd26_pensionSchemes_A.csv).'),
    dict(iso3='AUT', scheme='Ausgleichszulage, federal reimbursement to the pension insurers (all pension types)',
         year=2020, nc=1052.19792545, rec=198378, base='recipients, December 2020 (press)', raw='autBVA22',
         conf='high (spending); low (recipients)', q=['1.052,197.925,45'], thr=('pag', 'minimum'), s='alt',
         socx='part | old age [40.10.1.1.1.2] and survivors [40.10.2.1.1.2]; the part on invalidity pensions is SOCX '
              'incapacity [40.10.3.1.1.2]',
         note='Erfolg 2020, Ergebnishaushalt; the three SOCX lines sum to 1,045.2. Recipients from ORF Wien '
              '(198,378; 67.3% women), not checked against the Dachverband.', q2=('autORF', '(198.378)')),
    dict(iso3='AUT', scheme='Ausgleichszulage on old-age and survivors pensions (SOCX)', year=2019,
         socxcodes=['40.10.1.1.1.2', '40.10.2.1.1.2'], rec=None, base='', conf='high', thr=('pag', 'minimum'),
         s=None, socx='yes | [40.10.1.1.1.2] Supplement for those on low income: Statutory pension insurance + '
                      '[40.10.2.1.1.2] Widow(er)s pension (means-tested)',
         note='SOCX programme values; incapacity adds [40.10.3.1.1.2].'),
    dict(iso3='AUT', scheme='Ausgleichszulage on old-age and survivors pensions (SOCX)', year=2020,
         socxcodes=['40.10.1.1.1.2', '40.10.2.1.1.2'], rec=None, base='', conf='high', thr=('pag', 'minimum'),
         s='main', socx='yes | [40.10.1.1.1.2] Supplement for those on low income: Statutory pension insurance + '
                        '[40.10.2.1.1.2] Widow(er)s pension (means-tested)',
         note='SOCX programme values; incapacity adds [40.10.3.1.1.2]. The part of the AZ inside tau.'),
    # ---------------- DEU
    dict(iso3='DEU', scheme='Grundsicherung im Alter und bei Erwerbsminderung (SGB XII ch. 4), net expenditure',
         year=2019, nc=6868.460, rec=1085043, base='recipients, December 2019 (all ages from 18)', raw='deuBT',
         conf='high', q=['2019 556.864 6.868.460'], thr=('pag', 'targeted'), s=None,
         socx='part | old age: Other cash benefits (Social assistance) [276.10.1.1.3.5]; rest in incapacity '
              '[276.10.3.1.5.19] and other social policy areas [276.10.9.1.1.1]',
         note='Gross 7,425.324, receipts 556.864 (Laender returns under s. 46a SGB XII, as of 2 May 2024). Includes '
              'the reduced-earning-capacity part; 51.8% of December 2019 recipients were above the age limit '
              '(561,969).', q2=('deuDS', 'Dezember 2019 1 085 043')),
    dict(iso3='DEU', scheme='Grundsicherung im Alter und bei Erwerbsminderung (SGB XII ch. 4), net expenditure',
         year=2020, nc=7561.666, rec=1098625, base='recipients, December 2020 (all ages from 18)', raw='deuBT',
         conf='high', q=['2020 325.490 7.561.666'], thr=('pag', 'targeted'), s='alt',
         socx='part | old age: Other cash benefits (Social assistance) [276.10.1.1.3.5]; rest in incapacity '
              '[276.10.3.1.5.19] and other social policy areas [276.10.9.1.1.1]',
         note='Gross 7,887.157, receipts 325.490. Includes the reduced-earning-capacity part. Destatis: 564,110 of '
              'the December 2020 recipients above the age limit (51.3%; 3.1% of 65+); PaG table 3.2 prints 0.5% of 65+.',
         q2=('deuDS', '2020 2 1 098 625 534 520 305 115 229 405 564 110')),
    dict(iso3='DEU', scheme='Grundsicherung im Alter (age-limit part only), net expenditure', year=2020, nc=3300,
         rec=564110, base='recipients above the age limit, December 2020 (Destatis)', raw='deuKOM', conf='low',
         q=['rund 3,3 Milliarden Euro'], thr=('pag', 'targeted'), s='alt',
         socx='part | SOCX old age books 1,810.7 under [276.10.1.1.3.5]',
         note='"rund 3,3 Milliarden", a government answer as reported by BILD via kommunal.de; the Drucksache was not '
              'found. Recipients from Destatis PD21_188.'),
    dict(iso3='DEU', scheme='SOCX old age: Other cash benefits (Social assistance)', year=2019,
         socxcodes=['276.10.1.1.3.5'], rec=None, base='', conf='high', thr=('pag', 'targeted'), s=None,
         socx='yes | [276.10.1.1.3.5]',
         note='SOCX also books in old age Other cash benefits (means tested) [276.10.1.1.3.9] and, in kind, '
              'Accommodation for adults aged 65 and over (Social assistance) [276.10.1.2.1.3]; their content is not '
              'described.'),
    dict(iso3='DEU', scheme='SOCX old age: Other cash benefits (Social assistance)', year=2020,
         socxcodes=['276.10.1.1.3.5'], rec=None, base='', conf='high', thr=('pag', 'targeted'), s='main',
         socx='yes | [276.10.1.1.3.5]',
         note='The part of the old-age social assistance inside tau; incapacity [276.10.3.1.5.19] and other social '
              'policy areas [276.10.9.1.1.1] hold more social-assistance cash.'),
    # ---------------- BEL
    dict(iso3='BEL', scheme='GRAPA / IGO (income guarantee for the elderly)', year=2019, nc=593.7, rec=None, base='',
         raw='belRH', conf='high', q=['jaarlijkse uitgaven voor de IGO in 2019 nagenoeg 593,7 miljoen euro'],
         thr=('pag', 'targeted'), s=None, socx='yes | Minimum Guaranteed Income for the aged [56.10.1.1.1.15]',
         note='Court of Audit advice on bill 55 0416. SOCX 2019: see the SOCX row. PaG table 3.2: 5.05% of 65+ (2020).'),
    dict(iso3='BEL', scheme='GRAPA / IGO (SOCX)', year=2019, socxcodes=['56.10.1.1.1.15'], rec=None, base='',
         conf='high', thr=('pag', 'targeted'), s=None, socx='yes | [56.10.1.1.1.15]', note='SOCX programme value.'),
    dict(iso3='BEL', scheme='GRAPA / IGO (SOCX)', year=2020, socxcodes=['56.10.1.1.1.15'], rec=None, base='',
         conf='high', thr=('pag', 'targeted'), s='main', socx='yes | [56.10.1.1.1.15]',
         note='SOCX programme value; no national 2020 outturn found.'),
    dict(iso3='BEL', scheme='Minimum pension supplement (minimum garanti / gewaarborgd minimum), employees and '
         'self-employed', year=2020, nc=None, rec=None, base='', raw=None, conf='', thr=('pag', 'minimum'), s=None,
         socx='part | inside the SFP/ONP old-age pensions ([56.10.1.1.1.13], [56.10.1.1.1.22], [56.10.1.1.1.23])',
         note='GAP: no cost or count found (SFP annual reports and statistics not reached; Court of Audit, Planning '
              'Bureau and parliamentary documents searched). PaG table 3.2: 32.5% of 65+; PAR 2021: about 1 in 5 '
              'women, 1 in 20 men (2013). Threshold 32.5% AW against PaG\'s 33.7%: not reached, by 1.2 points.'),
    # ---------------- PRT
    dict(iso3='PRT', scheme='Complementos sociais (top-up of contributory pensions to the minimum: old age, '
         'survivors, invalidity)', year=2019, nc=1910, rec=None, base='', raw='prtCFP', conf='high',
         q=['Complementos sociais 1 910 1 946 1,9 1 927'], thr=('pag', 'minimum'), s=None,
         socx='part | inside old age pension (non means tested) [620.10.1.1.1.21] and survivors [620.10.2.1.1.19]; '
              'the invalidity part is SOCX incapacity',
         note='Conta da Seguranca Social 2019 (CSS/19 column), Subsistema de Solidariedade, financed by the State '
              'budget; covers invalidity too (no split found), so an upper bound. Recipients not found.'),
    dict(iso3='PRT', scheme='Complementos sociais (top-up of contributory pensions to the minimum: old age, '
         'survivors, invalidity)', year=2020, nc=1927, rec=None, base='', raw='prtCFP', conf='high',
         q=['Complementos sociais 1 910 1 946 1,9 1 927'], thr=('pag', 'minimum'), s='main',
         socx='part | inside old age pension (non means tested) [620.10.1.1.1.21] and survivors [620.10.2.1.1.19]; '
              'the invalidity part is SOCX incapacity',
         note='2020 execution, provisional (CFP, IGFSS data). Upper bound of the old-age and survivors part. CGA '
              '(civil servants) not covered. Recipients not found.'),
    dict(iso3='PRT', scheme='Complemento Solidario para Idosos (CSI)', year=2019, nc=216, rec=None, base='',
         raw='prtCFP', conf='high', q=['Complemento Solidário para Idosos 216 224 3,8 212'],
         thr=('pag', 'targeted'), s=None,
         socx='yes | Other cash periodic benefits (means tested) [620.10.1.1.3.14]', note='CSS/19 column.'),
    dict(iso3='PRT', scheme='Complemento Solidario para Idosos (CSI)', year=2020, nc=212, rec=171768,
         base='beneficiaries with at least one payment in the year (provisional)', raw='prtCFP', conf='high',
         q=['Complemento Solidário para Idosos 216 224 3,8 212'], thr=('pag', 'targeted'), s='main',
         socx='yes | Other cash periodic benefits (means tested) [620.10.1.1.3.14]',
         note='2020 execution, provisional. Recipients from Seguranca Social em Numeros 2023 (2020 column; 2021 '
              '167,399, 2022 162,646).',
         q2=('prtSS', '931 424 171 768')),
    dict(iso3='PRT', scheme='Pensoes sociais (non-contributory social pensions)', year=2019, nc=121, rec=None,
         base='', raw='prtCFP', conf='high', q=['Pensões sociais 121 119 -1,1 118'], thr=None, reach='no', s=None,
         socx='yes | old age pension (means tested) [620.10.1.1.1.22]',
         note='CSS/19 column; includes invalidity social pensions. Means-tested: income limit not collected, but a '
              'full-career contributory pension (PaG 0.5 AW: 38.1% AW) excludes it.'),
    dict(iso3='PRT', scheme='Pensoes sociais (non-contributory social pensions)', year=2020, nc=118, rec=None,
         base='', raw='prtCFP', conf='high', q=['Pensões sociais 121 119 -1,1 118'], thr=None, reach='no',
         s='main', socx='yes | old age pension (means tested) [620.10.1.1.1.22]',
         note='2020 execution, provisional; includes invalidity social pensions.'),
    # ---------------- NOR
    dict(iso3='NOR', scheme='Minimum pension level supplements (saertillegg, pensjonstillegg, minstenivaatillegg; old '
         'rules) and garantipensjon (new rules)', year=2020, nc=None, rec=138200,
         base='minimum pensioners among old-age pensioners, end of December 2020 (14.1%)', raw='norNAV', conf='high',
         q=['Ved utgangen av desember 2020 var det 138 200 minstepensjonister'], thr=('pag', 'targeted'), s=None,
         socx='part | inside Personal basic pension [578.10.1.1.1.1] and Supplementary pension [578.10.1.1.1.3]; '
              'no separate line',
         note='Spending GAP (Statsbudsjett at regjeringen.no returns 403). Reached at 0.5 AW: PaG\'s pension there '
              '(27.1% AW) is below the guarantee level (30.5% AW) and PaG\'s replacement rate at 0.5 AW (54.1) exceeds '
              'the one at 1 AW (39.4), so it is already inside the formula slope and s = 0 by definition.'),
    # ---------------- POL
    dict(iso3='POL', scheme='Doplata do emerytury minimalnej (State-budget top-up, new-system pensions)', year=2021,
         nc=180.8, rec=72200, base='new-system pensioners topped up in the year', raw='polPRAWO', conf='low',
         q=['72,2 tys. emerytów z nowego systemu pochłonęło prawie 180,8 mln zł'],
         thr=('pag', 'minimum'), s='main',
         socx='part | paid with the ZUS pension (Employees and self-employed [616.10.1.1.1.1]); no separate line',
         note='No 2020 value published; 2021 is the nearest year (2017: 101.5 for 38,100). New-system pensions only; '
              'the top-up of old-system pensions to the minimum is paid inside FUS pensions and not published. Reached '
              'at 0.5 AW with 2020 values (minimum 23.6% AW above PaG\'s 15.9% AW), but PaG\'s replacement rates '
              '(31.8 at 0.5 AW, 30.6 at 1 AW) show it not binding in PaG\'s projection, where the minimum is '
              'price-indexed: the definition of s does not settle this case.'),
    dict(iso3='POL', scheme='Old-age pensions paid at the minimum (all systems)', year=2020, nc=None, rec=221000,
         base='ZUS old-age pensions at the minimum after the March 2020 indexation (excl. international agreements '
              'and farmers)', raw='polZUS', conf='high', q=['Emerytura 221,0'], thr=('pag', 'minimum'), s=None,
         socx='part | inside [616.10.1.1.1.1] and [616.10.1.1.1.6]',
         note='Spending GAP. Table 2: 529.1 thousand ZUS benefits at the minimum, of which old-age 221.0, disability '
              '182.6, survivors 125.5 (raw text order). PaG table 3.2: 1.5% '
              'of 65+.'),
    # ---------------- LUX
    dict(iso3='LUX', scheme='Complement pour pension minimum (general scheme)', year=2023, nc=None, rec=None,
         base='', raw='luxPJ', conf='low', q=['de 19,5% en 2009 à 16,4% en 2023'], thr=('pag', 'minimum'), s=None,
         socx='part | inside Old-age pension [442.10.1.1.1.5] and Survivors\' pension [442.10.2.1.1.8]',
         note='GAP for 2019-2020: spending and count not found (IGSS rapport general pages return 404; the 2022 '
              'bilan technique handout gives the formula only). Share of pensions with the complement: 16.4% in 2023, '
              '19.5% in 2009 (ministry via Paperjam).'),
    # ---------------- HUN
    dict(iso3='HUN', scheme='Minimum old-age pension top-up (HUF 28,500 a month since 2008): upper bound', year=2020,
         hunbound=True, base='PaG table 3.2: 0.61% of 65+, times Eurostat 65+ on 1 January 2021', raw='hunPOP',
         conf='low (derived upper bound)', thr=('pag', 'minimum'), s='main',
         socx='part | inside Old age pension (non means-tested) [348.10.1.1.1.1]',
         note='Upper bound = recipients x full minimum x 12; the top-up is the gap to the computed pension, smaller. '
              'ONYF counts not reached.'),
    dict(iso3='HUN', scheme='SOCX old age: Other cash periodic benefits (means-tested)', year=2019,
         socxcodes=['348.10.1.1.3.4'], rec=None, base='', conf='medium', thr=('pag', 'targeted'), s=None,
         socx='yes | [348.10.1.1.3.4]', note='Presumably the old-age allowance (idoskoruak jaradeka); SOCX does not '
                                              'name it.'),
    dict(iso3='HUN', scheme='SOCX old age: Other cash periodic benefits (means-tested)', year=2020,
         socxcodes=['348.10.1.1.3.4'], rec=None, base='', conf='medium', thr=('pag', 'targeted'), s='main',
         socx='yes | [348.10.1.1.3.4]', note='Presumably the old-age allowance (idoskoruak jaradeka). SOCX old age '
                                              'also has Other cash lump sum benefits (means-tested) [348.10.1.1.3.2].'),
]
COLS = ['iso3', 'scheme', 'year', 'spending_nc', 'currency', 'spending_pctGDP', 'recipients', 'recipients_base',
        'inSOCXoldAge', 'reachedAtHalfAW', 'source', 'url', 'retrieved', 'confidence', 'note']


def _p(key):
    return os.path.join(RAWDIR, RAW[key][0])


def fetch(refetch = False):
    """ Download every raw file missing (or all, with refetch); a pdf also gets its pdftotext -layout .txt. """
    for key, (rel, url, _) in RAW.items():
        p = os.path.join(RAWDIR, rel)
        if os.path.exists(p) and not refetch:
            continue
        os.makedirs(os.path.dirname(p), exist_ok = True)
        print('fetching', rel)
        with urllib.request.urlopen(urllib.request.Request(url, headers = UA), timeout = 120) as r, \
                open(p, 'wb') as f:
            f.write(r.read())
    for key, (rel, url, _) in RAW.items():
        p = os.path.join(RAWDIR, rel)
        for mode, q in (('-layout', p[:-4] + '.txt'), ('-raw', p[:-4] + '_raw.txt')):
            if p.endswith('.pdf') and not os.path.exists(q):
                if not shutil.which('pdftotext'):
                    raise RuntimeError('pdftotext not found: needed once to write ' + q)
                subprocess.run(['pdftotext', '-enc', 'UTF-8', mode, p, q], check = True)


def _norm(s):
    return re.sub(r'\s+', ' ', s).strip()


def _text(key):
    """ The raw file's text, whitespace collapsed: the .txt (and _raw.txt) of a pdf, the tag-stripped html. """
    p = _p(key)
    if p.endswith('.pdf'):
        parts = [open(q, encoding = 'utf-8', errors = 'replace').read() for q in (p[:-4] + '.txt', p[:-4] + '_raw.txt')
                 if os.path.exists(q)]
        return _norm(' || '.join(parts))
    t = open(p, encoding = 'utf-8', errors = 'replace').read()
    return _norm(html.unescape(re.sub(r'<[^>]+>', ' ', t)))


def check():
    """ Every quote must appear in its raw file: the typed numbers are the published ones. """
    bad = []
    for r in ROWS:
        for key, q in [(r.get('raw'), x) for x in r.get('q', [])] + ([r['q2']] if 'q2' in r else []):
            if _norm(q) not in _text(key):
                bad.append((r['iso3'], r['year'], key, q))
    if bad:
        raise ValueError('quotes not found in the raw files:\n' + '\n'.join(map(str, bad)))


def gdp():
    j = json.load(open(_p('gdp'), encoding = 'utf-8'))
    dims, size = j['id'], j['size']
    out = {}
    for iso, g in ISO.items():
        for t in j['dimension']['time']['category']['index']:
            idx = {'geo': j['dimension']['geo']['category']['index'][g], 'time': j['dimension']['time']['category']['index'][t]}
            k = 0
            for d, s in zip(dims, size):
                k = k * s + idx.get(d, 0)
            v = j['value'].get(str(k))
            if v is not None:
                out[(iso, int(t))] = float(v)
    return out, j.get('updated', '')


def hunPop65():
    j = json.load(open(_p('hunPOP'), encoding = 'utf-8'))
    dims, size = j['id'], j['size']
    cat = {d: j['dimension'][d]['category']['index'] for d in dims}
    tot = 0.0
    for a in cat['age']:
        idx = {'age': cat['age'][a], 'time': cat['time']['2021'], 'geo': cat['geo']['HU']}
        k = 0
        for d, s in zip(dims, size):
            k = k * s + idx.get(d, 0)
        tot += float(j['value'][str(k)])
    return tot


def socxProgrammes():
    """ SOCX programme values (millions of national currency) by programme code, 2019 and 2020. """
    w = pd.read_excel(SOCX_XLSX, sheet_name = 'SOCX_NatCur', engine = 'openpyxl')
    w = w[w.Code.isin({c for r in ROWS for c in r.get('socxcodes', [])})]
    return {(r['Code'], y): float(r[y]) for r in w.to_dict('records') for y in (2019, 2020)}


def socxTP01():
    d = pd.read_csv(SOCX_AGG)
    d = d[(d.EXPEND_SOURCE == 'ES10') & (d.SPENDING_TYPE == 'C') & (d.PROGRAMME_TYPE == 'TP01')
          & (d.UNIT_MEASURE == 'XDC') & (d.PRICE_BASE == 'V')]
    return {(r.REF_AREA, int(r.TIME_PERIOD)): float(r.OBS_VALUE) * 10 ** (int(r.UNIT_MULT) - 6)
            for r in d.itertuples(index = False)}


def pag():
    """ PaG's public pension at 0.5 AW (% of AW: table 4.2, men, times 0.5) and table 3.2's first-tier levels. """
    d = pd.read_csv(TIERS, comment = '#', dtype = {'year': str})
    d = d[d.year == '2020']
    g = d[(d.item == 'grr_tier') & (d.tier == 'public_mandatory') & (d.sex == 'M') & (d.earnings_multiple == 0.5)]
    p05 = {r.iso3: 0.5 * float(r.value) for r in g.itertuples(index = False)}
    f = d[d.item == 'firstTier_value']
    lvl = {(r.iso3, r.tier): float(r.value) for r in f.itertuples(index = False)}
    return p05, lvl


def collect():
    check()
    G, gupd = gdp()
    S = socxProgrammes()
    p05, lvl = pag()
    out, roles = [], []
    for r in ROWS:
        iso, y = r['iso3'], r['year']
        if 'socxcodes' in r:
            nc = round(sum(S[(c, y)] for c in r['socxcodes']), 3)
            src, url = 'OECD SOCX programme-level data (SOCX2025-All-nat-cur.xlsx)', SOCX_URL
        elif r.get('hunbound'):
            rec = round(0.0061 * hunPop65())
            nc = round(rec * 28500 * 12 / 1e6, 1)
            r = dict(r, rec = rec)
            src, url = RAW['hunPOP'][2] + '; ' + PAG, RAW['hunPOP'][1]
        else:
            nc = r['nc']
            src, url = (RAW[r['raw']][2], RAW[r['raw']][1]) if r.get('raw') else ('', '')
        if 'q2' in r:
            src, url = src + '; ' + RAW[r['q2'][0]][2], url + ' ; ' + RAW[r['q2'][0]][1]
        thr = r['thr']
        if isinstance(thr, tuple):
            t = lvl[(iso, thr[1])]
            how = '{} {} {:.1f}% AW'.format(PAG, thr[1], t)
        elif thr is not None:
            t, how = thr, 'threshold {:.1f}% AW (note)'.format(thr)
        if thr is None:
            reach, rnote = r['reach'], ''
        else:
            reach = 'yes' if p05[iso] < t else 'no'
            rnote = ' Reach: PaG 2021 public pension at 0.5 AW {:.1f}% AW vs {}.'.format(p05[iso], how)
        pct = '' if nc is None or (iso, y) not in G else round(100 * nc / G[(iso, y)], 4)
        out.append(dict(iso3 = iso, scheme = r['scheme'], year = y, spending_nc = '' if nc is None else nc,
                        currency = CUR[iso] + ' million', spending_pctGDP = pct,
                        recipients = '' if r.get('rec') is None else int(r['rec']), recipients_base = r['base'],
                        inSOCXoldAge = r['socx'], reachedAtHalfAW = reach, source = src, url = url, retrieved = RET,
                        confidence = r['conf'], note = r['note'] + rnote))
        roles.append(r['s'])
    return out, roles, gupd


def sTable(out, roles):
    """ s by country, 2020 rows (POL: 2021): main numerator and the alternatives, over SOCX TP01 public cash. """
    T = socxTP01()
    lines = []
    for iso in ISO:
        rr = [(o, s) for o, s in zip(out, roles) if o['iso3'] == iso and s in ('main', 'alt') and o['spending_nc'] != '']
        main = sum(o['spending_nc'] for o, s in rr if s == 'main')
        tp = T[(iso, 2020)]
        alts = ['{} {:.2%}'.format(o['scheme'][:38], o['spending_nc'] / tp) for o, s in rr if s == 'alt']
        if not any(s == 'main' for o, s in rr):
            reached = any(o['reachedAtHalfAW'] == 'yes' for o in out if o['iso3'] == iso)
            lines.append('{} s = 0: its minimum is reached at 0.5 AW'.format(iso) if reached else
                         '{} s: gap, no top-up spending found'.format(iso))
            continue
        lines.append('{} s = {:.4f}  (numerator {:,.1f}, TP01 2020 {:,.1f} {}m){}'.format(
            iso, main / tp, main, tp, CUR[iso], ('; alt: ' + '; '.join(alts)) if alts else ''))
    return '\n'.join(lines)


def write(out, gupd, path = OUT):
    with open(path, 'w', newline = '', encoding = 'utf-8') as f:
        f.write('# DRAFT, not wired (2026-10-30): minimum-pension top-ups and tested old-age benefits, ITA FRA AUT DEU '
                'BEL PRT NOR POL LUX HUN; spending_nc in millions of national currency, spending_pctGDP against '
                'Eurostat nama_10_gdp CP_MNAC (updated {}); written by python/paper/oecd26Topups.py; see '
                'data/README.md. Read with comment="#".\n'.format(gupd))
        w = pd.DataFrame(out, columns = COLS)
        w.to_csv(f, index = False, lineterminator = '\n')


def main():
    p = argparse.ArgumentParser(description = __doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    p.add_argument('--force', action = 'store_true', help = 'rebuild the csv even if it exists')
    p.add_argument('--refetch', action = 'store_true', help = 'download every raw file again (implies --force)')
    a = p.parse_args()
    fetch(a.refetch)
    out, roles, gupd = collect()
    if os.path.exists(OUT) and not (a.force or a.refetch):
        print('up to date: {} (--force to rebuild)'.format(os.path.relpath(OUT, C.REPO)))
    else:
        write(out, gupd)
        print('{} rows -> {}'.format(len(out), os.path.relpath(OUT, C.REPO)))
    print(sTable(out, roles))


if __name__ == '__main__':
    main()
