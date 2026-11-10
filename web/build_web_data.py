"""OligoScout web arayuzu icin veri paketi.

reference_py/out altindaki GERCEK motor ciktilarini okur ve web/oligoscout_data.json
uretir. Hicbir deger burada hesaplanmaz/uydurulmaz; yalniz secilir ve sikistirilir.
Kullanim:  python web/build_web_data.py
"""
import csv, json, os, collections, hashlib, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'reference_py'))
OUT = os.path.join(ROOT, 'reference_py', 'out')
P = lambda *a: os.path.join(OUT, *a)

def tsv(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def jsonl(path):
    with open(path, encoding='utf-8') as f:
        return [json.loads(l) for l in f if l.strip()]

def num(x):
    try:
        v = float(x); return int(v) if v.is_integer() and 'e' not in str(x).lower() and '.' not in str(x) else v
    except Exception:
        return x

def maybe_json(x):
    if isinstance(x, str) and x[:1] in '{[':
        try: return json.loads(x)
        except Exception: return x
    return x

def text(path, limit=6000):
    try:
        with open(path, encoding='utf-8') as f: return f.read()[:limit]
    except Exception: return ''

RC = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def revcomp(s): return ''.join(RC.get(c,'N') for c in reversed(s.upper()))

# ---------- sablon ----------
fa = open(os.path.join(ROOT, 'hedef_chr22.fa'), encoding='utf-8').read().splitlines()
header = fa[0][1:]
template = ''.join(l.strip() for l in fa[1:]).upper()
OFFSET = 20000000  # 0-tabanli baslangic (baslikta beyan)

# ---------- katalog ----------
catalog = [{'id': r['kural'], 'fam': r['aile'], 'title': r['spec_basligi'], 'status': r['durum'], 'where': r['nerede']}
           for r in tsv(P('durum', 'qc_katalog_kapsami.tsv'))]
open_items = tsv(P('durum', 'acik_maddeler.tsv'))

# ---------- panel (ctDNA, tum genom taramasi) ----------
PD = P('panel_design_ctdna')
ranked = tsv(os.path.join(PD, 'final_ranked.tsv'))
panelqc = {r['pair_id']: r for r in tsv(os.path.join(PD, 'panel_qc.tsv'))}
contract = {r['field']: r['value'] for r in tsv(os.path.join(PD, 'panel.contract.tsv'))}
# guncel QC katmani (SPC-01..09, COND) varsa onu kullan: reference_py/out/web_panel_ctdna
QCDIR = P('web_panel_ctdna') if os.path.exists(P('web_panel_ctdna', 'panel_qc_records.jsonl')) else PD
_all = jsonl(os.path.join(QCDIR, 'panel_qc_records.jsonl'))
prim_recs = {r['primer_id']: r for r in _all if 'primer_id' in r}
pair_recs = collections.defaultdict(list); panel_recs = []
for r in _all:
    if 'primer_id' in r: continue
    (panel_recs if r.get('target') == '__panel__' else pair_recs[r.get('target')]).append(r)
panelqc = {r['pair_id']: r for r in tsv(os.path.join(QCDIR, 'panel_qc.tsv'))}
hetero = {(r['left_primer'], r['right_primer']): r for r in tsv(P('design_ctdna', 'design_pairs.tsv'))}
offt = collections.defaultdict(list)
for r in tsv(os.path.join(PD, 'offtarget_products.tsv')):
    offt[r['pair_id']].append({k: num(v) for k, v in r.items() if k != 'pair_id'})
products = collections.defaultdict(list)
for r in tsv(os.path.join(PD, 'panel.products.tsv')):
    products[r['pair_id']].append({k: num(v) for k, v in r.items() if k != 'pair_id'})

# baglanma yerleri: primer basina histogram + en yakin 40 site
bind = collections.defaultdict(list)
with open(os.path.join(PD, 'panel.bindings.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        bind[(r['pair_id'], r['role'])].append(r)

def pack_bindings(rows, primer):
    hist = collections.Counter()
    sites = []
    for r in rows:
        mm = int(r['mismatches']); mask = r['mask']
        L = len(mask)
        anchored = 'X' not in mask[max(0, L-12):] if L else True  # 3' 12-nt capa (sozlesme anchor_length)
        hist[(mm, anchored)] += 1
        site = r['site_plus_strand'].upper()
        oriented = site if r['strand'] == '+' else revcomp(site)
        sites.append({'chrom': r['chrom'], 'start': int(r['start']), 'strand': r['strand'], 'mm': mm,
                      'mask': mask, 'site': oriented, 'anchored': anchored})
    sites.sort(key=lambda s: (s['mm'], not s['anchored'], s['chrom'], s['start']))
    return {'n': len(rows), 'hist': [{'mm': k[0], 'anchored': k[1], 'n': v} for k, v in sorted(hist.items())],
            'top': sites[:40]}

# seqpair (dizi/cift QC) kayitlari: (F,R) ile eslestir
seqpair = {}
for r in jsonl(P('design_ctdna', 'seqpair_qc_records.jsonl')):
    if 'checks' in r: seqpair[(r['forward'], r['reverse'])] = r

def compact_check(c, role=None):
    o = {'id': c['id'], 'sev': c.get('severity'), 'value': c.get('value'), 'thr': c.get('threshold'),
         'note': c.get('note', '')}
    if role or c.get('role'): o['role'] = role or c.get('role')
    m = c.get('model')
    if isinstance(m, dict): o['model'] = m
    return o

# varyant kayitlari
var_recs = collections.defaultdict(list)
for r in jsonl(P('variant_design_ctdna', 'variant_qc_records.jsonl')):
    var_recs[r['target']].append(compact_check(r))
var_sites = {r['site_id']: {k: num(v) for k, v in r.items()} for r in tsv(P('variant_design_ctdna', 'variant_sites.tsv'))}
# variant_design pair id'leri design_ctdna/design_panel ile ayni; sirala dizilere gore
dpanel = {r['pair_id']: r for r in tsv(P('design_ctdna', 'design_panel.tsv'))}
seq2dpid = {(r['primer1'], r['primer2']): pid for pid, r in dpanel.items()}

detect = collections.defaultdict(list)
for r in tsv(P('design_ctdna', 'ctdna_detectability.tsv')):
    detect[r['pair_id']].append({k: num(v) for k, v in r.items() if k != 'pair_id'})

try:
    from oligoscout_thermo import align as _align
except Exception:
    _align = None
def cur_hetero(F, R):
    if _align is None: return None
    out = {}
    for mode in ('any', '3prime'):
        d = _align.best_dimer(F, R, mode)
        out[mode] = d.as_dict() if getattr(d, 'found', False) else {'found': False}
    return out

pairs = []
for r in ranked:
    pid = r['pair_id']; F = r['left_primer']; R = r['right_primer']
    fpos = template.find(F); rpos = template.find(revcomp(R))
    sp = seqpair.get((F, R))
    dpid = seq2dpid.get((F, R))
    pq = panelqc.get(pid, {})
    item = {k: num(v) for k, v in r.items()}
    item.update({
        'fStart': fpos, 'rStart': rpos, 'rEnd': rpos + len(R) - 1 if rpos >= 0 else -1,
        'panelqc': {k: num(v) for k, v in pq.items()},
        'seqpair': [compact_check(c) for c in sp['checks']] if sp else [],
        'seqpair_worst': sp['worst_severity'] if sp else None,
        'primers': {},
        'bindings': {'F': pack_bindings(bind[(pid, 'F')], F), 'R': pack_bindings(bind[(pid, 'R')], R)},
        'offtarget': offt.get(pid, []),
        'products': products.get(pid, []),
        'design_id': dpid,
        'variants': {'F': var_recs.get(f'{dpid}_F', []), 'R': var_recs.get(f'{dpid}_R', []), 'pair': var_recs.get(f'{dpid}_amplicon', []),
                     'sites': [var_sites.get(f'{dpid}_F'), var_sites.get(f'{dpid}_R')]} if dpid else None,
        'detect': detect.get(dpid, []) if dpid else [],
        'pair_checks': [compact_check(c) for c in pair_recs.get(pid, [])],
        'hetero_current': cur_hetero(F, R),
        'hetero': {k: num(v) for k, v in hetero.get((F, R), {}).items() if k in ('heterodimer_dG37', 'heterodimer3p_dG37', 'uncalibrated_score', 'pareto_rank', 'product_start', 'product_end')},
    })
    for role in ('F', 'R'):
        pr = prim_recs.get(f'{pid}_{role}')
        if pr:
            item['primers'][role] = {'thermo': pr['thermo'], 'worst': pr['worst_severity'],
                                     'not_evaluated': pr.get('not_evaluated', {}),
                                     'checks': [compact_check(c) for c in pr['checks']]}
    pairs.append(item)

# ---------- gnomAD varyant onbellegi (sablon bolgesi) ----------
variants = []; _seen = set()
for v in [x for x in jsonl(P('variant_cache.jsonl')) if 'pos' in x]:
    _k = (v['pos'], v['ref'], v['alt'])
    if _k in _seen: continue
    _seen.add(_k)
    pos = int(v['pos']); i1 = pos - 1 - OFFSET  # VCF 1-tabanli
    ref = v['ref']
    ok1 = template[i1:i1+len(ref)] == ref if 0 <= i1 < len(template) else False
    variants.append({'pos': pos, 'idx': i1, 'ref': ref, 'alt': v['alt'], 'af': v.get('af'), 'an': v.get('an'),
                     'rsid': v.get('rsid'), 'pop': {k: round(x, 6) for k, x in (v.get('pop_af') or {}).items()},
                     'ref_ok': ok1})

# ---------- multipleks ----------
MX = os.path.join(PD, 'multiplex')
multiplex = {
    'summary': tsv(os.path.join(MX, 'multiplex_summary.tsv')),
    'dimer': [{k: num(v) for k, v in r.items()} for r in tsv(os.path.join(MX, 'dimer_matrix.tsv'))],
    'combos': [{k: num(v) for k, v in r.items()} for r in tsv(os.path.join(MX, 'cross_combinations.tsv'))],
    'cross': [{k: num(v) for k, v in r.items()} for r in tsv(os.path.join(MX, 'cross_products.tsv'))],
}

# ---------- prob / blocker / erime ----------
probes = [{k: num(v) for k, v in r.items()} for r in tsv(P('design_ddpcr', 'probe_panel.tsv'))]
ddpcr_pairs = {r['pair_id']: r for r in tsv(P('design_ddpcr', 'design_panel.tsv'))}
for p in probes:
    d = ddpcr_pairs.get(p['pair_id']);
    if d: p['F'] = d['primer1']; p['R'] = d['primer2']; p['amplicon'] = num(d['expected_amplicon_length'])
blocker = [{k: maybe_json(v) for k, v in r.items()} for r in tsv(P('blocker', 'ct05_qc.tsv'))]
melt_sum = [{k: num(v) for k, v in r.items()} for r in tsv(P('melt_qpcr_sybr', 'amplicon_melt.tsv'))]
curves = collections.defaultdict(list)
for r in tsv(P('melt_qpcr_sybr', 'amplicon_melt_curves.tsv')):
    t = float(r['temperature_c'])
    if abs((t*2) - round(t*2)) < 1e-6:  # 0.5 C adim
        curves[r['id']].append([t, round(float(r['theta']), 4), round(float(r['minus_dtheta_dt']), 5)])
qpcr_panel = {r['pair_id']: r for r in tsv(P('design_qpcr_sybr', 'design_panel.tsv'))}
for m in melt_sum:
    d = qpcr_panel.get(m['id'])
    if d: m['F'] = d['primer1']; m['R'] = d['primer2']

# ---------- dogrulama ----------
validation = {
    'crosscheck': text(P('param_crosscheck', 'crosscheck_summary.md')),
    'diff_cpp': text(P('diff_cpp', 'diff_cpp_summary.md')),
    'terminal': text(P('terminal_audit', 'terminal_audit_summary.md')),
    'salt': text(P('salt_crosscheck', 'salt_crosscheck_summary.md')),
    'kontrol': text(P('kontrol', 'KONTROL_RAPORU.md')),
    'esik': text(P('kalibrasyon', 'esik_karar_notu.md')),
    'th0204': tsv(P('kalibrasyon', 'th02_th04_kanit.tsv')),
    'regresyon': {fn: text(P('regresyon', fn), 1600) for fn in sorted(os.listdir(P('regresyon'))) if fn.endswith('.log')},
}
profiles = text(P('profile_comparison.tsv'))
import sys
sys.path.insert(0, os.path.join(ROOT, 'reference_py'))
try:
    from oligoscout_thermo import profiles as _prof
    app_profiles = [v.as_dict() for v in _prof.PROFILES.values()]
except Exception as e:  # motor import edilemezse bos
    app_profiles = []

# motor parametre tablolari (tarayici motoru AYNI tablolari kullanir; kopya degil, disa aktarim)
def _prov(t):
    pv = getattr(t, 'provenance', None); out = {}
    if pv is not None:
        d = pv.__dict__ if hasattr(pv, '__dict__') else {}
        for k, v in d.items(): out[k] = v if isinstance(v, (str, int, float, bool, type(None))) else str(v)
    out['sha256'] = getattr(t, 'sha256', None); out['available'] = bool(getattr(t, 'is_available', False))
    return out
try:
    from oligoscout_thermo import params as _pm
    engine_params = {
        'wc': _pm.WC_SANTALUCIA_1998.require(), 'wc_prov': _prov(_pm.WC_SANTALUCIA_1998),
        'mismatch': _pm.MISMATCH_ALLAWI_PEYRET.require(), 'mismatch_prov': _prov(_pm.MISMATCH_ALLAWI_PEYRET),
        'loop': _pm.DNA_LOOP_SANTALUCIA_HICKS_2004.require(), 'loop_prov': _prov(_pm.DNA_LOOP_SANTALUCIA_HICKS_2004),
        'at_bulge_closure_dg37': _pm.AT_BULGE_CLOSURE.value,
        'R_cal': _pm.R_CAL, 'kelvin': _pm.KELVIN,
        'registry_digest': _pm.registry_digest() if hasattr(_pm, 'registry_digest') else None,
        'tables': [dict(name=n, **_prov(getattr(_pm, n))) for n in dir(_pm) if n.isupper() and hasattr(getattr(_pm, n), 'provenance')],
        'model_constants': {k: {kk: (vv if isinstance(vv, (str, int, float)) else str(vv)) for kk, vv in v.__dict__.items()} for k, v in getattr(_pm, 'MODEL_CONSTANTS', {}).items()},
    }
except Exception as e:
    engine_params = {'error': str(e)}

# tasarim denkligi parmak izi (tarayici motoru canli karsilastirir)
design_check = None
_dc = P('web_design_check')
if os.path.exists(os.path.join(_dc, 'design_pairs.tsv')):
    rows = tsv(os.path.join(_dc, 'design_pairs.tsv'))
    rejc = collections.Counter(r['rule'] for r in tsv(os.path.join(_dc, 'design_rejections.tsv')))
    design_check = {'source': 'reference_py/out/web_design_check (design_primers.py --profile ctdna, hedef_chr22.fa)',
                    'n_pairs': len(rows), 'rejections': dict(rejc),
                    'top': [[r['left_primer'], r['right_primer'], float(r['uncalibrated_score']), int(r['pareto_rank'])] for r in rows[:60]]}

data = {
    'meta': {'header': header, 'offset': OFFSET, 'template': template, 'chrom': 'chr22',
             'param_sha': pairs[0]['primers'].get('F', {}).get('checks', [{}])[0].get('model', {}).get('param_set_sha256', '') if pairs else '',
             'seqpair_param_sha': next(iter(seqpair.values()))['param_set_sha256'] if seqpair else '',
             'conditions': next(iter(seqpair.values()))['conditions'] if seqpair else {},
             'profile': 'ctdna', 'contract': contract, 'qc_dir': os.path.relpath(QCDIR, ROOT),
             'panel_checks': [compact_check(c) for c in panel_recs],
             'source_note': 'reference_py/out/panel_design_ctdna + design_ctdna + variant_design_ctdna'},
    'catalog': catalog, 'open_items': open_items, 'pairs': pairs, 'variants': variants,
    'multiplex': multiplex, 'probes': probes, 'blocker': blocker,
    'melt': {'summary': melt_sum, 'curves': curves}, 'validation': validation, 'profiles': profiles, 'app_profiles': app_profiles, 'engine_params': engine_params, 'design_check': design_check,
}
dst = os.path.join(ROOT, 'web', 'oligoscout_data.json')
with open(dst, 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, separators=(',', ':'))
print('yazildi', dst, os.path.getsize(dst), 'bayt;', len(pairs), 'cift;', len(variants), 'varyant;',
      sum(v['ref_ok'] for v in variants), 'ref uyumlu')
