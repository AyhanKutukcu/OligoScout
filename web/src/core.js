"use strict";
/* ================= veri & yardımcılar ================= */
const D = JSON.parse(document.getElementById('osdata').textContent);
OS.setParams(D.engine_params);
const TPL = D.meta.template, OFF = D.meta.offset, CHR = D.meta.chrom;
const COND0 = {na:50, mg:1.5, dntp:0.6, pnm:200, ta:60};
const $ = (s, r=document) => r.querySelector(s);
const $$ = (s, r=document) => Array.from(r.querySelectorAll(s));
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const nf = (x, d=1) => (x==null || x==='' || isNaN(+x)) ? '—' : (+x).toLocaleString('tr-TR', {minimumFractionDigits:d, maximumFractionDigits:d});
const ni = x => (x==null||isNaN(+x)) ? '—' : Math.round(+x).toLocaleString('tr-TR');
const gpos = i => (OFF + i + 1).toLocaleString('tr-TR');          // 0-tabanlı şablon indeksi → 1-tabanlı chr22 koordinatı
const afFmt = a => a==null ? '—' : a >= 0.01 ? (a*100).toLocaleString('tr-TR',{maximumFractionDigits:1})+' %' : a >= 1e-4 ? (a*100).toLocaleString('tr-TR',{maximumFractionDigits:3})+' %' : a.toExponential(1).replace('e-','×10⁻').replace(/⁻(\d+)/, (m,g)=>'⁻'+g.split('').map(c=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[c]).join(''));
const sevRank = s => ({FAIL:3, WARN:2, INFO:1, PASS:1, OK:1}[s] || 0);
const worst = arr => arr.reduce((w, s) => sevRank(s) > sevRank(w) ? s : w, 'INFO');
const tidyKey = k => String(k).replace(/_/g,' ');

const ICON = {
  home:'<path d="M3 11l9-7 9 7v9a1 1 0 0 1-1 1h-5v-6h-6v6H4a1 1 0 0 1-1-1z"/>',
  design:'<circle cx="10.5" cy="10.5" r="6.5"/><path d="M20 20l-4.6-4.6M8 10.5h5M10.5 8v5"/>',
  pair:'<path d="M4 7h11l-3-3M20 17H9l3 3"/>',
  mpx:'<rect x="3" y="3" width="7" height="7" rx="1.5"/><rect x="14" y="3" width="7" height="7" rx="1.5"/><rect x="3" y="14" width="7" height="7" rx="1.5"/><rect x="14" y="14" width="7" height="7" rx="1.5"/>',
  probe:'<path d="M4 12h16"/><circle cx="7" cy="12" r="2.5"/><circle cx="17" cy="12" r="2.5"/>',
  melt:'<path d="M3 18c4 0 4-12 9-12s5 12 9 12"/>',
  cat:'<path d="M5 4h14v16H5z"/><path d="M9 8h6M9 12h6M9 16h4"/>',
  check:'<path d="M12 3l8 4v5c0 5-3.4 8-8 9-4.6-1-8-4-8-9V7z"/><path d="M8.5 12l2.5 2.5 4.5-5"/>',
  map:'<path d="M3 6l6-2 6 2 6-2v14l-6 2-6-2-6 2z"/><path d="M9 4v14M15 6v14"/>',
  target:'<circle cx="12" cy="12" r="8"/><circle cx="12" cy="12" r="3"/>',
  snp:'<path d="M12 3v18M6 8l6-5 6 5"/><circle cx="12" cy="15" r="2"/>',
  struct:'<path d="M7 4v16M17 4v16M7 8h10M7 12h10M7 16h10"/>',
  thermo:'<path d="M10 4a2 2 0 0 1 4 0v10a4 4 0 1 1-4 0z"/><path d="M12 9v7"/>',
  list:'<path d="M8 6h13M8 12h13M8 18h13"/><circle cx="4" cy="6" r="1"/><circle cx="4" cy="12" r="1"/><circle cx="4" cy="18" r="1"/>',
  copy:'<rect x="9" y="9" width="11" height="11" rx="2"/><path d="M5 15V5a1 1 0 0 1 1-1h10"/>',
  sun:'<circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M2 12h2M20 12h2M5 5l1.5 1.5M17.5 17.5L19 19M19 5l-1.5 1.5M6.5 17.5L5 19"/>',
  info:'<circle cx="12" cy="12" r="9"/><path d="M12 11v5M12 8h.01"/>',
  ok:'<path d="M5 12.5l4.5 4.5L19 7"/>',
  warn:'<path d="M12 3l9.5 17h-19z"/><path d="M12 10v4M12 17h.01"/>',
  fail:'<path d="M6 6l12 12M18 6L6 18"/>',
  x:'<path d="M6 6l12 12M18 6L6 18"/>',
  arrow:'<path d="M5 12h14M13 6l6 6-6 6"/>',
  dl:'<path d="M12 4v11M7 10l5 5 5-5M5 20h14"/>'
};
const ic = (n, cls='') => `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="${cls}" aria-hidden="true">${ICON[n]||''}</svg>`;
const SEVLBL = {INFO:'INFO', PASS:'PASS', OK:'OK', WARN:'WARN', FAIL:'FAIL', NA:'—', PEND:'BEKLİYOR'};
function sev(s, label) {
  const k = s || 'NA'; const icn = {INFO:'info', PASS:'ok', OK:'ok', WARN:'warn', FAIL:'fail'}[k];
  return `<span class="sev ${k}">${icn ? ic(icn) : ''}${esc(label || SEVLBL[k] || k)}</span>`;
}
const seqC = s => String(s).split('').map(c => `<span class="n${c}">${c}</span>`).join('');
const copyBtn = (txt, title='Kopyala') => `<button class="copy" data-copy="${esc(txt)}" title="${title}" aria-label="${title}">${ic('copy')}</button>`;

/* kopyalama, ipucu */
document.addEventListener('click', e => {
  const b = e.target.closest('[data-copy]'); if (!b) return;
  const t = b.getAttribute('data-copy');
  const done = () => { b.innerHTML = ic('ok'); setTimeout(() => b.innerHTML = ic('copy'), 1200); };
  try { navigator.clipboard.writeText(t).then(done, () => fallbackCopy(t, done)); } catch (_) { fallbackCopy(t, done); }
});
function fallbackCopy(t, done) { const ta = document.createElement('textarea'); ta.value = t; document.body.appendChild(ta); ta.select(); try { document.execCommand('copy'); } catch (_) {} ta.remove(); done(); }
const TIP = document.createElement('div'); TIP.className = 'tip'; TIP.hidden = true; document.body.appendChild(TIP);
document.addEventListener('mousemove', e => {
  const el = e.target.closest && e.target.closest('[data-tip]');
  if (!el) { TIP.hidden = true; return; }
  TIP.innerHTML = el.getAttribute('data-tip'); TIP.hidden = false;
  const w = TIP.offsetWidth, h = TIP.offsetHeight;
  let x = e.clientX + 14, y = e.clientY + 14;
  if (x + w > innerWidth - 8) x = e.clientX - w - 14; if (y + h > innerHeight - 8) y = e.clientY - h - 14;
  TIP.style.left = Math.max(6, x) + 'px'; TIP.style.top = Math.max(6, y) + 'px';
});
document.addEventListener('scroll', () => TIP.hidden = true, true);

/* tema */
(function () {
  const root = document.documentElement; let s = null; try { s = localStorage.getItem('os-theme'); } catch (_) {}
  if (s) root.setAttribute('data-theme', s);
  document.addEventListener('click', e => {
    if (!e.target.closest('#themeBtn')) return;
    const cur = root.getAttribute('data-theme') || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    const nx = cur === 'dark' ? 'light' : 'dark'; root.setAttribute('data-theme', nx); try { localStorage.setItem('os-theme', nx); } catch (_) {}
    rerender();
  });
})();
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();

/* ================= katalog ================= */
const CAT = D.catalog.map(r => {
  const parts = (r.title || '').split(' / ');
  return {...r, scope: parts[0] || '', name: parts[1] || (r.id === 'REP-03' ? 'Kalibrasyon durumu beyanı' : r.id === 'VAR-08' ? 'Varyant özeti (çift düzeyi)' : r.id), method: parts.slice(2).join(' / ')};
});
const CATI = Object.fromEntries(CAT.map(r => [r.id, r]));
const FAMORDER = ['Dizi butunlugu','Kompozisyon','Termodinamik','Ikincil yapi','Dimer','Cift uyumu','Ozgulluk','Amplikon','Prob','ctDNA','mRNA/RT','qPCR','Bisulfit','Multipleks','Varyant','Kosul','NGS','Sanger','Rapor'];
const FAMTR = {'Dizi butunlugu':'Dizi bütünlüğü','Kompozisyon':'Kompozisyon','Termodinamik':'Termodinamik','Ikincil yapi':'İkincil yapı','Dimer':'Dimer','Cift uyumu':'Çift uyumu','Ozgulluk':'Özgüllük','Amplikon':'Amplikon','Prob':'Prob','ctDNA':'ctDNA','mRNA/RT':'mRNA / RT','qPCR':'qPCR','Bisulfit':'Bisülfit','Multipleks':'Multipleks','Varyant':'Varyant','Kosul':'Reaksiyon koşulu','NGS':'NGS','Sanger':'Sanger','Rapor':'Rapor'};
const STATUS = {
  etiketli_uygulandi: {t:'Uygulandı', c:'PASS'}, etiketli_kismi: {t:'Kısmi', c:'WARN'},
  etiketli_kismi_tablo_pending: {t:'Kısmi · tablo bekliyor', c:'WARN'}, wsl_script_hazir: {t:'WSL indeksi bekliyor', c:'PEND'},
  veri_gerekir_transkriptom_index: {t:'Transkriptom indeksi gerekir', c:'PEND'}
};
const isImpl = r => r.status === 'etiketli_uygulandi' || r.status === 'etiketli_kismi';
const CAT_IMPL = CAT.filter(isImpl).length;            // 86 = 85 tam + SPC-08 kısmi
const CAT_FULL = CAT.filter(r => r.status === 'etiketli_uygulandi').length;           // 85

/* ================= çift modeli ================= */
const PROFILES = Object.fromEntries((D.app_profiles || []).map(p => [p.profile, p]));
const PROJ = D.pairs.map(p => ({
  key: 'P:' + p.pair_id, kind: 'project', id: p.pair_id, label: p.pair_id, rank: p.final_rank,
  F: p.left_primer, R: p.right_primer, fStart: p.fStart, fEnd: p.fStart + p.left_primer.length - 1,
  rStart: p.rStart, rEnd: p.rEnd, prod: p.product_length, tmF: p.tm_left_c, tmR: p.tm_right_c, dtm: p.delta_tm_c,
  thetaF: p.theta_f_at_ta, thetaR: p.theta_r_at_ta, verdict: p.verdict, spec: p.specificity_verdict, varv: p.variant_verdict,
  score: p.thermo_score, tpl: TPL, tplOff: OFF, chrom: CHR, hasVar: true, profile: 'ctdna', cond: {...COND0}, raw: p
}));
const STATE = { view: 'genel', pair: PROJ[0].key, ptab: 'harita', local: [], localTpl: null, design: null, mpxSel: null };
const getPair = k => PROJ.find(p => p.key === k) || STATE.local.find(p => p.key === k) || PROJ[0];

/* varyantlar (gnomAD v4 önbelleği, şablon indeksli) */
const VARS = (() => { const s = new Set(); return D.variants.filter(v => { const k = v.pos + v.ref + '>' + v.alt; if (s.has(k)) return false; s.add(k); return true; }).sort((a, b) => a.idx - b.idx); })();
const varClass = af => af >= 0.01 ? 'common' : af >= 0.001 ? 'low' : 'rare';
const varColor = af => af >= 0.01 ? css('--fail') : af >= 0.001 ? css('--warn') : css('--ink3');
function varsIn(lo, hi) { return VARS.filter(v => v.idx <= hi && v.idx + Math.max(1, v.ref.length) - 1 >= lo); }
function vtype(v) { return v.ref.length === 1 && v.alt.length === 1 ? 'SNV' : v.ref.length > v.alt.length ? 'Delesyon' : 'İnsersiyon'; }

/* ================= kanıt (91 kural) ================= */
// her kural için: {status:'eval'|'na'|'pend'|'panel'|'profile', sev, recs:[{role,src,value,thr,note,model}], reason}
const PROFILE_FAMS = {ctdna: ['ctDNA'], ddpcr: ['Prob'], rt_qpcr: ['mRNA/RT'], qpcr_sybr: ['qPCR'], sanger: ['Sanger'], ngs_amplicon: ['NGS'], bisulfite: ['Bisulfit'], standard_pcr: []};
const PROFILE_SPECIFIC = ['ctDNA','Prob','mRNA/RT','qPCR','Sanger','NGS','Bisulfit'];
const ROLEN = {forward:'F', reverse:'R'};
function rec(role, src, c, extra={}) { role = ROLEN[role] || role; return {role, src, id: c.id, sev: c.sev ?? c.severity, value: c.value, thr: c.thr ?? c.threshold, note: c.note, model: c.model, ...extra}; }

function evidenceFor(pm) {
  const E = {}; const add = (id, r) => { (E[id] = E[id] || {status:'eval', recs:[]}).recs.push(r); };
  const setNA = (id, status, reason) => { if (!E[id]) E[id] = {status, recs:[], reason}; };
  const browser = (id, role, sev, value, thr, note, model) => add(id, {role, src:'tarayıcı', id, sev, value, thr, note, model});
  const prof = pm.profile || 'ctdna';
  if (pm.kind === 'project') {
    const p = pm.raw;
    for (const c of p.seqpair) add(c.id, rec(c.role, 'sequence_qc.py', c));
    for (const role of ['F','R']) { const pr = p.primers[role]; if (!pr) continue;
      for (const c of pr.checks) add(c.id, rec(role, 'panel_qc.py', c));
      for (const [id, why] of Object.entries(pr.not_evaluated || {})) if (!E[id] || !E[id].recs.some(r => r.role === role)) (E[id] = E[id] || {status:'eval', recs:[]}).recs.push({role, src:'panel_qc.py', id, sev:null, value:null, note:'değerlendirilmedi: ' + why, ne:true});
    }
    for (const c of p.pair_checks || []) add(c.id, rec('pair', 'panel_qc.py', c));
    for (const c of D.meta.panel_checks || []) add(c.id, rec('panel', 'panel_qc.py', c, {panelLevel:true}));
    if (p.variants) { for (const role of ['F','R']) for (const c of p.variants[role]) add(c.id, rec(role, 'variant_qc.py', c)); for (const c of p.variants.pair || []) add(c.id, rec('amplicon', 'variant_qc.py', c)); }
    const h = p.hetero_current;
    if (h) {
      const a = h.any, t = h['3prime'];
      const dA = a && a.found !== false ? a.dG37_kcal_mol : 0, d3 = t && t.found !== false ? t.dG37_kcal_mol : 0;
      add('DIM-03', {role:'pair', src:'align.best_dimer (veri paketi)', id:'DIM-03', sev: dA < -9 ? 'FAIL' : dA < -6 ? 'WARN' : 'INFO', value: a, thr:{warn:'dG37 < -6.0', fail:'dG37 < -9.0'}, note:'F×R hetero-dimer; eşikler DIM-01 ile aynı'});
      add('DIM-04', {role:'pair', src:'align.best_dimer (veri paketi)', id:'DIM-04', sev: d3 < -7 ? 'FAIL' : d3 < -5 ? 'WARN' : 'INFO', value: t, thr:{warn:'dG37 < -5.0', fail:'dG37 < -7.0'}, note:"3'-uçlu, uzatılabilir F×R dimeri; eşikler DIM-02 ile aynı"});
      const rev = OS.bestDimer(pm.R, pm.F, 'any'); const fw = OS.bestDimer(pm.F, pm.R, 'any');
      browser('DIM-08', 'pair', Math.abs((rev?rev.dG:0) - (fw?fw.dG:0)) < 1e-9 ? 'INFO' : 'FAIL', {dimer_FR: +(fw?fw.dG:0).toFixed(4), dimer_RF: +(rev?rev.dG:0).toFixed(4)}, {fail:'dimer(A,B) ≠ dimer(B,A)'}, 'Simetri değişmezliği tarayıcı portunda canlı sınandı', {model:'align.py portu'});
    }
    // ΔTm ve TH-07 (motor değerlerinden türetim)
    browser('PAIR-01', 'pair', pm.dtm > 2 ? 'WARN' : 'INFO', {delta_tm_c: pm.dtm, tm_left_c: pm.tmF, tm_right_c: pm.tmR}, {warn:'> 2.0 °C (DesignConstraints.max_delta_tm)'}, 'Değer final_ranked.tsv; eşik design.py sert filtresi', null);
    const p04 = p.seqpair.find(c => c.id === 'PAIR-04');
    if (p04) { const tp = p04.value.product_tm_c, tmin = Math.min(pm.tmF, pm.tmR); const ta = 0.3*tmin + 0.7*tp - 14.9;
      browser('TH-07', 'pair', Math.abs(ta - pm.cond.ta) > 5 ? 'WARN' : 'INFO', {rychlik_ta_c: +ta.toFixed(2), anneal_c: pm.cond.ta, product_tm_c: tp, tm_primer_min_c: tmin}, {warn:'|Ta(Rychlik) − Ta| > 5 °C'}, 'Rychlik 1990 formülü motorun PAIR-04 ürün Tm değeriyle tarayıcıda uygulandı', null); }
    // CT ailesi
    if (p.detect && p.detect.length) { const d0 = p.detect[0];
      add('CT-02', {role:'amplicon', src:'design_primers.py (ctdna_detectability)', id:'CT-02', sev:'INFO', value:{amplicon_bp:d0.amplicon_bp, p_intact_fragment:d0.p_intact_fragment}, thr:{}, note:'cfDNA fragman dağılımı üzerinden bütün fragman olasılığı'});
      add('CT-03', {role:'pair', src:'design_primers.py (ctdna_detectability)', id:'CT-03', sev: p.detect.find(x => x.vaf === 0.001)?.p_detect < 0.5 ? 'WARN' : 'INFO', value:Object.fromEntries(p.detect.map(x => ['VAF '+x.vaf, x.p_detect])), thr:{warn:'P(tespit | VAF 0,1 %) < 0,5'}, note:'λ = GE × VAF × P_bütün; GE ' + d0.genome_equivalents});
    }
    browser('CT-01', 'amplicon', pm.prod > 100 ? 'FAIL' : 'INFO', {amplicon_bp: pm.prod}, {fail:'> 100 bp (ctdna profili)'}, 'ctDNA profili: cfDNA ~166 bp tepe', null);
    setNA('CT-04', 'na', 'ARMS alel-spesifik primer tasarlanmadı; 3′ terminal uyumsuzluk tablosu kaynak farkı (GG/CG) kullanıcı kararı bekliyor');
    setNA('CT-05', 'panel', 'Proje reaktif seti (blocker/prob) için ayrı koşu — Prob & Blocker sayfası');
    setNA('CT-06', 'na', 'Hedef gen beyan edilmedi (sequence_qc.py --gene)');
    // multipleks panel düzeyi
    for (const m of D.multiplex.summary) { const id = m.check === 'MPX-02f' ? 'MPX-02' : m.check;
      add(id, {role:'panel', src:'multiplex_qc.py', id: m.check, sev: m.severity, value: m.value, thr: m.threshold, note: m.note, panelLevel:true}); }
  } else {
    localChecks(pm, add, browser);
    for (const id of ['SPC-01','SPC-02','SPC-03','SPC-04','SPC-05','SPC-06','SPC-07','SPC-09']) setNA(id, 'scan', 'Genom-geneli tarama gerekir (C++ oligoscout_panel → panel_qc.py)');
    for (const id of ['MPX-01','MPX-02','MPX-03','DIM-06']) setNA(id, 'panel', 'Panel düzeyi; Multipleks sayfasında havuz olarak değerlendirilir');
    if (prof === 'ctdna') browser('CT-01', 'amplicon', pm.prod > 100 ? 'FAIL' : 'INFO', {amplicon_bp: pm.prod}, {fail:'> 100 bp'}, 'ctDNA profili', null);
  }
  // ortak, beyan edilen kayıtlar
  if (!E['SEQ-01']) for (const [role, s] of [['F', pm.F], ['R', pm.R]]) browser('SEQ-01', role, /^[ACGT]+$/.test(s) ? 'INFO' : 'FAIL', {pure_acgt: /^[ACGT]+$/.test(s), length: s.length}, {fail:'A/C/G/T dışı tanınmayan karakter'}, 'oligomod ayrıştırması: modifikasyon yok', null);
  if (!E['CMP-01']) for (const [role, s] of [['F', pm.F], ['R', pm.R]]) { const g = OS.gcf(s), cl = OS.clamp(s), ho = OS.longestHomo(s), rp = Math.max(OS.repUnits(s,2), OS.repUnits(s,3));
    browser('CMP-01', role, (g < .4 || g > .6) ? 'FAIL' : 'INFO', {gc_pct: +(100*g).toFixed(2)}, {fail:'%40–%60 dışı (design.py sert filtre)'}, 'Tasarım sert filtresi; tarayıcıda yeniden hesaplandı', null);
    browser('CMP-02', role, (cl < 1 || cl > 3) ? 'FAIL' : 'INFO', {gc_in_last5: cl}, {fail:'1–3 dışı'}, '', null);
    browser('CMP-03', role, ho >= 4 || s.includes('GGGG') ? 'FAIL' : 'INFO', {longest_homopolymer: ho}, {fail:'≥ 4 veya GGGG'}, '', null);
    browser('CMP-04', role, rp >= 4 ? 'FAIL' : 'INFO', {max_repeat_units: rp}, {fail:'≥ 4 ardışık birim'}, '', null); }
  if (!E['SEQ-06']) setNA('SEQ-06', 'eval0', null);
  if (pm.kind === 'project') { E['SEQ-06'] = {status:'eval', recs:[{role:'template', src:'veri paketi', id:'SEQ-06', sev:'INFO', value:{header: D.meta.header, build:'GRCh38', patch:'p14', bases: TPL.length}, thr:{warn:'FASTA SHA-256 hesaplanamadı'}, note:'FASTA SHA-256 sequence_qc.py çıktısında; bu pakete yalnız başlık alındı'}]}; }
  else E['SEQ-06'] = {status:'eval', recs:[{role:'template', src:'tarayıcı', id:'SEQ-06', sev: pm.tplName ? 'INFO' : 'WARN', value:{header: pm.tplName || null, bases: pm.tpl.length}, thr:{warn:'referans kimliği beyan edilmedi'}, note:'Yapıştırılan şablonun başlığı'}]};
  // rapor
  E['REP-01'] = {status:'eval', recs:[{role:'panel', src:'veri paketi', id:'REP-01', sev:'INFO', value:{param_set_sha256: D.meta.seqpair_param_sha.slice(0,16)+'…', registry_digest: (D.engine_params.registry_digest||'').slice(0,16)+'…', same: D.meta.seqpair_param_sha === D.engine_params.registry_digest}, thr:{}, note:'Kayıtların parametre digest’i ile params.registry_digest() aynı'}]};
  E['REP-02'] = {status:'eval', recs:[{role:'panel', src:'veri paketi', id:'REP-02', sev:'INFO', value:{alanlar:'id, value, threshold, severity, model, conditions, param_set_sha256, note'}, thr:{}, note:'Her kayıt kanıt şemasında'}]};
  E['REP-03'] = {status:'eval', recs:[{role:'panel', src:'veri paketi', id:'REP-03', sev:'WARN', value:{calibration_status:'REP-03 / KALİBRE DEĞİL', parameter_tables_verified:false, ranking_score_calibrated:false}, thr:{}, note:'Eşikler literatür başlangıç değeri; ıslak-lab kalibrasyonu bekliyor (22 yer)'}]};
  // kalan kurallar: durum ataması
  for (const r of CAT) {
    if (E[r.id]) { const e = E[r.id]; if (e.recs.length) { const s = e.recs.filter(x => !x.ne).map(x => x.sev).filter(Boolean); e.sev = s.length ? worst(s) : null; if (!s.length) { e.status = 'ne'; e.reason = e.recs[0].note; } } continue; }
    if (!isImpl(r)) { E[r.id] = {status:'pend', recs:[], reason: STATUS[r.status]?.t + ' · ' + r.where}; continue; }
    if (PROFILE_SPECIFIC.includes(r.fam) && !(PROFILE_FAMS[prof] || []).includes(r.fam)) { E[r.id] = {status:'profile', recs:[], reason:`${FAMTR[r.fam]} ailesi ${prof} profilinde uygulanmaz`}; continue; }
    E[r.id] = {status:'ne', recs:[], reason: NE_REASON[r.id] || 'Bu çift için koşulmadı'};
  }
  return E;
}
const NE_REASON = {
  'SEQ-03':'Dejenere baz yok', 'DIM-05':'Bu çift için prob tanımlanmadı (Prob & Blocker sayfası)', 'DIM-07':'Adaptör/M13 kuyruğu tanımlanmadı (sequence_qc.py --tails)',
  'SPC-08':'--contigs verilmedi (chrM Python; ALT/rDNA C++ bekliyor)', 'TH-06':'Katkı maddesi beyan edilmedi', 'AMP-02':'Erime profili qPCR koşusunda (qPCR Erime sayfası)',
  'VAR-05':'Varyant-farkında indeks gerekir', 'SS-05':'Prob yok'
};

/* ---- yerel (tarayıcıda tasarlanan) çift için kontroller ---- */
function localChecks(pm, add, browser) {
  const c = pm.cond;
  for (const [role, s] of [['F', pm.F], ['R', pm.R]]) {
    const t = OS.thermo(s, c, pm.salt || 'owczarzy_2008'), sp = OS.saltSpread(s, c), d3 = OS.dG3(s), hv = OS.halves(s);
    const pr = PROFILES[pm.profile] || PROFILES.ctdna; const [tlo, thi] = pr ? pr.tm_window_c : [58, 62];
    browser('SEQ-02', role, (s.length < 15 || s.length > 35) ? 'FAIL' : (s.length < 18 || s.length > 30) ? 'WARN' : 'INFO', {length_nt: s.length}, {warn:'18-30 dışı', fail:'15-35 dışı'}, '', null);
    browser('SEQ-03', role, 'INFO', {degeneracy: 1}, {warn:'>32', fail:'>256'}, 'dejenere baz yok', null);
    const e1 = OS.shannon(s, 1), e2 = OS.shannon(s, 2);
    browser('CMP-05', role, e1 < 1.2 ? 'FAIL' : e1 < 1.5 ? 'WARN' : 'INFO', {shannon_k1_bit: +e1.toFixed(4), shannon_k2_bit: +e2.toFixed(4)}, {warn:'k1 < 1.5 bit', fail:'k1 < 1.2 bit'}, '', null);
    browser('CMP-06', role, OS.g4(s) ? 'FAIL' : 'INFO', {g_quadruplex: OS.g4(s)}, {fail:'G≥3 N1-7 dört tekrar'}, '', null);
    browser('CMP-07', role, (hv.d3 - hv.d5) < -2 ? 'WARN' : 'INFO', {dG5_half_kcal_mol: +hv.d5.toFixed(3), dG3_half_kcal_mol: +hv.d3.toFixed(3), delta_kcal_mol: +(hv.d3 - hv.d5).toFixed(3)}, {warn:"3' yarı > 2 kcal/mol daha kararlı"}, '', null);
    browser('CMP-08', role, s.endsWith('T') ? 'WARN' : 'INFO', {three_prime_base: s.slice(-1)}, {warn:"3'-T"}, '', null);
    browser('TH-01', role, (t.tm < tlo - 3 || t.tm > thi + 3) ? 'FAIL' : (t.tm < tlo || t.tm > thi) ? 'WARN' : 'INFO', {dH_kcal_mol: +t.dH.toFixed(2), dS_cal_mol_K: +t.dS.toFixed(2), tm_1M_c: +t.tm1M.toFixed(3), tm_c: +t.tm.toFixed(3)}, {warn:`${tlo}-${thi} °C dışı`, fail:`${tlo-3}-${thi+3} °C dışı`}, 'concentration.build(primer_excess) + ' + (pm.salt || 'owczarzy_2008'), {model:'motor portu'});
    browser('TH-02', role, sp.spread > 4 ? 'FAIL' : sp.spread > 2 ? 'WARN' : 'INFO', {owczarzy_2004: +sp.owczarzy_2004.toFixed(3), owczarzy_2008: +sp.owczarzy_2008.toFixed(3), von_ahsen_2001: +sp.von_ahsen_2001.toFixed(3), spread_c: +sp.spread.toFixed(3)}, {warn:'>2.0 C', fail:'>4.0 C'}, '', null);
    browser('TH-03', role, t.thetaTa < 0.8 ? 'FAIL' : t.thetaTa < 0.95 ? 'WARN' : 'INFO', {theta_at_ta: +t.thetaTa.toFixed(6), dG_ta_kcal_mol: +t.dGta.toFixed(4), ta_c: c.ta}, {warn:'<0.95', fail:'<0.8'}, '', null);
    browser('TH-04', role, d3 < -11 ? 'FAIL' : (d3 > -6 || d3 < -9) ? 'WARN' : 'INFO', {dG37_last_5_nt_kcal_mol: +d3.toFixed(4)}, {warn:'>-6.0 veya <-9.0', fail:'<-11.0'}, '', null);
    browser('TH-05', role, Math.abs(t.theta(t.tm) - 0.5) > 1e-9 ? 'FAIL' : 'INFO', {theta_at_tm: +t.theta(t.tm).toFixed(12)}, {fail:'|θ(Tm) − 0.5| > 1e-9'}, 'analitik özdeşlik testi', null);
    const da = OS.bestDimer(s, s, 'any'), d3p = OS.bestDimer(s, s, '3prime');
    browser('DIM-01', role, da && da.dG < -9 ? 'FAIL' : da && da.dG < -6 ? 'WARN' : 'INFO', da ? {...dimToVal(da)} : {found:false}, {warn:'dG37 < -6.0', fail:'dG37 < -9.0'}, '', {model:'nn_single_internal_mismatch_or_bulge'});
    browser('DIM-02', role, d3p && d3p.dG < -7 ? 'FAIL' : d3p && d3p.dG < -5 ? 'WARN' : 'INFO', d3p ? {...dimToVal(d3p)} : {found:false}, {warn:'dG37 < -5.0', fail:'dG37 < -7.0'}, '', {model:'nn_single_internal_mismatch_or_bulge'});
    const hp = OS.hairpinStem(s);
    browser('SS-02', role, hp && hp.three_prime_paired && hp.dG < -2 ? 'WARN' : 'INFO', {stem_dG37_kcal_mol: hp ? +hp.dG.toFixed(3) : 0, three_prime_extendable: !!(hp && hp.three_prime_paired)}, {warn:"3' uzatılabilir ve dG37 < -2.0"}, 'align.hairpin_stem (yalnız sap yığını)', null);
    E_SS01(browser, role, hp);
    if (pm.hasVar) localVar(pm, role, browser);
  }
  const fw = OS.bestDimer(pm.F, pm.R, 'any'), f3 = OS.bestDimer(pm.F, pm.R, '3prime'), rv = OS.bestDimer(pm.R, pm.F, 'any');
  browser('DIM-03', 'pair', fw && fw.dG < -9 ? 'FAIL' : fw && fw.dG < -6 ? 'WARN' : 'INFO', fw ? dimToVal(fw) : {found:false}, {warn:'dG37 < -6.0', fail:'dG37 < -9.0'}, '', null);
  browser('DIM-04', 'pair', f3 && f3.dG < -7 ? 'FAIL' : f3 && f3.dG < -5 ? 'WARN' : 'INFO', f3 ? dimToVal(f3) : {found:false}, {warn:'dG37 < -5.0', fail:'dG37 < -7.0'}, '', null);
  browser('DIM-08', 'pair', Math.abs((fw?fw.dG:0) - (rv?rv.dG:0)) < 1e-9 ? 'INFO' : 'FAIL', {dimer_FR: fw ? +fw.dG.toFixed(4) : 0, dimer_RF: rv ? +rv.dG.toFixed(4) : 0}, {fail:'asimetri'}, '', null);
  browser('PAIR-01', 'pair', pm.dtm > 2 ? 'WARN' : 'INFO', {delta_tm_c: +pm.dtm.toFixed(3)}, {warn:'> 2.0 °C'}, '', null);
  const pr = PROFILES[pm.profile] || {product_bp:[60,100]}; const [plo, phi] = pr.product_bp; const edge = (phi - plo) * 0.1;
  browser('PAIR-02', 'pair', (pm.prod < plo || pm.prod > phi) ? 'FAIL' : (pm.prod < plo + edge || pm.prod > phi - edge) ? 'WARN' : 'INFO', {product_bp: pm.prod}, {warn:'pencere kenarına < %10', fail:`${plo}-${phi} bp dışı`}, '', null);
  const gap = pm.rStart - pm.fEnd - 1;
  browser('PAIR-03', 'pair', gap < 0 ? 'FAIL' : gap < 10 ? 'WARN' : 'INFO', {inner_gap_nt: gap}, {warn:'ara < 10 nt', fail:'örtüşme'}, '', null);
  const amp = pm.tpl.slice(pm.fStart, pm.rEnd + 1); const gg = OS.gcf(amp), gw = OS.gcWindowMax(amp);
  browser('AMP-01', 'amplicon', gw > .85 ? 'FAIL' : (gg < .35 || gg > .65 || gw > .75) ? 'WARN' : 'INFO', {gc_global: +gg.toFixed(4), gc_window_max: +gw.toFixed(4), window_nt: 50}, {warn:'global %35-%65 dışı veya pencere > %75', fail:'pencere > %85'}, '', null);
  browser('AMP-03', 'amplicon', OS.shannon(amp,1) < 1.5 || OS.g4(amp) ? 'WARN' : 'INFO', {amplicon_shannon_k1_bit: +OS.shannon(amp,1).toFixed(4), amplicon_g_quadruplex: OS.g4(amp)}, {warn:'düşük karmaşıklık ya da G-quadruplex'}, '', null);
  const nN = (amp.match(/N/g) || []).length;
  browser('SEQ-04', 'template', nN ? 'WARN' : 'INFO', {amplicon_N: nN}, {warn:'amplikonda N > 0', fail:'primer bölgesinde N'}, '', null);
  browser('SEQ-05', 'template', 'INFO', {forward_matches_template: pm.tpl.slice(pm.fStart, pm.fEnd + 1) === pm.F, reverse_matches_template: OS.rc(pm.tpl.slice(pm.rStart, pm.rEnd + 1)) === pm.R, product_start: pm.fStart, product_end: pm.rEnd + 1}, {fail:'koordinattaki dizi primerle uyuşmuyor'}, '', null);
  add('TH-07', {role:'pair', src:'tarayıcı', id:'TH-07', sev:null, ne:true, note:'Ürün Tm modeli (amplicon.py) tarayıcıya taşınmadı; RUN_MELT.cmd ile hesaplanır'});
}
function E_SS01(browser, role, hp) { browser('SS-01', role, hp && hp.dG < -3 ? 'WARN' : 'INFO', hp ? {stem_dG37_kcal_mol: +hp.dG.toFixed(3), dot_bracket: hp.db, loop_nt: hp.loop} : {dot_bracket:'yok'}, {warn:'sap dG37 < -3.0 (tarayıcı göstergesi)'}, 'Tarayıcı: hairpin_stem yaklaşık; motor SS-01 Zuker MFE + döngü terimleri kullanır', null); }
function dimToVal(d) { return {alignment: d.alignment, dG37_kcal_mol: +d.dG.toFixed(4), dH_kcal_mol: +d.dH.toFixed(2), dS_cal_mol_K: +d.dS.toFixed(2), helix_len: d.helix_len, three_prime_paired: d.three_prime_paired, mismatch_count: d.mismatch_count, bulge_count: d.bulge_count}; }
function footprint(pm, role) { return role === 'F' ? [pm.fStart, pm.fEnd] : [pm.rStart, pm.rEnd]; }
function dist3(pm, role, idx) { return role === 'F' ? pm.fEnd - idx : idx - pm.rStart; }
function localVar(pm, role, browser) {
  const [lo, hi] = footprint(pm, role); const vs = varsIn(lo, hi);
  const max = vs.reduce((m, v) => Math.max(m, v.af || 0), 0);
  const w5 = vs.filter(v => dist3(pm, role, v.idx) < 5); const m5 = w5.reduce((m, v) => Math.max(m, v.af || 0), 0);
  browser('VAR-01', role, max >= 0.001 ? 'WARN' : 'INFO', {n_variants: vs.length, max_af: max}, {warn:'AF >= 0.001'}, 'bağlanma yerinin tamamı tarandı (gnomAD v4 önbelleği)', null);
  browser('VAR-02', role, m5 >= 0.01 ? 'FAIL' : m5 >= 1e-4 ? 'WARN' : 'INFO', {window_nt: 5, n_variants: w5.length, max_af: m5}, {warn:'AF >= 0.0001', fail:'AF >= 0.01'}, "3' uca yakın uyumsuzluk uzamayı doğrudan engeller", null);
  const cm = vs.filter(v => v.af >= 0.01);
  browser('VAR-03', role, max >= 0.05 ? 'FAIL' : max >= 0.01 ? 'WARN' : 'INFO', {n_common: cm.length, max_af: max}, {warn:'AF >= 0.01', fail:'AF >= 0.05'}, 'yaygın varyant kohortun büyük kısmında alel düşmesi demek', null);
}

/* ================= yönlendirme ================= */
const VIEWS = [
  {id:'genel', t:'Genel Bakış', i:'home', grp:'Proje'},
  {id:'tasarla', t:'Primer Tasarla', i:'design', grp:'Proje'},
  {id:'cift', t:'Çift İnceleme', i:'pair', grp:'Proje'},
  {id:'multipleks', t:'Multipleks', i:'mpx', grp:'Panel'},
  {id:'prob', t:'Prob & Blocker', i:'probe', grp:'Panel'},
  {id:'erime', t:'qPCR Erime', i:'melt', grp:'Panel'},
  {id:'katalog', t:'QC Kataloğu', i:'cat', grp:'Motor'},
  {id:'dogrulama', t:'Doğrulama', i:'check', grp:'Motor'}
];
const RENDER = {};
function route() {
  const h = (location.hash || '').replace('#', '');
  STATE.view = VIEWS.some(v => v.id === h) ? h : 'genel';
  rerender(true);
}
function go(view, opts={}) { Object.assign(STATE, opts); if (location.hash !== '#' + view) location.hash = view; else rerender(true); }
function rerender(scrollTop) {
  $$('.nav a').forEach(a => a.setAttribute('aria-current', a.dataset.v === STATE.view ? 'page' : 'false'));
  const v = VIEWS.find(x => x.id === STATE.view);
  $('#vtitle').textContent = v.t; $('#vsub').textContent = SUBT[v.id] || '';
  const el = $('#view'); el.innerHTML = ''; el.className = 'wrap fadein';
  try { RENDER[v.id](el); } catch (err) { el.innerHTML = `<div class="note warn">${ic('warn')}<div><b>Görünüm çizilemedi.</b> ${esc(err.message)}</div></div>`; console.error(err); }
  if (scrollTop) window.scrollTo(0, 0);
}
const SUBT = {
  genel:`${CHR}:${(OFF+1).toLocaleString('tr-TR')}–${(OFF+TPL.length).toLocaleString('tr-TR')} · GRCh38.p14 · ctDNA profili`,
  tasarla:'Primer-BLAST akışı · referans motorun birebir tarayıcı portu',
  cift:'Harita, özgüllük, varyant, yapı, termodinamik ve 91 kuralın kanıtı',
  multipleks:'Çiftler arası dimer matrisi, çapraz ürünler ve havuz kurucu',
  prob:'ddPCR hidroliz probları ve CT-05 alel ayrım reaktifleri',
  erime:'Amplikon erime profili · Poland tipi model · QPCR-02',
  katalog:`19 aile · 91 kural · ${CAT_IMPL} uygulanmış`,
  dogrulama:'Tarayıcı ↔ Python ↔ C++ ↔ primer3 denetimleri'
};
window.addEventListener('hashchange', route);
