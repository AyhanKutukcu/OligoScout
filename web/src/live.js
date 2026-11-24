/* Local Python and C++ operations. The original report also opens as a file. */
const LIVE = {
  status: null,
  mode: location.protocol === 'file:' ? 'file' : 'connecting',
  pending: new Map(),
  errors: new Map(),
  tool: 'scan',
  source: '',
  last: null,
};

function liveReady() { return LIVE.mode === 'ready'; }

function analysisStatus(job) {
  if (job.status === 'queued') return 'Sırada';
  return job.status === 'running' ? 'Analiz sürüyor' : 'Analiz hazırlanıyor';
}

function analysisLoading(job) {
  return `<div class="analysis-loading"><span class="analysis-spinner" aria-hidden="true"></span><span>${analysisStatus(job)}</span></div>`;
}

function refreshConnection() {
  const badge = $('#liveConnection');
  if (!badge) return;
  const job = [...LIVE.pending.values()].at(-1);
  const label = job ? analysisStatus(job) : '';
  badge.hidden = !label;
  const text = $('#liveStatusText', badge);
  if (text.textContent !== label) text.textContent = label;
}

async function liveRequest(path, body) {
  const options = body ? {method: 'POST', headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(body)} : {};
  const response = await fetch(path, {...options, cache: 'no-store'});
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || 'Yerel motor isteği tamamlanamadı.');
  return data;
}

async function connectLocal() {
  if (LIVE.mode === 'file') return;
  try {
    LIVE.status = await liveRequest('/api/status');
    LIVE.mode = 'ready';
  } catch {
    LIVE.mode = 'unavailable';
  }
  refreshConnection();
  if (['genel', 'tasarla', 'cift', 'yerel', 'multipleks'].includes(STATE.view)) rerender();
}

function liveNote(text, warning = false) {
  return `<div class="note ${warning ? 'warn' : ''}">${ic(warning ? 'warn' : 'info')}<div>${esc(text)}</div></div>`;
}

function localAvailability() {
  return 'Analiz başlatmak için RUN_WEB.cmd dosyasını çalıştırıp yerel siteyi açın.';
}

function jobNotice(key) {
  const pending = LIVE.pending.get(key);
  const error = LIVE.errors.get(key);
  return `<div class="live-notice" data-live-progress="${esc(key)}" aria-live="polite">${pending
    ? analysisLoading(pending) : error ? liveNote(error, true) : ''}</div>`;
}

function updateJobNotice(key) {
  $$('[data-live-progress]').filter(el => el.dataset.liveProgress === key).forEach(el => {
    const state = LIVE.pending.get(key);
    const error = LIVE.errors.get(key);
    el.innerHTML = state ? analysisLoading(state) : error ? liveNote(error, true) : '';
  });
  refreshConnection();
}

async function runLocalJob(key, data, onDone) {
  if (!liveReady() || LIVE.pending.has(key)) return;
  LIVE.pending.set(key, {status:'starting'});
  LIVE.errors.delete(key);
  updateJobNotice(key);
  try {
    let job = await liveRequest('/api/jobs', data);
    while (job.status === 'queued' || job.status === 'running') {
      const statusChanged = LIVE.pending.get(key)?.status !== job.status;
      LIVE.pending.set(key, job);
      if (statusChanged) updateJobNotice(key);
      await new Promise(resolve => setTimeout(resolve, job.status === 'queued' ? 900 : 550));
      job = await liveRequest('/api/jobs/' + job.id);
    }
    if (job.status !== 'done') throw new Error(job.error || 'Motor sonuç üretemedi.');
    LIVE.pending.delete(key);
    onDone(job.result, job.id);
  } catch (error) {
    LIVE.errors.set(key, error.message);
  } finally {
    LIVE.pending.delete(key);
    updateJobNotice(key);
  }
}

function readLiveRegion(value) {
  const clean = String(value).trim().replace(/[.,]/g, '');
  const parts = clean.match(/^(?:chr)?([1-9]|1[0-9]|2[0-2]|X|Y):([0-9]+)-([0-9]+)$/i);
  if (!parts) throw new Error('Bölgeyi kromozom:başlangıç-bitiş biçiminde girin.');
  const first = Number(parts[2]), last = Number(parts[3]), size = last - first + 1;
  if (!Number.isSafeInteger(last) || first < 1 || size < 1 || size > 5000) {
    throw new Error('Bölge 1 tabanlı ve 1–5.000 baz uzunluğunda olmalı.');
  }
  return {region: 'chr' + parts[1].toUpperCase() + ':' + first + '-' + last, size};
}

function liveDesignTemplate() {
  if (DSG.src === 'bolge') return readLiveRegion(DSG.region);
  const fasta = DSG.fasta.trim();
  const lines = fasta.split(/\r?\n/);
  const header = lines[0].match(/^>\s*(\S+)(?:\s|$)/);
  if (!header || lines.slice(1).some(line => line.startsWith('>'))) {
    throw new Error('Tek FASTA kaydı ve başında GRCh38 koordinatı gerekli.');
  }
  const template = readLiveRegion(header[1]);
  const sequence = lines.slice(1).join('').replace(/\s/g, '').toUpperCase();
  if (!/^[ACGTN]+$/.test(sequence) || sequence.length !== template.size) {
    throw new Error('FASTA dizisinin uzunluğu ve alfabesi başlıktaki bölgeyle uyumlu olmalı.');
  }
  return {...template, sequence};
}

function liveConditions(form) {
  const value = name => form.elements[name].valueAsNumber;
  return {na: value('na'), mg: value('mg'), dntp: value('dntp'), pnm: value('pnm'), ta: value('ta')};
}

function submitLiveDesign(el) {
  const form = $('#dform', el);
  const button = $('#goBtn', el);
  try {
    const t = liveDesignTemplate();
    const read = id => $('#' + id, el).valueAsNumber;
    const range = (first, last, requireBoth = false) => {
      const a = read(first), b = read(last);
      if (!Number.isFinite(a) && !Number.isFinite(b)) return null;
      if (requireBoth && (!Number.isFinite(a) || !Number.isFinite(b))) {
        throw new Error('Kapsanacak hedefin başını ve sonunu birlikte girin.');
      }
      return [Number.isFinite(a) ? a : 1, Number.isFinite(b) ? b : t.size];
    };
    const constraints = readCons(el);
    const payload = {tool:'design', region:t.region, reference_sequence:t.sequence,
      profile:$('#prof', el).value, conversion:$('#conversion', el).value,
      constraints, salt_model:constraints.salt_model, top:+$('#nRet', el).value,
      conditions:{na:read('na'), mg:read('mg'), dntp:read('dntp'), pnm:read('pnm'), ta:read('ta')},
      target_range:range('tFrom', 'tTo'), left_region:range('fFrom', 'fTo'),
      right_region:range('rFrom', 'rTo'), must_cover:range('cFrom', 'cTo', true),
      variant_filter:$('#vOn', el).checked, variant_af:+$('#vAf', el).value,
      variant_scope:$('#vScope', el).value, search_mode:$('#sMode', el)?.value || 'anchored'};
    if (!form.reportValidity()) return;
    button.disabled = true;
    STATE.design = null;
    $('#dres', el).replaceChildren();
    runLocalJob('design', payload, (result, id) => {
      acceptLiveDesign(result, id);
      if (STATE.view === 'tasarla') showDesign($('#view'), STATE.design);
      else if (STATE.view === 'cift' || STATE.view === 'yerel') rerender();
    }).finally(() => {
      if (STATE.view !== 'tasarla') return;
      const currentButton = $('#goBtn');
      if (currentButton) currentButton.disabled = !liveReady();
      if (LIVE.errors.has('design')) $('#dres').innerHTML = designIntro();
    });
  } catch (error) {
    LIVE.errors.set('design', error.message);
    updateJobNotice('design');
  }
}

function makeLivePairs(result, source) {
  const t = result.template;
  return result.pairs.map((pair, index) => {
    const p = pair.design;
    return {key:'L:' + source + ':' + pair.id, kind:'live', source, id:pair.id,
      label:'Tasarım #' + (index + 1), short:'#' + (index + 1), rank:index + 1,
      F:pair.F, R:pair.R, fStart:p.left.start, fEnd:p.left.end - 1,
      rStart:p.right.start, rEnd:p.right.end - 1, prod:p.product_length,
      tmF:p.left.tm_c, tmR:p.right.tm_c, dtm:p.delta_tm, thetaF:p.left.theta_ta,
      thetaR:p.right.theta_ta, score:p.score, pareto:p.pareto_rank, d:p,
      tpl:t.sequence, tplOff:t.offset, tplName:t.region + ' GRCh38.p14',
      chrom:t.chrom, conversion:t.conversion, hasVar:false, profile:result.profile,
      cond:result.conditions, salt:result.salt_model};
  });
}

function acceptLiveDesign(result, source) {
  const t = result.template;
  const pms = makeLivePairs(result, source);
  const live = {result, source};
  const recentSources = [...new Set(STATE.local.map(pm => pm.source).filter(Boolean))].slice(-2);
  STATE.local = [...STATE.local.filter(pm => recentSources.includes(pm.source)), ...pms];
  STATE.design = {live, pms, ms:result.elapsed_ms, cons:result.constraints,
    cond:result.conditions, t:{seq:t.sequence, off:t.offset, name:t.region, hasVar:false},
    prof:result.profile, vOn:result.variant_filter};
  LIVE.source = pms.length ? source : (STATE.local.at(-1)?.source || '');
  LIVE.selected = pms[0]?.key || null;
  LIVE.selectedKeys = null;
  if (pms.length) STATE.pair = pms[0].key;
  if (result.specificity) {
    applyScan(source, result.specificity);
    LIVE.last = {tool:'scan', result:result.specificity};
  }
  for (const pm of pms) {
    const qc = result.qc?.[pm.id];
    if (qc) resultForPair(pm, 'qc', qc);
  }
}

const localToolLabels = {
  scan: 'Özgüllük / Multipleks', qc: 'Dizi, dimer ve termodinamik QC',
  variants: 'gnomAD varyantları', isoform: 'GENCODE izoformları',
  probe: 'Hidroliz probu', melt: 'Amplikon erime profili',
  thermo: 'Tek primer termodinamiği', blocker: 'Prob ve blocker QC',
};

function localToolLink(tool, label, pm) {
  return `<button class="btn sm" data-live-tool="${tool}" ${pm ? 'data-live-pair="' + esc(pm.key) + '"' : ''}>${esc(label)}</button>`;
}

document.addEventListener('click', e => {
  const link = e.target.closest('[data-live-tool]');
  if (!link) return;
  LIVE.tool = link.dataset.liveTool;
  LIVE.selectedKeys = null;
  LIVE.selected = link.dataset.livePair || null;
  if (LIVE.selected) {
    const pm = getPair(LIVE.selected);
    LIVE.source = pm?.source || (pm?.kind === 'project' ? 'project' : '');
  } else if (LIVE.source === 'project') LIVE.source = STATE.local.at(-1)?.source || '';
  go('yerel');
});

function liveValue(value) {
  if (value == null) return '—';
  return typeof value === 'object' ? evidenceJson(value) : fmtv(value);
}

function localRows(rows, columns) {
  if (!rows?.length) return '<p class="small faint">Kayıt yok.</p>';
  const cols = columns || Object.keys(rows[0]).map(key => [key, tidyKey(key)]);
  return `<div class="scrollx"><table class="t"><thead><tr>${cols.map(([, label]) =>
    `<th scope="col">${esc(label)}</th>`).join('')}</tr></thead><tbody>${rows.map(row =>
    '<tr>' + cols.map(([key]) => `<td>${esc(liveValue(row[key]))}</td>`).join('') + '</tr>').join('')}</tbody></table></div>`;
}

function liveResult(title, body) {
  return `<section class="card"><div class="hd"><h2>${esc(title)}</h2></div><div class="bd">${body}</div></section>`;
}

function localRecords(records) {
  if (!records?.length) return '<p class="small faint">Kanıt kaydı yok.</p>';
  return `<div class="live-records">${records.map(r => `<details><summary>${esc(r.id)} · ${esc(r.role || CATI[r.id]?.name || '')} ${sev(r.severity || 'NA')}</summary>
    <dl><dt>Ölçüm</dt><dd>${esc(liveValue(r.value))}</dd><dt>Eşik</dt><dd>${esc(liveValue(r.threshold))}</dd><dt>Model</dt><dd>${esc(liveValue(r.model))}</dd></dl>
    ${r.note ? `<p class="small muted">${esc(r.note)}</p>` : ''}</details>`).join('')}</div>`;
}

function saveLocalResult(data) {
  const blob = new Blob([JSON.stringify(data, null, 2)], {type:'application/json'});
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = 'oligoscout-yerel-sonuc.json';
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function localSourcePairs(source) {
  return source === 'project' ? PROJ : STATE.local.filter(pm => pm.kind === 'live' && pm.source === source);
}

function resultForPair(pm, action, result) {
  if (action === 'qc') {
    pm.liveQC = result;
    pm.liveVariants = result.variants;
    pm.hasVar = result.variants.evaluated;
  } else if (action === 'variants') {
    pm.liveVariants = result;
    pm.hasVar = result.evaluated;
  } else if (action === 'melt' || action === 'probe' || action === 'isoform') {
    pm.liveExtras ??= {};
    pm.liveExtras[action] = result;
  }
}

function applyScan(source, result) {
  for (const pm of localSourcePairs(source)) {
    if (result.bisulfite) {
      if (result.summary.some(row => row.pair_id === pm.id)) pm.liveBis = result;
    } else if (result.pairs[pm.id]) {
      const scan = result.pairs[pm.id];
      pm.liveScan = {result: scan, contract: result.contracts.primary,
        records: result.records, multiplex: result.multiplex};
      if (scan.variants) resultForPair(pm, 'variants', scan.variants);
    }
  }
}

function startPairJob(pm, action) {
  const key = action + ':' + pm.key;
  runLocalJob(key, {tool:action, source:pm.source || 'project', pairs:[pm.id]}, result => {
    if (action === 'scan') applyScan(pm.source || 'project', result);
    else resultForPair(pm, action, result);
    if ((STATE.view === 'cift' && STATE.pair === pm.key) ||
        (STATE.view === 'dokumantasyon' && STATE.docSection === 'cift' && STATE.docPair === pm.key)) rerender();
  });
}

function pairAction(pm, tool, label) {
  const key = tool + ':' + pm.key;
  return liveReady()
    ? `<button class="btn sm" data-pair-job="${tool}" ${LIVE.pending.has(key) ? 'disabled' : ''}>${esc(label)}</button>`
    : '';
}

function bindPairActions(host, pm) {
  $$('[data-pair-job]', host).forEach(button => button.onclick = () => {
    button.disabled = true;
    const key = button.dataset.pairJob + ':' + pm.key;
    startPairJob(pm, button.dataset.pairJob);
    const wait = () => {
      if (!button.isConnected) return;
      if (!LIVE.pending.has(key)) button.disabled = false;
      else setTimeout(wait, 600);
    };
    wait();
  });
}

function liveEvidence(pm) {
  const evidence = {};
  const add = (record, source, role) => {
    if (!record?.id) return;
    const group = evidence[record.id] ??= {status:'eval', recs:[]};
    group.recs.push(rec(role || record.role, source, record));
  };
  pm.liveQC?.records.forEach(r => add(r, 'Kalite kontrolü'));
  pm.liveVariants?.records?.forEach(r => add(r, pm.liveVariants.source || 'gnomAD v4'));
  const scan = pm.liveScan;
  if (scan) {
    for (const role of ['F', 'R']) {
      scan.result.primers[role].forEach(out => out.checks?.filter(r => r.id.startsWith('SPC-') || !pm.liveQC)
        .forEach(r => add(r, 'Genom özgüllüğü', role)));
    }
    scan.result.records.forEach(r => add(r, 'Genom özgüllüğü', 'pair'));
    scan.records.forEach(r => add(r, 'Genom özgüllüğü', 'panel'));
    scan.multiplex?.summary.forEach(r => add({...r, id:r.check}, 'multiplex_qc.py', 'panel'));
  }
  if (pm.liveBis) {
    const check = pm.liveBis.records.find(r => r.pair_id === pm.id);
    if (check) add(check.check || check, 'Bisülfit özgüllüğü');
  }
  if (pm.liveExtras) {
    for (const [key, result] of Object.entries(pm.liveExtras)) {
      add(result.check, key === 'isoform' ? 'GENCODE v47 · isoform_qc.py' : 'amplicon.py');
      result.result?.checks?.forEach(r => add(r, 'probe.py'));
    }
  }
  for (const group of Object.values(evidence)) group.sev = worst(group.recs.map(r => r.sev).filter(Boolean));
  for (const c of CAT) evidence[c.id] ??= {status:'ne', recs:[], reason:'Bu çift için henüz değerlendirilmedi.'};
  return evidence;
}

function qualityCheckCards(E, prefixes) {
  const order = new Map(CAT.map((r, i) => [r.id, i]));
  const ids = Object.keys(E).filter(id => prefixes.some(p => id.startsWith(p)) && E[id].recs.some(r => !r.ne))
    .sort((a, b) => (order.get(a) ?? 999) - (order.get(b) ?? 999) || a.localeCompare(b));
  if (!ids.length) return '';
  return liveResult('Kalite kontrolü · kanıt kayıtları', `<p class="small muted" style="margin:0 0 12px">Ölçüm, eşik ve yapı çizimi için bir karta tıklayın.</p><div class="ccgrid">${ids.map(id => ckCard(id, E[id])).join('')}</div>`);
}

function renderLocalOutput(tool, result) {
  if (tool === 'scan') {
    if (result.bisulfite) return liveResult('BIS-04 · CT ve GA taraması',
      localRows(result.summary) + localRecords(result.records.map(r => r.check || r)));
    const rows = Object.values(result.pairs).map(p => p.summary).filter(Boolean);
    return liveResult('GRCh38 · Özgüllük taraması', localRows(rows, [
      ['pair_id','Çift'], ['products_intended','Hedef ürün'], ['offtarget_total','Hedef dışı'],
      ['offtarget_exact','Tam hedef dışı'], ['sites_f','F yeri'], ['sites_r','R yeri'],
      ['spc01','SPC-01'], ['spc03','SPC-03'],
    ]) + localRecords(result.records)) + (result.ranking ? liveResult('Özgüllük ve varyant sonrası sıralama',
      localRows(result.ranking, [['final_rank','Sıra'],['pair_id','Çift'],['verdict','Karar'],
        ['specificity_verdict','Özgüllük'],['variant_verdict','Varyant'],['thermo_rank','Pareto katmanı']])) : '')
      + (result.multiplex ? liveResult(
      'Multipleks · çapraz ürünler ve dimerler', localRows(result.multiplex.summary) +
      `<details><summary>Dimer matrisi</summary>${localRows(result.multiplex.dimer)}</details>
      <details><summary>Çapraz ürünler · ilk 200</summary>${localRows(result.multiplex.cross)}</details>`) : '');
  }
  if (tool === 'variants') return result.evaluated
    ? liveResult(result.source || 'gnomAD · varyant kanıtları', localRecords(result.records))
    : liveNote(result.reason, true);
  if (tool === 'melt') {
    const plot = result.curve.map(row => [row[0], row[2]]);
    return liveResult('Erime eğrisi · −dθ/dT', `<div class="scrollx" id="liveMelt"></div>${kvHTML(result.summary)}${localRecords([result.check])}
      <p class="tiny faint">${plot.length} sıcaklık noktası · amplikon modeli</p>`);
  }
  if (tool === 'probe') return result.result
    ? liveResult('Hidroliz probu değerlendirmesi', localRecords(result.result.checks))
    : liveResult('Uygun prob adayı bulunamadı', kvHTML(result.rejections));
  if (tool === 'isoform') return liveResult('RT-03 · GENCODE izoform kapsamı', localRecords([result.check]));
  return liveResult('Analiz kayıtları', localRecords(result.records || result.result?.checks));
}

function localConditionsHTML() {
  const fields = [['na','Na⁺ / K⁺ (mM)','50','.1'],['mg','Mg²⁺ (mM)','1.5','.01'],
    ['dntp','dNTP (mM)','0.6','.01'],['pnm','Primer (nM)','200','.1'],['ta','Ta (°C)','60','.1']];
  return `<div class="fgrid c3">${fields.map(([id,label,value,step]) =>
    `<label class="f"><span>${label}</span><input name="${id}" type="number" min="0" step="${step}" value="${value}" required></label>`).join('')}</div>`;
}

RENDER.yerel = el => {
  const tool = LIVE.tool;
  const pairTool = tool !== 'thermo' && tool !== 'blocker';
  const source = LIVE.source;
  const localSources = new Map(STATE.local.filter(pm => pm.kind === 'live')
    .map(pm => [pm.source, pm.sourceLabel || pm.tplName]));
  const pairs = localSourcePairs(source);
  const choice = pairs.find(pm => pm.key === LIVE.selected) || pairs[0];
  const key = 'workbench';
  const selectedKeys = LIVE.selectedKeys || [choice?.key];
  const control = (name, label, required = true) =>
    `<label class="f"><span>${label}</span><input name="${name}" spellcheck="false" maxlength="500" ${required ? 'required' : ''}></label>`;
  el.innerHTML = `<div class="stack"><section class="card"><div class="hd"><h2>Yerel analiz</h2>
    </div>
    <form class="bd stack" id="localForm" autocomplete="off">
      <label class="f"><span>Analiz türü</span><select id="localTool">${Object.entries(localToolLabels).map(([id,label]) =>
        `<option value="${id}" ${tool === id ? 'selected' : ''}>${label}</option>`).join('')}</select></label>
      ${pairTool ? `<label class="f"><span>Primer kaynağı</span><select id="localSource">
        <option value="">Primer kaynağı seçin</option>
        ${source === 'project' ? '<option value="project" selected>Dokümantasyon · seçilen örnek</option>' : ''}
        ${[...localSources].reverse().map(([id, label]) => `<option value="${esc(id)}" ${source === id ? 'selected' : ''}>
        Yerel tasarım · ${esc(label)} · ${esc(id.slice(0, 5))}</option>`).join('')}
        </select></label>
        ${!pairs.length ? `<div class="note">${ic('info')}<div>Bu analiz için bir primer çifti gerekir. <a href="#tasarla">Primer Tasarla</a> sayfasında kendi hedefinizi çalıştırın.</div></div>` :
        tool === 'scan' ? `<fieldset class="live-fieldset"><legend>Taranacak çiftler</legend><button type="button" class="btn sm" id="selectAllPairs">Tümünü seç</button>
          <div class="live-pairs">${pairs.map(pm => `<label class="switch"><input type="checkbox" name="pair" value="${esc(pm.key)}" ${selectedKeys.includes(pm.key) ? 'checked' : ''}>
          ${esc(pm.id)} · ${pm.prod} bp</label>`).join('')}</div></fieldset>` :
        `<label class="f"><span>Primer çifti</span><select name="pair">${pairs.map(pm =>
          `<option value="${esc(pm.key)}" ${choice?.key === pm.key ? 'selected' : ''}>${esc(pm.id)} · ${pm.prod} bp</option>`).join('')}</select></label>`}` : ''}
      ${tool === 'scan' ? `<p class="small muted">${choice?.conversion === 'CT' || choice?.conversion === 'GA' ? 'Her iki bisülfit indeksi taranır.' : '24 kromozom taranır. Birden fazla çift seçilirse multipleks sonuçları da üretilir.'} Tarama dizin yükleme süresine bağlı olarak birkaç dakika sürebilir.</p>
      <div class="fgrid c3">${[['anchor','3′ çapa (nt)',12,8,24],['mismatches','En fazla uyumsuzluk',3,0,3],['min_product','En kısa ürün (bp)',50,20,5000],['max_product','En uzun ürün (bp)',3000,20,5000]].map(([id,label,val,min,max]) =>
        `<label class="f"><span>${label}</span><input name="${id}" type="number" min="${min}" max="${max}" value="${val}" required></label>`).join('')}
        <label class="f"><span>Arama modu</span><select name="search_mode"><option value="anchored" selected>Hızlı · tam 3′ çapa</option><option value="sensitive">Hassas · Primer-BLAST kuralı</option></select></label></div>` : ''}
      ${tool === 'isoform' ? control('gene','GENCODE gen simgesi') : ''}
      ${tool === 'probe' ? control('probe','Var olan prob (boş bırakılırsa aday aranır)',false) : ''}
      ${tool === 'thermo' ? control('sequence','Primer dizisi (5′ → 3′; LNA: +N)') + localConditionsHTML() : ''}
      ${tool === 'blocker' ? `<div class="fgrid c3">${control('mutant','Mutant prob')}${control('wt','WT prob')}${control('blocker','Blocker (isteğe bağlı)',false)}${control('forward','İleri primer (isteğe bağlı)',false)}${control('reverse','Geri primer (isteğe bağlı)',false)}</div>${localConditionsHTML()}` : ''}
      ${!liveReady() ? liveNote(localAvailability()) : ''}
      ${jobNotice(key)}
      <div class="row"><button id="runLocalTool" type="submit" class="btn pri" ${!liveReady() || LIVE.pending.has(key) || (pairTool && !pairs.length) ? 'disabled' : ''}>Analizi çalıştır</button></div>
    </form></section>
    ${LIVE.last?.tool === tool ? `<div class="row"><button class="btn sm" id="localDownload">Sonucu JSON olarak indir</button>
      ${pairTool && choice ? `<a class="btn sm" href="${choice.kind === 'project' ? '#dokumantasyon/cift/' + choice.key : '#cift'}">Çift incelemeyi aç</a>` : ''}</div>${renderLocalOutput(tool, LIVE.last.result)}` : ''}
  </div>`;
  const form = $('#localForm', el);
  $('#localTool', el).onchange = e => { LIVE.tool = e.target.value; rerender(); };
  const sourceSelect = $('#localSource', el);
  if (sourceSelect) sourceSelect.onchange = e => { LIVE.source = e.target.value; LIVE.selected = null; LIVE.selectedKeys = null; rerender(); };
  const all = $('#selectAllPairs', el);
  if (all) all.onclick = () => { $$('input[name="pair"]', form).forEach(b => b.checked = true); };
  form.onsubmit = event => {
    event.preventDefault();
    if (!form.reportValidity()) return;
    const fd = new FormData(form);
    const keys = fd.getAll('pair');
    const selected = pairs.filter(pm => keys.includes(pm.key));
    if (pairTool && !selected.length) {
      LIVE.errors.set(key, 'En az bir primer çifti seçin.');
      updateJobNotice(key);
      return;
    }
    LIVE.selectedKeys = keys;
    LIVE.selected = selected[0]?.key;
    const body = {tool, source, pairs:selected.map(pm => pm.id)};
    if (tool === 'scan') { for (const name of ['anchor','mismatches','min_product','max_product']) body[name] = +fd.get(name); body.search_mode = fd.get('search_mode') || 'anchored'; }
    if (tool === 'isoform') body.gene = fd.get('gene');
    if (tool === 'probe') body.probe = fd.get('probe');
    if (tool === 'thermo' || tool === 'blocker') {
      for (const name of ['sequence','mutant','wt','blocker','forward','reverse']) if (fd.has(name)) body[name] = fd.get(name);
      body.conditions = liveConditions(form);
    }
    const button = $('#runLocalTool', form);
    button.disabled = true;
    runLocalJob(key, body, result => {
      LIVE.last = {tool, result};
      if (tool === 'scan') applyScan(body.source, result);
      else if (pairTool) resultForPair(selected[0], tool, result);
      if (selected.length) STATE.pair = selected[0].key;
      if (STATE.view === 'yerel' && LIVE.tool === tool) rerender();
    }).finally(() => { if (button.isConnected) button.disabled = false; });
  };
  const down = $('#localDownload', el);
  if (down) down.onclick = () => saveLocalResult(LIVE.last.result);
  const plot = $('#liveMelt', el);
  if (plot) {
    plot.innerHTML = lineChart(Math.max(320, W_OF(plot)), 220, [{
      pts:LIVE.last.result.curve.map(r => [r[0], r[2]]), c:css('--accent'), area:true,
    }], {xlabel:'Sıcaklık (°C)', ylabel:'−dθ/dT'});
  }
};
