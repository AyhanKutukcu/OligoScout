/* ================= ÇİFT İNCELEME ================= */
const PTABS = [
  {id:'harita', t:'Harita', i:'map'}, {id:'ozgulluk', t:'Özgüllük', i:'target'}, {id:'varyant', t:'Varyant', i:'snp'},
  {id:'yapi', t:'Yapı & Dimer', i:'struct'}, {id:'termo', t:'Termodinamik', i:'thermo'}, {id:'kanit', t:'Kanıt', i:'list'}];
RENDER.cift = el => {
  const example = STATE.view === 'dokumantasyon';
  const pool = example ? PROJ : STATE.local;
  const pm = pool.find(p => p.key === (example ? STATE.docPair : STATE.pair)) ||
    (example ? pool[0] : STATE.design?.pms?.[0] || pool.at(-1));
  if (!pm) {
    el.innerHTML = liveResult('Henüz bir primer çifti yok', `<p class="muted">Tasarım ve otomatik özgüllük taraması tamamlandığında çiftinizi burada inceleyebilirsiniz.</p>
      <div class="row"><a class="btn pri" href="#tasarla">Primer tasarla</a><a class="btn" href="#dokumantasyon/panel">Örnek çıktıyı gör</a></div>`);
    return;
  }
  if (example) STATE.docPair = pm.key; else STATE.pair = pm.key;
  const E = evidenceFor(pm); pm._E = E;
  const famSev = fams => {
    const assessed = Object.entries(E)
      .filter(([id]) => fams.some(f => id.startsWith(f)))
      .map(([, group]) => group.sev).filter(Boolean);
    return assessed.length ? worst(assessed) : null;
  };
  const tabSev = {ozgulluk: pm.kind === 'project' || pm.liveScan ? famSev(['SPC-']) : null, varyant: pm.hasVar ? famSev(['VAR-']) : null, yapi: famSev(['DIM-0','SS-']), termo: famSev(['TH-']), kanit: famSev([''])};
  const opts = pool.map(p => `<option value="${esc(p.key)}" ${p.key === pm.key ? 'selected' : ''}>${esc(example ? p.id : p.label + ' · ' + p.tplName)} · ${p.prod} bp</option>`).join('');
  const amp = pm.tpl.slice(pm.fStart, pm.rEnd + 1);
  const ampFa = `>${pm.id} ${pm.chrom}:${gposAny(pm, pm.fStart)}-${gposAny(pm, pm.rEnd)} ${pm.prod} bp\n${(amp.match(/.{1,60}/g) || []).join('\n')}`;
  el.innerHTML = `<div class="stack">
    <div class="card"><div class="bd pairhead">
      <div style="min-width:0">
        <div class="row" style="gap:8px;margin-bottom:8px"><select id="pairSel" style="width:auto;max-width:100%;font-weight:600">${opts}</select>
          <span class="tag ${pm.kind === 'project' ? 'eng' : 'br'}">${pm.kind === 'project' ? 'motor çıktısı · ' + esc(D.meta.qc_dir) : (pm.kind === 'live' ? 'Primer tasarımı · ' : 'tarayıcı ön tasarımı · ') + esc(pm.profile)}</span></div>
        <div style="display:grid;gap:6px">
          <div class="primerline"><span class="r F">F</span><span class="seq">${seqC(pm.F)}</span><span class="faint small">${pm.F.length} nt · ${pm.chrom}:${gposAny(pm, pm.fStart)}(+)</span>${copyBtn(pm.F)}</div>
          <div class="primerline"><span class="r R">R</span><span class="seq">${seqC(pm.R)}</span><span class="faint small">${pm.R.length} nt · ${pm.chrom}:${gposAny(pm, pm.rEnd)}(−)</span>${copyBtn(pm.R)}</div>
          <div class="small muted">Amplikon <b class="mono">${pm.chrom}:${gposAny(pm, pm.fStart)}–${gposAny(pm, pm.rEnd)}</b> · ${pm.prod} bp · GC %${nf(OS.gcf(amp)*100,1)} ${copyBtn(ampFa, 'Amplikonu FASTA olarak kopyala')}</div>
        </div>
      </div>
      <div style="display:grid;gap:10px;justify-items:end">
        <div class="row" style="gap:6px;justify-content:flex-end">${pm.kind === 'project' ? `<span class="chip">karar ${sev(pm.verdict)}</span><span class="chip">özgüllük ${sev(pm.spec)}</span><span class="chip">varyant ${sev(pm.varv)}</span>` : `<span class="chip">Pareto <b>${pm.pareto}</b></span><span class="chip">özgüllük ${pm.liveScan ? sev(worst([pm.liveScan.result.summary?.spc01, pm.liveScan.result.summary?.spc03])) : pm.liveBis ? sev('INFO','BIS-04 hazır') : sev('NA','taranmadı')}</span>`}</div>
        <div class="mini-stats"><div><div class="k">Ürün</div><div class="v">${pm.prod}<small class="faint" style="font-size:11px"> bp</small></div></div><div><div class="k">Tm F / R</div><div class="v">${nf(pm.tmF,1)} / ${nf(pm.tmR,1)}</div></div><div><div class="k">ΔTm</div><div class="v">${nf(pm.dtm,2)}</div></div><div><div class="k">θ(Ta) F / R</div><div class="v">${nf(pm.thetaF,2)} / ${nf(pm.thetaR,2)}</div></div></div>
      </div></div>
      <div class="tabs pair-tabs" role="tablist" style="padding:0 10px">${PTABS.map(t => `<button role="tab" data-ptab="${t.id}" aria-selected="${STATE.ptab === t.id}">${ic(t.i)}${t.t}${tabSev[t.id] ? ' ' + sev(tabSev[t.id]) : ''}</button>`).join('')}</div>
    </div>
    <div id="ptab"></div></div>`;
  $('#pairSel', el).onchange = e => {
    if (example) go('dokumantasyon/cift/' + e.target.value + '/' + STATE.ptab);
    else { STATE.pair = e.target.value; rerender(); }
  };
  $$('[data-ptab]', el).forEach(b => b.onclick = () => {
    if (example) go('dokumantasyon/cift/' + pm.key + '/' + b.dataset.ptab);
    else { STATE.ptab = b.dataset.ptab; rerender(); }
  });
  const host = $('#ptab', el);
  ({harita: tabMap, ozgulluk: tabSpec, varyant: tabVar, yapi: tabStruct, termo: tabThermo, kanit: tabEvidence})[STATE.ptab](host, pm, E);
  if (STATE.ptab === 'kanit' && liveReady() && !example && !pm.liveQC) {
    host.insertAdjacentHTML('afterbegin', `<div class="row live-actions">${pairAction(pm, 'qc', 'Kalite kontrolü yap')}</div>${jobNotice('qc:' + pm.key)}`);
    bindPairActions(host, pm);
  }
};

/* ---- harita ---- */
function tabMap(host, pm) {
  const others = pm.kind === 'project' ? PROJ.filter(p => p.key !== pm.key) : STATE.local.filter(p => p.key !== pm.key && p.source === pm.source);
  const pad = STATE.mapPad || 24;
  host.innerHTML = `<div class="stack">
    <div class="card"><div class="hd"><h2>Şablondaki konum</h2><span class="faint small">${esc(pm.tplName || D.meta.header)}</span></div><div class="bd scrollx" id="mTrack"></div></div>
    <div class="card"><div class="hd"><h2>Dizi haritası</h2><span class="faint small">iki iplik · primer bağlanma yerleri · gnomAD varyantları</span><span class="sp"></span>
      <div class="seg" id="padSeg"><button data-pad="12" aria-pressed="${pad===12}">±12 nt</button><button data-pad="24" aria-pressed="${pad===24}">±24 nt</button><button data-pad="60" aria-pressed="${pad===60}">±60 nt</button></div></div>
      <div class="bd" style="padding-bottom:6px"><div class="legend" style="margin-bottom:10px"><span><i class="dot" style="background:var(--fwd)"></i>Forward (+ iplik, 5′→3′)</span><span><i class="dot" style="background:var(--rev)"></i>Reverse (− iplik)</span><span><i class="dot" style="background:var(--amp);border:1px solid var(--accent)"></i>Amplikon</span>${pm.hasVar ? `<span>▼ <span style="color:var(--fail)">AF ≥ %1</span> · <span style="color:var(--warn)">AF ≥ %0,1</span> · <span class="faint">nadir</span></span>` : ''}<span class="faint">Varyant üzerine gelin: rsID, alel, AF, popülasyon</span></div>
      <div class="scrollx" id="sMap"></div></div></div>
    <div class="g2">
      <div class="card"><div class="hd"><h2>Amplikon dizisi</h2><span class="sp"></span>${copyBtn(pm.tpl.slice(pm.fStart, pm.rEnd + 1), 'Amplikonu kopyala')}</div><div class="bd"><div class="seq small" style="line-height:1.9">${ampHTML(pm)}</div></div></div>
      <div class="card"><div class="hd"><h2>Amplikon GC profili</h2><span class="faint small">10 nt kayan pencere · AMP-01</span></div><div class="bd" id="gcChart"></div></div>
    </div></div>`;
  const w = W_OF($('#mTrack', host)) - 36;
  $('#mTrack', host).innerHTML = overviewTrack(Math.max(640, w), {tpl: pm.tpl, off: pm.tplOff, pairs: [pm, ...others.slice(0, 8)].map(p => ({...p, short: p.key === pm.key ? '▶ ' + (p.short || p.id.replace('DESIGN_','#')) : (p.short || p.id.replace('DESIGN_','#'))})), sel: pm.key, showVars: pm.hasVar && pm.tpl === TPL, region: [pm.fStart, pm.rEnd], compact: true, win: [Math.max(0, Math.min(pm.fStart, ...others.slice(0,8).map(o => o.fStart)) - 30), Math.min(pm.tpl.length, Math.max(pm.rEnd, ...others.slice(0,8).map(o => o.rEnd)) + 30)]});
  $('#sMap', host).innerHTML = seqMap(Math.max(560, W_OF($('#sMap', host))), pm, {pad});
  $$('#padSeg button', host).forEach(b => b.onclick = () => { STATE.mapPad = +b.dataset.pad; rerender(); });
  const amp = pm.tpl.slice(pm.fStart, pm.rEnd + 1); const pts = [];
  for (let i = 0; i + 10 <= amp.length; i++) pts.push([i + 5, OS.gcf(amp.slice(i, i + 10)) * 100]);
  $('#gcChart', host).innerHTML = lineChart(Math.max(300, W_OF($('#gcChart', host)) - 36), 170, [{pts, c: css('--accent'), area: true}], {y0: 0, y1: 100, x0: 0, x1: amp.length, yticks: [0, 25, 50, 75, 100], ylabel: 'GC %', xlabel: 'amplikon konumu (nt)',
    bands: [{a: 0, b: pm.F.length, c: css('--fwd'), op: .12}, {a: amp.length - pm.R.length, b: amp.length, c: css('--rev'), op: .12}], vlines: []});
}
function ampHTML(pm) {
  const a = pm.tpl.slice(pm.fStart, pm.rEnd + 1), fL = pm.F.length, rL = pm.R.length; const vv = pm.hasVar ? new Set(varsIn(pm.fStart, pm.rEnd, pm).filter(v => v.af >= 0.001).map(v => v.idx - pm.fStart)) : new Set();
  return a.split('').map((c, i) => { const inF = i < fL, inR = i >= a.length - rL; const st = inF ? 'color:var(--fwd);font-weight:700' : inR ? 'color:var(--rev);font-weight:700' : 'color:var(--ink2)';
    return `<span style="${st}${vv.has(i) ? ';background:var(--warn-soft);border-bottom:2px solid var(--warn)' : ''}">${c}</span>${(i + 1) % 10 === 0 ? ' ' : ''}`; }).join('');
}

/* ---- özgüllük ---- */
const CHRLEN = {chr1:248956422,chr2:242193529,chr3:198295559,chr4:190214555,chr5:181538259,chr6:170805979,chr7:159345973,chr8:145138636,chr9:138394717,chr10:133797422,chr11:135086622,chr12:133275309,chr13:114364328,chr14:107043718,chr15:101991189,chr16:90338345,chr17:83257441,chr18:80373285,chr19:58617616,chr20:64444167,chr21:46709983,chr22:50818468,chrX:156040895,chrY:57227415};
function tabSpec(host, pm, E) {
  const scan = pm.liveScan;
  const actions = `<div class="row live-actions">${pairAction(pm, 'scan', scan || pm.liveBis ? 'Yeniden tara' : 'Özgüllük taraması')}
    ${localToolLink('scan', 'Panel ve arama ayarları', pm)}</div>${jobNotice('scan:' + pm.key)}`;
  if (pm.liveBis) {
    const result = pm.liveBis;
    const rows = result.summary.filter(row => row.pair_id === pm.id);
    host.innerHTML = actions + liveResult('BIS-04 · CT ve GA özgüllüğü',
      localRows(rows) + localRecords(result.records.filter(r => r.pair_id === pm.id).map(r => r.check || r))) +
      liveResult('CT / GA arama sözleşmeleri',
        Object.entries(result.contracts).map(([conv, contract]) => `<h3>${conv}</h3>${kvHTML(contract)}`).join(''));
    bindPairActions(host, pm);
    return;
  }
  if (pm.kind !== 'project' && !scan) {
    host.innerHTML = `<div class="stack">${actions}${liveNote(liveReady()
      ? 'Bu çift henüz taranmadı.'
      : localAvailability())}</div>`;
    bindPairActions(host, pm);
    return;
  }
  const contract = scan?.contract || D.meta.contract;
  const p = scan ? {...scan.result, panelqc:scan.result.summary} : pm.raw;
  const q = p.panelqc;
  const ex = [+q.exact_f, +q.exact_r], ot = +q.offtarget_total;
  const intended = +q.products_intended;
  const specific = intended === 1 && ot === 0 && ex[0] === 1 && ex[1] === 1;
  const summarySeverity = scan ? worst([q.spc01, q.spc03, E['SPC-01']?.sev]) : p.specificity_verdict;
  const heading = specific ? 'Arama kapsamında tek hedef ürün' : intended !== 1
    ? 'Beklenen hedef ürün doğrulanamadı' : ot ? ot + ' hedef dışı ürün bulundu'
    : 'Primer tam eşleşme sayıları incelenmeli';
  const mmHist = role => { const h = {0:0,1:0,2:0,3:0}; p.bindings[role].hist.forEach(x => h[x.mm] = (h[x.mm] || 0) + x.n); return h; };
  const hF = mmHist('F'), hR = mmHist('R');
  host.innerHTML = `<div class="stack">${actions}
    <div class="verdict ${specific ? 'ok' : summarySeverity === 'FAIL' ? 'fail' : 'warn'}"><div class="ic">${ic(specific ? 'ok' : 'warn')}</div><div style="min-width:0"><h3>${heading}</h3>
      <div class="small">${scan ? 'Genom taraması sonucu' : 'Kayıtlı panel sonucu'} · ${q.products_total} toplam ürün · ${intended} hedef · ${ot} hedef dışı. ${contract.shards} kromozom, primer başına en fazla ${contract.max_mismatches} uyumsuzluk, ${contract.min_amplicon_length}–${contract.max_amplicon_length} bp ürün.</div>
      </div><span class="sp" style="flex:1"></span>${sev(summarySeverity)}</div>
    <div class="ccgrid">${['SPC-01','SPC-02','SPC-03','SPC-05','SPC-06','SPC-07'].map(id => ckCard(id, E[id])).join('')}</div>
    <div class="g2">
      <div class="card"><div class="hd"><h2>Uyumsuzluk sınıfı histogramı</h2><span class="faint small">SPC-02 · tüm bağlanma yerleri 3′ ${contract.anchor_length} nt çapalı</span></div><div class="bd" id="mmHist"></div></div>
      <div class="card"><div class="hd"><h2>Genomdaki bağlanma yerleri</h2><span class="faint small">en yakın 40 yer / primer</span></div><div class="bd scrollx" id="ideo"></div></div>
    </div>
    <div class="card"><div class="hd"><h2>Bağlanma yeri hizalamaları</h2><span class="faint small">Primer-BLAST gösterimi · nokta = eşleşme, harf = uyumsuz baz</span><span class="sp"></span>
      <div class="seg" id="roleSeg"><button data-role="F" aria-pressed="${(STATE.specRole||'F')==='F'}">Forward · ${p.bindings.F.n}</button><button data-role="R" aria-pressed="${STATE.specRole==='R'}">Reverse · ${p.bindings.R.n}</button></div>
      <div class="seg" id="mmSeg"><button data-mm="9" aria-pressed="${(STATE.specMM??9)===9}">tümü</button><button data-mm="1" aria-pressed="${STATE.specMM===1}">≤1</button><button data-mm="2" aria-pressed="${STATE.specMM===2}">≤2</button></div></div>
      <div class="bd flush" id="alnList"></div></div>
    <div class="card"><div class="hd"><h2>Hedef dışı ürünler</h2><span class="faint small">SPC-03 · F–R, R–F ve tek primer (F–F, R–R) ürünleri · ${contract.min_amplicon_length}–${contract.max_amplicon_length} bp</span></div><div class="bd ${p.offtarget.length ? 'flush' : ''}">${offtHTML(pm, p, contract, q)}</div></div>
    <div class="card"><div class="hd"><h2>Arama sözleşmesi</h2><span class="faint small">SPC-09 · bu sonuç şu indeks ve sözleşmeyle üretildi</span></div><div class="bd">${kvHTML(scan ? {...contract, ...scan.records.find(r => r.id === 'SPC-09')?.value} :
      (D.meta.panel_checks.find(c => c.id === 'SPC-09') || {}).value || contract)}</div></div>
  </div>`;
  // histogram
  const mx = Math.max(...[0,1,2,3].map(k => Math.max(hF[k], hR[k])), 1); const Wd = Math.max(300, W_OF($('#mmHist', host)) - 36), Hh = 190, bw = (Wd - 60) / 4;
  const Y = v => 160 - Math.sqrt(v / mx) * 140;
  $('#mmHist', host).innerHTML = `<svg width="${Wd}" height="${Hh}" viewBox="0 0 ${Wd} ${Hh}">${[0,1,2,3].map(k => { const x0 = 50 + k * bw;
    return `<rect x="${x0 + bw*.14}" y="${Y(hF[k])}" width="${bw*.34}" height="${160 - Y(hF[k])}" fill="${css('--fwd')}" rx="2" data-tip="F · ${k} uyumsuzluk: ${hF[k]} yer"/><rect x="${x0 + bw*.52}" y="${Y(hR[k])}" width="${bw*.34}" height="${160 - Y(hR[k])}" fill="${css('--rev')}" rx="2" data-tip="R · ${k} uyumsuzluk: ${hR[k]} yer"/>
      <text x="${x0 + bw*.31}" y="${Y(hF[k]) - 4}" class="axm" text-anchor="middle">${hF[k]}</text><text x="${x0 + bw*.69}" y="${Y(hR[k]) - 4}" class="axm" text-anchor="middle">${hR[k]}</text>
      <text x="${x0 + bw/2}" y="178" class="lab" text-anchor="middle">${k === 0 ? 'tam eşleşme' : k + ' uyumsuz'}</text>`; }).join('')}<line x1="46" x2="${Wd-6}" y1="160" y2="160" stroke="${css('--line2')}"/><text x="6" y="20" class="ax">yer (√)</text></svg>
    <div class="legend"><span><i class="dot" style="background:var(--fwd)"></i>Forward</span><span><i class="dot" style="background:var(--rev)"></i>Reverse</span><span class="faint">Yükseklik karekök ölçekli</span></div>`;
  // ideogram
  const chrs = Object.keys(CHRLEN); const iw = Math.max(320, W_OF($('#ideo', host)) - 36), maxL = CHRLEN.chr1; const rowH = 15;
  const sites = [...p.bindings.F.top.map(s => ({...s, role:'F'})), ...p.bindings.R.top.map(s => ({...s, role:'R'}))];
  $('#ideo', host).innerHTML = `<svg width="${iw}" height="${chrs.length * rowH + 20}" viewBox="0 0 ${iw} ${chrs.length * rowH + 20}">${chrs.map((c, i) => { const y = 8 + i * rowH, w = (iw - 60) * CHRLEN[c] / maxL;
    return `<text x="40" y="${y + 8}" class="axm" text-anchor="end">${c.replace('chr','')}</text><rect x="48" y="${y + 2}" width="${w}" height="7" rx="3.5" fill="${css('--sunk')}" stroke="${css('--line')}"/>` + sites.filter(s => s.chrom === c).map(s => { const x = 48 + (s.start / CHRLEN[c]) * w; const hit = s.mm === 0 && s.chrom === pm.chrom && s.strand === (s.role === 'F' ? '+' : '-') && s.start === pm.tplOff + (s.role === 'F' ? pm.fStart : pm.rStart);
      return `<rect x="${x - 1.2}" y="${y - (s.role === 'F' ? 1 : -3)}" width="${hit ? 3.4 : 2.4}" height="${hit ? 13 : 8}" fill="${hit ? css('--ok') : s.role === 'F' ? css('--fwd') : css('--rev')}" opacity="${hit ? 1 : s.mm <= 1 ? .95 : s.mm === 2 ? .7 : .4}" data-tip="${esc(`<b>${s.role} · ${s.chrom}:${(s.start+1).toLocaleString('tr-TR')} (${s.strand})</b><br>${s.mm} uyumsuzluk${hit ? ' · HEDEF' : ''}<br><span class='m'>${s.mask}</span>`)}"/>`; }).join(''); }).join('')}</svg>
    <div class="legend"><span><i class="dot" style="background:var(--ok)"></i>hedef (${esc(pm.chrom)})</span><span><i class="dot" style="background:var(--fwd)"></i>F yeri</span><span><i class="dot" style="background:var(--rev)"></i>R yeri</span><span class="faint">soluk = daha çok uyumsuzluk</span></div>`;
  const drawAln = () => { const role = STATE.specRole || 'F', mm = STATE.specMM ?? 9; const primer = role === 'F' ? pm.F : pm.R; const list = p.bindings[role].top.filter(s => s.mm <= mm);
    const tgt = role === 'F' ? pm.fStart : pm.rStart;
    const head = siteAlignHTML(primer, primer, contract.anchor_length);
    $('#alnList', host).innerHTML = `<div class="siterow" style="background:var(--sunk)"><div class="loc"><b>${role === 'F' ? 'Forward' : 'Reverse'} primer</b> 5′→3′</div><div class="aln">${head.top}</div><div class="tiny faint">3′ çapa: son ${contract.anchor_length} nt</div></div>` +
      list.map(s => { const al = siteAlignHTML(primer, s.site, contract.anchor_length); const hit = s.mm === 0 && s.chrom === pm.chrom && s.strand === (role === 'F' ? '+' : '-') && s.start === pm.tplOff + tgt;
        return `<div class="siterow ${hit ? 'hit' : ''}"><div class="loc"><b>${s.chrom}:${(s.start + 1).toLocaleString('tr-TR')}</b> (${s.strand})${hit ? ' · <span style="color:var(--ok);font-weight:700">hedef</span>' : ''}</div><div class="aln">${al.bot}</div><div class="row" style="gap:6px;justify-content:flex-end">${s.mm === 0 ? sev(hit ? 'PASS' : 'FAIL', hit ? 'tam · hedef' : 'tam eşleşme') : `<span class="chip mono">${s.mm} mm</span>`}${s.anchored ? '<span class="tag">3′ sağlam</span>' : ''}</div></div>`; }).join('') +
      `<div class="siterow"><div class="tiny faint" style="grid-column:1/-1">Toplam ${p.bindings[role].n} yerden en yakın ${p.bindings[role].top.length} tanesi pakette. ${scan ? 'En düşük uyumsuzluklu yerler gösteriliyor.' : 'Kaynak: kayıtlı proje paneli.'}</div></div>`; };
  bindPairActions(host, pm);
  drawAln();
  $$('#roleSeg button', host).forEach(b => b.onclick = () => { STATE.specRole = b.dataset.role; $$('#roleSeg button', host).forEach(x => x.setAttribute('aria-pressed', x === b)); drawAln(); });
  $$('#mmSeg button', host).forEach(b => b.onclick = () => { STATE.specMM = +b.dataset.mm; $$('#mmSeg button', host).forEach(x => x.setAttribute('aria-pressed', x === b)); drawAln(); });
}
function offtHTML(pm, p, contract, summary) {
  if (!p.offtarget.length) return `<div class="row" style="gap:10px">${sev('PASS', 'yok')}<span class="small muted">Arama sözleşmesi içinde hedef dışı ürün raporlanmadı. Beklenen hedef ürün sayısı ${summary.products_intended}.</span></div>`;
  const maskRow = (primer, mask) => primer.split('').map((c, i) => mask[i] === 'X' ? `<span class="mm">${c}</span>` : '<span>.</span>').join('');
  const findSite = (role, chrom, pos, strand) => p.bindings[role].top.find(s => s.chrom === chrom && s.start === pos && s.strand === strand);
  return p.offtarget.map(o => { const lr = o.left_role || 'F', rr = o.right_role || 'R', single = lr === rr;
    const sf = findSite(lr, o.chrom, o.left_position ?? o.amplicon_start, '+'), sr = findSite(rr, o.chrom, o.right_position ?? o.amplicon_start + o.amplicon_length - pm[rr].length, '-');
    const fa = sf ? siteAlignHTML(pm[lr], sf.site, contract.anchor_length).bot : maskRow(pm[lr], o.left_mask), ra = sr ? siteAlignHTML(pm[rr], sr.site, contract.anchor_length).bot : maskRow(pm[rr], o.right_mask);
    return `<div style="padding:12px 18px;border-bottom:1px solid var(--line)"><div class="row" style="gap:8px;margin-bottom:8px"><b class="mono">${o.chrom}:${(o.amplicon_start+1).toLocaleString('tr-TR')}</b><span class="chip"><b>${o.amplicon_length}</b> bp</span><span class="chip">${o.total_mismatches} toplam uyumsuzluk</span><span class="chip" title="sol site + iplikte, sağ site − iplikte">${lr}→ ←${rr}${single ? ' · tek primer' : ''}</span><span class="tag">${esc(o.class || (o.total_mismatches === 0 ? 'tam eşleşmeli' : 'uyumsuzluklu'))}</span>${o.left_anchor_intact === 'True' || o.left_anchor_intact === true ? `<span class="tag">sol ${lr} 3′ sağlam</span>` : ''}${o.right_anchor_intact === 'True' || o.right_anchor_intact === true ? `<span class="tag">sağ ${rr} 3′ sağlam</span>` : ''}<span class="sp" style="flex:1"></span>${o.total_mismatches === 0 ? sev('FAIL','tam hedef dışı') : sev('WARN','hedef dışı')}</div>
      <div class="aln small"><span class="lbl">${lr}  5′ </span>${seqC(pm[lr])}<span class="lbl"> 3′</span>\n<span class="lbl">   yer </span>${fa}\n<span class="lbl">${rr}  5′ </span>${seqC(pm[rr])}<span class="lbl"> 3′</span>\n<span class="lbl">   yer </span>${ra}</div></div>`; }).join('');
}
function kvHTML(obj) { return `<div class="kv">${Object.entries(obj || {}).filter(([k]) => !technicalField(k)).map(([k, v]) => `<div>${esc(tidyKey(k))}</div><div class="mono">${esc(typeof v === 'object' ? evidenceJson(v) : v)}</div>`).join('')}</div>`; }
function ckCard(id, e) {
  const c = CATI[id] || {name: id}; if (!e) return '';
  const recs = (e.recs || []).filter(r => !r.ne);
  const val = recs.map(r => `${r.role && r.role !== 'pair' && r.role !== 'panel' ? `<b>${r.role}</b> ` : ''}${shortVal(id, r.value)}`).join('<br>');
  return `<button class="checkcard" data-rule="${id}" style="text-align:left;cursor:pointer"><div class="top"><span class="id">${id}</span>${e.sev ? sev(e.sev) : sev(e.status === 'ne' ? 'NA' : 'PEND', e.status === 'ne' ? 'koşulmadı' : 'bekliyor')}</div><div class="tt">${esc(c.name)}</div><div class="val small">${val || `<span class="faint">${esc(e.reason || '')}</span>`}</div></button>`;
}
function shortVal(id, v) {
  if (v == null) return '—'; if (typeof v !== 'object') return esc(v);
  const pick = {'SPC-01':['exact_match_sites'],'SPC-02':['exact_sites','anchored_mismatch_sites'],'SPC-03':['products_total','offtarget_total','offtarget_exact'],'SPC-05':['terminal3_tolerant','terminal3_disruptive'],'SPC-06':['pseudogene_overlap_sites'],'SPC-07':['segdup_overlap_sites'],'SPC-04':['evaluated','n_unscored'],
    'TH-01':['tm_c'],'TH-02':['spread_c'],'TH-03':['theta_at_ta'],'TH-04':['dG37_last_5_nt_kcal_mol'],'TH-05':['theta_at_tm'],'DIM-01':['dG37_kcal_mol'],'DIM-02':['dG37_kcal_mol'],'DIM-03':['dG37_kcal_mol'],'DIM-04':['dG37_kcal_mol'],
    'VAR-01':['n_variants','max_af'],'VAR-02':['n_variants','max_af'],'VAR-03':['n_common','max_af'],'VAR-04':['n_scored','n_unscored'],'VAR-06':['n_indels','length_change'],'VAR-07':['n_population_specific'],'VAR-08':['verdict','checks_triggered']}[id];
  const ks = pick || Object.keys(v).filter(k => typeof v[k] !== 'object').slice(0, 2);
  return ks.filter(k => k in v).map(k => `<span class="faint">${esc(tidyKey(k))}</span> <span class="mono">${esc(fmtv(v[k]))}</span>`).join(' · ') || '<span class="faint">ayrıntı →</span>';
}
const fmtv = x => typeof x === 'number' ? (Math.abs(x) > 0 && Math.abs(x) < 1e-3 ? x.toExponential(2) : Number.isInteger(x) ? x : +x.toFixed(4)) : typeof x === 'boolean' ? (x ? 'evet' : 'hayır') : x;

/* ---- varyant (SNPCheck) ---- */
function tabVar(host, pm, E) {
  const liveVar = pm.kind === 'live' || !!pm.liveVariants;
  if (liveVar && !pm.liveVariants?.evaluated) {
    const result = pm.liveVariants;
    const retry = liveReady() ? `<div class="row live-actions">
      ${pairAction(pm, 'variants', 'Varyant verisini yeniden sorgula')}</div>${jobNotice('variants:' + pm.key)}` : '';
    host.innerHTML = `<div class="stack">${result ? renderLocalOutput('variants', result) : liveNote('Bu çift için gnomAD kaydı henüz alınmadı.')}${retry}</div>`;
    bindPairActions(host, pm);
    return;
  }
  if (!pm.hasVar) { host.innerHTML = `<div class="verdict na"><div class="ic">${ic('snp')}</div><div><h3>Bu şablon için varyant verisi yok</h3><div class="small muted">gnomAD v4 önbelleği yalnız proje bölgesini (chr22:20.000.001–20.001.000) kapsar. Başka bölgeler için <code>RUN_VARIANT.cmd</code> ile varyant önbelleği oluşturun.</div></div></div>`; return; }
  const vr = ['VAR-01','VAR-02','VAR-03','VAR-04','VAR-06','VAR-07','VAR-08'];
  const vsev = worst(vr.map(id => E[id] && E[id].sev).filter(Boolean));
  const ampV = varsIn(pm.fEnd + 1, pm.rStart - 1, pm);
  host.innerHTML = `<div class="stack">
    <div class="verdict ${vsev === 'FAIL' ? 'fail' : vsev === 'WARN' ? 'warn' : 'ok'}"><div class="ic">${ic(vsev === 'INFO' ? 'ok' : 'warn')}</div><div><h3>${vsev === 'FAIL' ? 'Yaygın varyant primer altında' : vsev === 'WARN' ? 'Primer altında incelenmesi gereken varyant var' : 'Primer altında anlamlı varyant yok'}</h3><div class="small">gnomAD v4 genom (r4) · ${varsIn(pm.fStart, pm.rEnd, pm).length} varyant amplikonda · ${varsIn(pm.fStart, pm.fEnd, pm).length} F altında · ${varsIn(pm.rStart, pm.rEnd, pm).length} R altında</div></div><span class="sp" style="flex:1"></span>${sev(vsev)}</div>
    ${['F','R'].map(role => `<div class="card"><div class="hd"><h2>${role === 'F' ? 'Forward' : 'Reverse'} primer bağlanma yeri</h2><span class="faint small">${pm.chrom}:${gposAny(pm, footprint(pm, role)[0])}–${gposAny(pm, footprint(pm, role)[1])} · ${role === 'F' ? '+' : '−'} iplik</span><span class="sp"></span>${(pm.kind === 'project' ? ['VAR-01','VAR-02','VAR-03'] : ['VAR-01','VAR-02','VAR-03']).map(id => { const r = (E[id]?.recs || []).find(x => x.role === role); return r ? `<span class="chip">${id} ${sev(r.sev)}</span>` : ''; }).join('')}</div>
      <div class="bd"><div class="scrollx" id="snp${role}"></div><div class="g2" style="margin-top:14px"><div id="lol${role}"></div><div class="scrollx" id="vt${role}"></div></div></div></div>`).join('')}
    <div class="card"><div class="hd"><h2>Motor kayıtları</h2><span class="faint small">variant_qc.py · VAR ailesi (kodun kural numaralandırması)</span></div><div class="bd"><div class="ccgrid">${vr.map(id => E[id] ? ckCard(id, E[id]) : '').join('')}</div></div></div>
    ${liveVar ? liveResult('gnomAD · amplikon ve primer bölgeleri', `<div class="scrollx" id="liveVarTrack"></div><p class="small muted">${pm.liveVariants.variants.length ? 'İşaretlerin üzerine gelerek alel ve frekans bilgisine bakabilirsiniz.' : 'Sorgu tamamlandı. Bu amplikonda kayıtlı varyant bulunmadı.'}</p>`) : ''}
    <div class="card"><div class="hd"><h2>Primerler arası varyantlar</h2><span class="faint small">amplifikasyonu engellemez; prob/blocker tasarımı için bilgi</span></div><div class="bd small">${ampV.length ? `${ampV.length} varyant · AF ≥ %0,1: ${ampV.filter(v => v.af >= 0.001).map(v => `<span class="chip mono">${esc(v.rsid)} ${esc(v.ref)}>${esc(v.alt)} ${afFmt(v.af)}</span>`).join(' ') || 'yok'}` : 'Yok'}</div></div>
    ${liveVar ? liveResult('Amplikonla örtüşen varyantlar', localRows(pm.liveVariants.variants,
      [['pos','Konum'],['rsid','rsID'],['ref','Referans'],['alt','Alternatif'],['genome_af','Genom AF'],['exome_af','Ekzom AF'],['af','Kullanılan AF']])) : ''}
  </div>`;
  const liveTrack = $('#liveVarTrack', host);
  if (liveTrack) liveTrack.innerHTML = overviewTrack(Math.max(345, W_OF(liveTrack)), {tpl:pm.tpl, off:pm.tplOff,
    pairs:[pm], win:[pm.fStart, pm.rEnd + 1], showGC:false, variantPair:pm, interactive:false});
  for (const role of ['F','R']) {
    const [lo, hi] = footprint(pm, role); const flank = 6; const a = Math.max(0, lo - flank), b = Math.min(pm.tpl.length - 1, hi + flank);
    const vs = varsIn(a, b, pm); const primer = role === 'F' ? pm.F : pm.R;
    // SNPCheck kutuları: primer 5'->3' yönünde göster
    const idxs = []; for (let i = a; i <= b; i++) idxs.push(i); if (role === 'R') idxs.reverse();
    const baseAt = i => role === 'F' ? pm.tpl[i] : OS.COMP[pm.tpl[i]];
    const BW = 30, LM = 70; const W = LM + idxs.length * BW + 20; const byIdx = {}; vs.forEach(v => (byIdx[v.idx] = byIdx[v.idx] || []).push(v));
    const maxStack = Math.max(1, ...Object.values(byIdx).map(x => x.length)); const H = 92 + Math.min(maxStack, 6) * 22 + 10;
    const o = [`<text x="8" y="36" class="ax">${role === 'F' ? '5′' : '5′'} primer</text><text x="8" y="70" class="ax">3′ mesafe</text>`];
    idxs.forEach((i, k) => { const x = LM + k * BW; const inP = i >= lo && i <= hi; const d = dist3(pm, role, i); const c = baseAt(i);
      o.push(`<rect x="${x+1}" y="20" width="${BW-2}" height="28" rx="5" fill="${inP ? (d < 5 ? css('--rev-soft') : css('--accent-soft')) : css('--sunk')}" stroke="${inP && d < 5 ? css('--rev') : 'none'}" stroke-width="1"/><text x="${x + BW/2}" y="39" text-anchor="middle" class="seqt" style="font-size:14px;font-weight:${inP?700:400};fill:${css('--n' + c) || css('--ink')};opacity:${inP?1:.5}">${c}</text>`);
      if (inP) o.push(`<text x="${x + BW/2}" y="68" text-anchor="middle" class="axm" style="fill:${d < 5 ? css('--rev') : css('--ink3')};font-weight:${d<5?700:400}">${d}</text>`);
      (byIdx[i] || []).sort((p, q) => q.af - p.af).slice(0, 6).forEach((v, j) => { const y = 82 + j * 22; const col = varColor(v.af);
        o.push(`<rect x="${x+2}" y="${y}" width="${BW-4}" height="18" rx="4" fill="${col}" opacity="${v.af >= 0.001 ? .95 : .35}" data-tip="${esc(varTip(v) + `<br>3′ ucuna ${d} nt`)}"/><text x="${x + BW/2}" y="${y + 13}" text-anchor="middle" style="font:700 9.5px var(--mono);fill:${v.af >= 0.001 ? css('--panel') : css('--ink')};pointer-events:none">${v.ref.length === 1 && v.alt.length === 1 ? (role === 'F' ? v.alt : OS.COMP[v.alt]) : (v.ref.length > v.alt.length ? 'del' : 'ins')}</text>`); }); });
    const p3 = role === 'F' ? LM + (idxs.indexOf(hi) + 1) * BW : LM + (idxs.indexOf(lo) + 1) * BW;
    o.push(`<text x="${p3 - 2}" y="14" text-anchor="end" style="font:700 10px var(--sans);fill:${css('--rev')}">3′ ucu · son 5 nt (VAR-02 penceresi)</text>`);
    $('#snp' + role, host).innerHTML = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" role="img" aria-label="SNPCheck görünümü">${o.join('')}</svg>`;
    // lolipop: konum vs log AF
    const inV = varsIn(lo, hi, pm); const Lw = Math.max(280, W_OF($('#lol' + role, host)));
    const pts = inV.map(v => ({x: primer.length - 1 - dist3(pm, role, v.idx), y: Math.log10(Math.max(v.af || 1e-7, 1e-7)), v}));
    const m = {l:46, r:10, t:12, b:30}, iw = Lw - m.l - m.r, ih = 170 - m.t - m.b; const X = k => m.l + (k + .5) / primer.length * iw, Yv = y => m.t + ih - (y + 6.5) / 6.5 * ih;
    const g = [`<rect x="${X(primer.length - 5.5)}" y="${m.t}" width="${X(primer.length - .5) - X(primer.length - 5.5)}" height="${ih}" fill="${css('--rev')}" opacity=".08"/>`];
    for (const [t, lab] of [[-6,'10⁻⁶'],[-4,'10⁻⁴'],[-3,'0,1 %'],[-2,'1 %'],[0,'100 %']]) g.push(`<line x1="${m.l}" x2="${Lw-m.r}" y1="${Yv(t)}" y2="${Yv(t)}" stroke="${t === -3 || t === -2 ? css('--warn') : css('--grid')}" ${t === -3 || t === -2 ? 'stroke-dasharray="4 3" opacity=".7"' : ''}/><text x="${m.l-6}" y="${Yv(t)+3}" class="axm" text-anchor="end">${lab}</text>`);
    primer.split('').forEach((c, k) => g.push(`<text x="${X(k)}" y="${m.t + ih + 14}" text-anchor="middle" class="seqt" style="font-size:9.5px;fill:${css('--n' + c)}">${c}</text>`));
    for (const p of pts) g.push(`<line x1="${X(p.x)}" x2="${X(p.x)}" y1="${Yv(-6.5)}" y2="${Yv(p.y)}" stroke="${varColor(p.v.af)}" stroke-width="1.5"/><circle cx="${X(p.x)}" cy="${Yv(p.y)}" r="${p.v.af >= 0.001 ? 5 : 3.5}" fill="${varColor(p.v.af)}" stroke="${css('--panel')}" stroke-width="1.5" data-tip="${esc(varTip(p.v))}"/>`);
    $('#lol' + role, host).innerHTML = `<div class="sub-h">Alel frekansı · primer boyunca (5′→3′)</div><svg width="${Lw}" height="170" viewBox="0 0 ${Lw} 170">${g.join('')}<text x="${m.l}" y="${m.t-2}" class="ax">AF (log)</text></svg>`;
    // tablo + popülasyon şeridi
    const POPS = ['afr','amr','asj','eas','fin','mid','nfe','sas','ami','remaining'];
    $('#vt' + role, host).innerHTML = `<div class="sub-h">Varyantlar · ${inV.length}</div>` + (inV.length ? `<table class="t" style="min-width:480px"><thead><tr><th>rsID</th><th>Alel</th><th class="n">3′</th><th class="n">AF</th><th>Popülasyonlar</th></tr></thead><tbody>${inV.sort((p, q) => q.af - p.af).map(v => `<tr><td class="mono small"><b>${esc(v.rsid)}</b><div class="tiny faint">${vtype(v)} · ${v.pos.toLocaleString('tr-TR')}</div></td><td class="mono small">${esc(v.ref.length > 8 ? v.ref.slice(0,6)+'…' : v.ref)}&gt;${esc(v.alt.length > 8 ? v.alt.slice(0,6)+'…' : v.alt)}</td><td class="n">${dist3(pm, role, v.idx)}</td><td class="n" style="color:${varColor(v.af)};font-weight:${v.af>=0.001?700:400}">${afFmt(v.af)}</td><td><div style="display:flex;gap:2px">${POPS.map(k => { const a = (v.pop || {})[k] || 0; const op = a <= 0 ? 0 : Math.min(1, .18 + (Math.log10(a) + 6) / 6);
        return `<i style="width:13px;height:13px;border-radius:3px;background:${a >= 0.01 ? css('--fail') : a >= 0.001 ? css('--warn') : css('--accent')};opacity:${a <= 0 ? .08 : op}" data-tip="${k}: ${afFmt(a)}"></i>`; }).join('')}</div></td></tr>`).join('')}</tbody></table><div class="tiny faint" style="margin-top:6px">Şerit sırası: ${POPS.join(' · ')}</div>` : '<div class="faint small">Bağlanma yerinde gnomAD varyantı yok.</div>');
  }
}

/* ---- yapı & dimer ---- */
function tabStruct(host, pm, E) {
  const motor = pm.kind === 'live' || !!pm.liveQC;
  if (motor && !pm.liveQC) { host.innerHTML = liveNote('Bu çift için yapı ve dimer QC kaydı bulunamadı.'); return; }
  const hasStructureRecords = motor || pm.kind === 'project';
  const recVal = (id, role) => { const r = (E[id]?.recs || []).find(x => x.role === role && !x.ne); return r ? r.value : null; };
  const dCard = (title, sub, anyV, p3V, list, sevs) => `<div class="card"><div class="hd"><h2>${title}</h2><span class="faint small">${sub}</span><span class="sp"></span>${sevs}</div><div class="bd">
      <div class="sub-h">En kararlı yapı (her konum)</div>${anyV && anyV.alignment ? '<div class="scrollx">' + dimerHTML(anyV.alignment) + '</div>' : '<div class="faint small">Kararlı yapı yok</div>'}
      ${anyV && anyV.alignment ? `<div class="row small" style="gap:14px;margin-top:8px"><span>ΔG₃₇ <b class="mono">${nf(anyV.dG37_kcal_mol,2)}</b> kcal/mol</span><span>ΔH <b class="mono">${nf(anyV.dH_kcal_mol,1)}</b></span><span>ΔS <b class="mono">${nf(anyV.dS_cal_mol_K,1)}</b></span><span>helis <b class="mono">${anyV.helix_len}</b> bp</span>${anyV.mismatch_count ? '<span class="tag">1 iç uyumsuzluk</span>' : ''}${anyV.bulge_count ? '<span class="tag">1 bulge</span>' : ''}</div>` : ''}
      <div class="sub-h" style="margin-top:14px">3′ uçlu, uzatılabilir yapı</div>${p3V && p3V.alignment ? '<div class="scrollx">' + dimerHTML(p3V.alignment) + '</div>' + `<div class="small" style="margin-top:6px">ΔG₃₇ <b class="mono">${nf(p3V.dG37_kcal_mol,2)}</b> kcal/mol · helis ${p3V.helix_len} bp</div>` : '<div class="faint small">3′ ucu eşleşen yapı yok — uzama riski düşük</div>'}
      ${list && list.length ? `<details style="margin-top:12px"><summary class="small" style="cursor:pointer;font-weight:600">Diğer olası yapılar (${list.length})</summary>${list.map(d => `<div style="border-top:1px dashed var(--line);padding:8px 0">${dimerHTML(d.alignment)}<div class="tiny faint">ΔG₃₇ ${nf(d.dG,2)} kcal/mol · ${d.helix_len} bp${d.three_prime_paired ? ' · 3′ eşleşmiş' : ''}</div></div>`).join('')}</details>` : ''}
    </div></div>`;
  const sevPair = (a, b) => [a, b].map((id, k) => { const r = E[id]; return r ? `<span class="chip">${id} ${sev(r.recs.filter(x => x.role === (k === 0 ? 'X' : 'X')).length ? 'INFO' : r.sev)}</span>` : ''; }).join('');
  const sevRole = (ids, role) => ids.map(id => { const r = (E[id]?.recs || []).find(x => x.role === role && !x.ne); return r ? `<span class="chip">${id} ${sev(r.sev)}</span>` : ''; }).join('');
  const alts = (a, b, best) => OS.dimerList(a, b, 5).filter(d => !best || Math.abs(d.dG - best.dG37_kcal_mol) > 1e-6).slice(0, 4);
  const hetA = recVal('DIM-03', 'pair'), het3 = recVal('DIM-04', 'pair');
  host.innerHTML = `<div class="stack">
    <div class="note">${ic('info')}<div>Dimer modeli: <b>${hasStructureRecords ? 'Analiz kaydı' : 'motor portu'}</b> · WC yığınları + en fazla bir iç uyumsuzluk (Allawi-Peyret) ya da tek bazlık bulge (SantaLucia-Hicks 2004; A·T kapanış sabiti ${D.engine_params.at_bulge_closure_dg37} kcal/mol, proje modelleme kararı). <span class="bd" style="color:var(--accent)">┃</span> G·C, <span style="color:var(--ink3)">│</span> A·T bağı, <span style="color:var(--fail)">×</span> uyumsuzluk, <span style="color:var(--fail)">^</span> bulge.</div></div>
    <div class="g3">
      ${dCard('Forward self-dimer', 'DIM-01 · DIM-02', recVal('DIM-01','F'), recVal('DIM-02','F'), alts(pm.F, pm.F, recVal('DIM-01','F')), sevRole(['DIM-01','DIM-02'],'F'))}
      ${dCard('Reverse self-dimer', 'DIM-01 · DIM-02', recVal('DIM-01','R'), recVal('DIM-02','R'), alts(pm.R, pm.R, recVal('DIM-01','R')), sevRole(['DIM-01','DIM-02'],'R'))}
      ${dCard('F × R hetero-dimer', 'DIM-03 · DIM-04', hetA, het3, alts(pm.F, pm.R, hetA), sevRole(['DIM-03','DIM-04'],'pair'))}
    </div>
    <div class="g2">${['F','R'].map(role => { const s = role === 'F' ? pm.F : pm.R; const ss01 = recVal('SS-01', role); const hp = OS.hairpinStem(s); const ss02 = recVal('SS-02', role);
      const db = ss01 && ss01.dot_bracket ? ss01.dot_bracket : null;
      return `<div class="card"><div class="hd"><h2>${role === 'F' ? 'Forward' : 'Reverse'} · hairpin</h2><span class="faint small">SS-01 · SS-02</span><span class="sp"></span>${sevRole(['SS-01','SS-02'], role)}</div><div class="bd">
        ${hasStructureRecords && db ? `<div class="sub-h">Motor: Zuker MFE (fold.zuker_mfe)</div><div class="scrollx">${hairpinSVG(s, db)}</div><div class="small">MFE <b class="mono">${nf(ss01.mfe_kcal_mol,2)}</b> kcal/mol · nokta-parantez <span class="mono">${esc(db)}</span></div>` : ''}
        ${`
        <div class="sub-h" style="margin-top:12px">En güçlü sap adayı (align.hairpin_stem)</div>
        ${hp ? `<div class="row" style="gap:16px;align-items:flex-start"><div>${stemLoopSVG(s, hp)}</div><div class="small" style="display:grid;gap:4px;min-width:150px"><span>Sap ΔG₃₇ <b class="mono">${nf(hp.dG,2)}</b> kcal/mol</span><span>sap ${hp.helix_len} bp · ilmek ${hp.loop} nt</span><span>3′ ucu sapta: <b>${hp.three_prime_paired ? 'evet' : 'hayır'}</b></span><span class="faint tiny">Yalnız yığın enerjisi; ilmek cezası eklenmez. ΔG > 0 ise yapı kararsızdır.</span>${ss02 ? `<span>SS-02 sap ΔG <b class="mono">${nf(ss02.stem_dG37_kcal_mol,2)}</b> · uzatılabilir: ${ss02.three_prime_extendable ? 'evet' : 'hayır'}</span>` : ''}</div></div>` : '<div class="faint small">En az 2 bp sap yok</div>'}`}
      </div></div>`; }).join('')}</div>
    ${hasStructureRecords ? `<div class="g2"><div class="card"><div class="hd"><h2>Primer bazlarının erişilebilirliği</h2><span class="faint small">SS-03 · bazın eşleşmemiş olma olasılığı</span></div><div class="bd" id="ss03"></div></div>
      <div class="card"><div class="hd"><h2>Amplikon ikincil yapısı</h2><span class="faint small">SS-04 · tek iplik amplikon MFE</span><span class="sp"></span>${E['SS-04'] ? sev(E['SS-04'].sev) : ''}</div><div class="bd" id="ss04"></div></div></div>` : ''}
  </div>`;
  if (hasStructureRecords) {
    const host3 = $('#ss03', host); const w = Math.max(300, W_OF(host3) - 36);
    host3.innerHTML = ['F','R'].map(role => { const v = recVal('SS-03', role); if (!v || !v.p_unpaired_by_base) return ''; const s = role === 'F' ? pm.F : pm.R; const n = v.p_unpaired_by_base.length, bw = (w - 40) / n;
      return `<div class="sub-h">${role} · 3′ ortalaması ${nf(v.mean_p_unpaired_3p,3)}</div><svg width="${w}" height="96" viewBox="0 0 ${w} 96">${v.p_unpaired_by_base.map((p, i) => `<rect x="${34 + i*bw + 1}" y="${8 + (1-p)*60}" width="${bw-2}" height="${p*60}" rx="2" fill="${p < .5 ? css('--warn') : role === 'F' ? css('--fwd') : css('--rev')}" opacity=".85" data-tip="${role} ${i+1}. baz ${s[i] || ''}: p(eşleşmemiş) ${nf(p,3)}"/><text x="${34 + i*bw + bw/2}" y="84" text-anchor="middle" class="seqt" style="font-size:9px;fill:${css('--n'+(s[i]||'N'))}">${s[i] || ''}</text>`).join('')}<text x="2" y="12" class="axm">1,0</text><text x="2" y="70" class="axm">0</text><line x1="34" x2="${w}" y1="38" y2="38" stroke="${css('--grid')}" stroke-dasharray="3 3"/></svg>`; }).join('') + `<div class="tiny faint">Primer kendi içinde eşleştiğinde bu olasılık düşer (sarı &lt; 0,5).</div>`;
    const s4 = (E['SS-04']?.recs || []).find(r => Number.isFinite(r.value?.structure_tm_c)); const h4 = $('#ss04', host);
    if (s4) { const v = s4.value; h4.innerHTML = `<div class="mini-stats" style="grid-template-columns:repeat(3,1fr)"><div><div class="k">MFE</div><div class="v">${nf(v.amplicon_mfe_kcal_mol,2)}</div></div><div><div class="k">Tm (yapı)</div><div class="v">${nf(v.structure_tm_c,1)} °C</div></div><div><div class="k">Baz çifti</div><div class="v">${v.n_base_pairs}</div></div></div>
      ${gaugeHTML(v.structure_tm_c, 20, 90, [{a:20,b:50,c:css('--ok'),op:.3},{a:50,b:60,c:css('--warn'),op:.4},{a:60,b:90,c:css('--fail'),op:.35}], x => nf(x,0), ' °C')}
      <div class="small muted" style="margin-top:8px">Eşik: WARN Tm_yapı &gt; 50 °C, FAIL &gt; 60 °C · Ta ${v.anneal_c} °C. ${esc(s4.note || '')}</div>`; }
    else h4.innerHTML = '<p class="small muted">Amplikon ikincil yapısı için ölçüm yok.</p>';
  }
  if (motor) host.firstElementChild.insertAdjacentHTML('beforeend', qualityCheckCards(E, ['DIM-', 'SS-', 'AMP-']));
}

/* ---- termodinamik ---- */
function tabLiveThermo(host, pm, E) {
  const qc = pm.liveQC;
  if (!qc) { host.innerHTML = liveNote('Bu çift için termodinamik QC kaydı bulunamadı.'); return; }
  const get = (id, role) => E[id]?.recs.find(r => r.role === role && r.value)?.value;
  const roles = ['F', 'R'];
  const ta = get('TH-03', 'F')?.ta_c ?? pm.cond.ta;
  const tmValues = roles.map(role => qc.thermo?.[role]?.tm_c).filter(Number.isFinite);
  const series = roles.map(role => ({role, pts:qc.melting_curves?.[role] || [], c:css(role === 'F' ? '--fwd' : '--rev'), dash:role === 'R' ? '6 4' : undefined}));
  const canPlot = series.every(s => s.pts.length > 1);
  host.innerHTML = `<div class="stack">
    ${canPlot ? liveResult('Primer erime eğrileri · θ(T)', `
      <div class="row doc-legend"><span><i class="dot" style="background:var(--fwd)"></i>İleri primer</span>
      <span><i class="dot" style="background:var(--rev)"></i>Geri primer (kesikli)</span><span>Kesikli çizgi: Ta ${nf(ta,1)} °C</span></div>
      <div id="liveThermoPlot" class="scrollx"></div>
      <p class="small muted">θ(T), modelin o sıcaklıkta hedefe bağlı olduğunu öngördüğü primer oranıdır. Tm’de θ ≈ 0,5.</p>`) :
      liveNote('Bu QC kaydında eğri noktaları yok. Yeni bir tasarımda eğriler otomatik oluşturulur.')}
    <div class="g2">${roles.map(role => {
      const t = qc.thermo?.[role];
      const bound = get('TH-03', role);
      const salt = get('TH-02', role);
      return liveResult(role === 'F' ? 'İleri primer (F)' : 'Geri primer (R)', `
        <div class="mini-stats thermo-stats"><div><div class="k">Tm</div><div class="v">${nf(t?.tm_c,2)} °C</div></div>
        <div><div class="k">Bağlanma oranı · θ(Ta)</div><div class="v">${nf(bound?.theta_at_ta,3)}</div></div>
        <div><div class="k">Tuz modelleri arası Tm farkı</div><div class="v">${nf(salt?.spread_c,2)} °C</div></div></div>`);
    }).join('')}</div>
    ${qualityCheckCards(E, ['TH-', 'COND-', 'PAIR-'])}
  </div>`;
  if (!canPlot) return;
  const plot = $('#liveThermoPlot', host);
  const range = [...tmValues, ta].filter(Number.isFinite);
  const full = series[0].pts;
  const x0 = range.length ? Math.max(full[0][0], Math.floor(Math.min(...range) - 8)) : full[0][0];
  const x1 = range.length ? Math.min(full.at(-1)[0], Math.ceil(Math.max(...range) + 8)) : full.at(-1)[0];
  for (const s of series) s.pts = s.pts.filter(([x]) => x >= x0 && x <= x1);
  plot.innerHTML = lineChart(W_OF(plot), 245, series, {
    x0, x1, y0:0, y1:1, xlabel:'Sıcaklık (°C)', yticks:[0,.25,.5,.75,1], yfmt:y=>nf(y,2),
    vlines:Number.isFinite(ta) ? [{x:ta, t:'Ta', end:ta > (x0+x1)/2}] : [],
    hover:i => `<b>${nf(series[0].pts[i][0],1)} °C</b><br>${series.map(s => s.role + ' · θ ' + nf(s.pts[i]?.[1],3)).join('<br>')}`
  });
  $('svg', plot).setAttribute('aria-label', 'İleri ve geri primerin sıcaklığa göre bağlanma oranı');
}

function tabThermo(host, pm, E) {
  if (pm.kind === 'live' || pm.liveQC) { tabLiveThermo(host, pm, E); return; }
  const c = pm.cond; const tF = OS.thermo(pm.F, c, pm.salt), tR = OS.thermo(pm.R, c, pm.salt);
  const prof = PROFILES[pm.profile] || PROFILES.ctdna; const [tlo, thi] = prof.tm_window_c;
  const col = role => { const s = role === 'F' ? pm.F : pm.R, t = role === 'F' ? tF : tR, sp = OS.saltSpread(s, c), d3 = OS.dG3(s), hv = OS.halves(s);
    const rs = id => { const r = (E[id]?.recs || []).find(x => x.role === role && !x.ne); return r ? sev(r.sev) : ''; };
    const models = [['owczarzy_2004','Owczarzy 2004'],['owczarzy_2008','Owczarzy 2008'],['von_ahsen_2001','von Ahsen 2001']]; const mn = Math.min(sp.owczarzy_2004, sp.owczarzy_2008, sp.von_ahsen_2001) - 1.5, mxx = Math.max(sp.owczarzy_2004, sp.owczarzy_2008, sp.von_ahsen_2001) + 1.5;
    return `<div class="card"><div class="hd"><h2>${role === 'F' ? 'Forward' : 'Reverse'}</h2><span class="seq small">${seqC(s)}</span></div><div class="bd" style="display:grid;gap:16px">
      <div><div class="row" style="justify-content:space-between"><b class="small">TH-01 · Tm (NN, ${esc(pm.salt || 'owczarzy_2008')})</b>${rs('TH-01')}</div><div style="font:700 26px var(--display)">${nf(t.tm,2)} °C</div>
        ${gaugeHTML(t.tm, tlo - 6, thi + 6, [{a:tlo,b:thi,c:css('--ok'),op:.35}], x => nf(x,0), ' °C')}<div class="tiny faint">Profil penceresi ${tlo}–${thi} °C · Tm(1 M) ${nf(t.tm1M,2)} °C · rejim ${esc(t.regime)}</div></div>
      <div><div class="row" style="justify-content:space-between"><b class="small">TH-02 · tuz modeli duyarlılığı</b>${rs('TH-02')}</div>
        <div style="position:relative;height:58px;margin:4px 6px 0">${'<div style="position:absolute;left:0;right:0;top:27px;height:2px;background:var(--line2)"></div>'}${models.map(([k, lab], i) => { const pc = (sp[k] - mn) / (mxx - mn) * 100; const on = sp.used.includes(k);
          return `<div style="position:absolute;left:${pc}%;top:18px;width:4px;height:20px;margin-left:-2px;border-radius:2px;background:${on ? 'var(--accent)' : 'var(--ink3)'}"></div><div class="tiny mono" style="position:absolute;left:${pc}%;top:${i === 1 ? 40 : 0}px;transform:translateX(${pc < 15 ? '0' : pc > 85 ? '-100%' : '-50%'});white-space:nowrap;color:${on ? 'var(--ink)' : 'var(--ink3)'}">${lab} · ${nf(sp[k],1)}</div>`; }).join('')}</div>
        <div class="tiny faint">Karşılaştırılan modeller arası yayılım <b class="mono">${nf(sp.spread,2)} °C</b> (WARN &gt; 2, FAIL &gt; 4) · TH-02/TH-04 eşik karar desteği Doğrulama sayfasında</div></div>
      <div class="g2" style="gap:12px">
        <div><div class="row" style="justify-content:space-between"><b class="small">TH-03 · θ(Ta)</b>${rs('TH-03')}</div><div style="font:700 22px var(--display)">${nf(t.thetaTa,3)}</div>${gaugeHTML(t.thetaTa, 0, 1, [{a:0,b:.8,c:css('--fail'),op:.25},{a:.8,b:.95,c:css('--warn'),op:.3},{a:.95,b:1,c:css('--ok'),op:.35}], x => nf(x,1))}<div class="tiny faint">ΔG(Ta) ${nf(t.dGta,2)} kcal/mol · primer fazlası</div></div>
        <div><div class="row" style="justify-content:space-between"><b class="small">TH-04 · 3′ son 5 nt ΔG</b>${rs('TH-04')}</div><div style="font:700 22px var(--display)">${nf(d3,2)}</div>${gaugeHTML(d3, -13, -3, [{a:-13,b:-11,c:css('--fail'),op:.3},{a:-11,b:-9,c:css('--warn'),op:.3},{a:-9,b:-6,c:css('--ok'),op:.35},{a:-6,b:-3,c:css('--warn'),op:.3}], x => nf(x,0))}<div class="tiny faint">INFO bandı −9…−6 kcal/mol</div></div>
      </div>
      <div class="kv"><div>ΔH (kcal/mol)</div><div class="mono">${nf(t.dH,2)}</div><div>ΔS 1 M (cal/mol·K)</div><div class="mono">${nf(t.dS,2)}</div><div>ΔS etkin</div><div class="mono">${nf(t.dS_eff,3)}</div><div>ΔG₃₇ 1 M</div><div class="mono">${nf(t.dG37,3)}</div><div>CMP-07 ΔG 5′ / 3′ yarı</div><div class="mono">${nf(hv.d5,2)} / ${nf(hv.d3,2)}</div><div>GC · uzunluk</div><div class="mono">%${nf(OS.gcf(s)*100,1)} · ${s.length} nt</div></div>
    </div></div>`; };
  const prodTm = pm.kind === 'project' ? (pm.raw.seqpair.find(x => x.id === 'PAIR-04') || {}).value : null;
  const det = pm.kind === 'project' ? pm.raw.detect : [];
  host.innerHTML = `<div class="stack"><div class="g2">${col('F')}${col('R')}</div>
    <div class="card"><div class="hd"><h2>Erime eğrisi θ(T)</h2><span class="faint small">TH-05 · iki-hal modeli · motorun concentration.build portu</span><span class="sp"></span><div class="legend"><span><i class="dot" style="background:var(--fwd)"></i>F</span><span><i class="dot" style="background:var(--rev)"></i>R</span></div></div><div class="bd" id="melt"></div></div>
    <div class="g2">
      <div class="card"><div class="hd"><h2>Ürün ve önerilen Ta</h2><span class="faint small">PAIR-04 · TH-07</span></div><div class="bd">${prodTm ? `<div class="mini-stats" style="grid-template-columns:repeat(3,1fr)"><div><div class="k">Tm ürün</div><div class="v">${nf(prodTm.product_tm_c,1)}</div></div><div><div class="k">fark</div><div class="v">${nf(prodTm.gap_c,1)}</div></div><div><div class="k">Ta Rychlik</div><div class="v">${nf(0.3*Math.min(pm.tmF,pm.tmR) + 0.7*prodTm.product_tm_c - 14.9,1)}</div></div></div><div class="small muted" style="margin-top:10px">Ta = 0,3·Tm(primer) + 0,7·Tm(ürün) − 14,9 (Rychlik 1990). Kullanılan Ta ${c.ta} °C; θ-tabanlı bakışta Ta'da θ ≈ ${nf(Math.min(tF.thetaTa,tR.thetaTa),2)}.</div>` : '<div class="faint small">Ürün Tm değeri motor amplikon modülünden gelir; bu çift için mevcut değil.</div>'}</div></div>
      <div class="card"><div class="hd"><h2>ctDNA tespit edilebilirliği</h2><span class="faint small">CT-02 · CT-03 · Poisson</span></div><div class="bd" id="det">${det.length ? '' : '<div class="faint small">Yalnız ctDNA profilindeki proje çiftleri için.</div>'}</div></div>
    </div></div>`;
  const w = Math.max(320, W_OF($('#melt', host)) - 36); const T = []; for (let x = 40; x <= 80.001; x += 0.25) T.push(x);
  const sF = T.map(x => [x, tF.theta(x)]), sR = T.map(x => [x, tR.theta(x)]);
  $('#melt', host).innerHTML = lineChart(w, 250, [{pts: sF, c: css('--fwd'), w: 2.4}, {pts: sR, c: css('--rev'), w: 2.4}], {x0: 40, x1: 80, y0: 0, y1: 1, yticks: [0, .25, .5, .75, 1], ylabel: 'θ (bağlı kesir)', xlabel: 'sıcaklık °C',
    vlines: [{x: c.ta, t: 'Ta ' + c.ta + ' °C', c: css('--ink2')}, {x: tF.tm, t: 'Tm F ' + nf(tF.tm,1), c: css('--fwd'), dy: 16, end: true}, {x: tR.tm, t: 'Tm R ' + nf(tR.tm,1), c: css('--rev'), dy: 32, end: true}],
    points: [{x: c.ta, y: tF.thetaTa, c: css('--fwd')}, {x: c.ta, y: tR.thetaTa, c: css('--rev')}],
    hover: i => `<b>${nf(T[i],2)} °C</b><br>θ F ${nf(sF[i][1],3)} · θ R ${nf(sR[i][1],3)}`});
  if (det.length) { const d0 = det[0]; const dw = Math.max(300, W_OF($('#det', host)) - 36); const pts = []; for (let e = -5; e <= -1.001; e += 0.05) { const vaf = Math.pow(10, e); pts.push([e, 1 - Math.exp(-d0.genome_equivalents * vaf * d0.p_intact_fragment)]); }
    $('#det', host).innerHTML = lineChart(dw, 200, [{pts, c: css('--accent'), area: true, w: 2.2}], {x0: -5, x1: -1, y0: 0, y1: 1, xticks: [-5,-4,-3,-2,-1], xfmt: v => ({'-5':'0,001 %','-4':'0,01 %','-3':'0,1 %','-2':'1 %','-1':'10 %'})[v], yticks: [0,.5,1], ylabel: 'P(≥1 kopya)', xlabel: 'VAF',
      points: det.map(x => ({x: Math.log10(x.vaf), y: x.p_detect, c: css('--accent2'), tip: `VAF ${x.vaf} → P ${nf(x.p_detect,3)} (λ ${nf(x.lambda,2)})`, t: nf(x.p_detect,2)}))}) + `<div class="small muted">Amplikon ${d0.amplicon_bp} bp · bütün fragman olasılığı <b class="mono">${nf(d0.p_intact_fragment,3)}</b> · ${nf(d0.genome_equivalents,0)} genom eşdeğeri · noktalar motor çıktısı, eğri aynı formül.</div>`; }
}

/* ---- kanıt (91 kural) ---- */
function tabEvidence(host, pm, E) {
  const f = STATE.evFilter || 'all', q = (STATE.evQ || '').toLowerCase();
  const cls = e => e.sev ? e.sev : e.status === 'profile' ? 'NA' : e.status === 'pend' ? 'PEND' : 'NA';
  const all = CAT.map(r => ({r, e: E[r.id] || {status:'ne', recs:[]}}));
  const cnt = {eval: all.filter(x => x.e.sev).length, FAIL: all.filter(x => x.e.sev === 'FAIL').length, WARN: all.filter(x => x.e.sev === 'WARN').length, INFO: all.filter(x => x.e.sev === 'INFO' || x.e.sev === 'PASS').length,
    ne: all.filter(x => !x.e.sev && (x.e.status === 'ne' || x.e.status === 'na' || x.e.status === 'panel' || x.e.status === 'scan')).length, profile: all.filter(x => x.e.status === 'profile').length, pend: all.filter(x => x.e.status === 'pend').length};
  const pass = x => (f === 'all' || (f === 'eval' && x.e.sev) || x.e.sev === f || (f === 'ne' && !x.e.sev && x.e.status !== 'profile' && x.e.status !== 'pend') || (f === 'profile' && x.e.status === 'profile') || (f === 'pend' && x.e.status === 'pend')) && (!q || (x.r.id + ' ' + x.r.name + ' ' + x.r.method).toLowerCase().includes(q));
  host.innerHTML = `<div class="stack"><div class="card"><div class="bd">
    <div class="row" style="justify-content:space-between;gap:16px"><div><div class="eyebrow">Bu çift için QC kataloğu</div><h2 style="font-size:19px;margin-top:4px">${cnt.eval} kural kanıtla değerlendirildi · ${cnt.FAIL} FAIL · ${cnt.WARN} WARN</h2>
      <div class="small muted" style="margin-top:4px">${cnt.ne} kural bu çift için koşulmadı ya da panel düzeyinde · ${cnt.profile} kural ${esc(pm.profile)} profilinde uygulanmaz · ${cnt.pend} kural altyapı bekliyor. Ölçümü, eşiği ve modeli görmek için bir kurala tıklayın.</div></div>
      <div style="display:flex;height:12px;border-radius:6px;overflow:hidden;width:min(320px,100%);background:var(--sunk)">${[['FAIL','--fail'],['WARN','--warn'],['INFO','--info']].map(([k, v]) => `<i style="width:${cnt[k]/CAT.length*100}%;background:var(${v})"></i>`).join('')}<i style="width:${cnt.ne/CAT.length*100}%;background:var(--line2)"></i><i style="width:${cnt.pend/CAT.length*100}%;background:var(--probe);opacity:.6"></i></div></div>
    <div class="row" style="margin-top:14px;gap:10px"><div class="seg" id="evSeg">${[['all','Tümü ' + CAT.length],['FAIL','FAIL ' + cnt.FAIL],['WARN','WARN ' + cnt.WARN],['eval','Değerlendirilen ' + cnt.eval],['ne','Koşulmadı ' + cnt.ne],['profile','Profil dışı ' + cnt.profile],['pend','Bekliyor ' + cnt.pend]].map(([k, t]) => `<button data-f="${k}" aria-pressed="${f === k}">${t}</button>`).join('')}</div>
      <input type="text" id="evQ" placeholder="Kural ara: TH-02, dimer, varyant…" value="${esc(STATE.evQ || '')}" style="max-width:260px"></div></div></div>
    <div class="card"><div class="bd" id="evBody">${FAMORDER.map(fam => { const rows = all.filter(x => x.r.fam === fam && pass(x)); if (!rows.length) return '';
      const fs = all.filter(x => x.r.fam === fam); const m = k => fs.filter(x => cls(x.e) === k || (k === 'INFO' && x.e.sev === 'PASS')).length;
      return `<div class="famhd"><h3>${FAMTR[fam]}</h3><span class="mono tiny faint">${fs.length} kural</span><div class="meter">${[['FAIL','--fail'],['WARN','--warn'],['INFO','--info'],['NA','--line2'],['PEND','--probe']].map(([k, v]) => `<i style="width:${m(k)/fs.length*100}%;background:var(${v})"></i>`).join('')}</div></div>
        <div class="rulegrid">${rows.map(({r, e}) => `<button class="rule s-${cls(e)}" data-rule="${r.id}"><div class="top"><span class="id">${r.id}</span>${e.sev ? sev(e.sev) : sev(cls(e), e.status === 'profile' ? 'profil dışı' : e.status === 'pend' ? 'bekliyor' : e.status === 'scan' ? 'tarama' : e.status === 'panel' ? 'panel' : 'koşulmadı')}</div><div class="nm">${esc(r.name)}</div><div class="vv">${e.sev ? esc(roleSummary(e)) : esc(e.reason || '')}</div></button>`).join('')}</div>`; }).join('') || '<div class="empty">Filtreye uyan kural yok.</div>'}</div></div></div>`;
  $$('#evSeg button', host).forEach(b => b.onclick = () => { STATE.evFilter = b.dataset.f; rerender(); });
  const qi = $('#evQ', host); qi.oninput = () => { STATE.evQ = qi.value; clearTimeout(qi._t); qi._t = setTimeout(() => { rerender(); const n = $('#evQ'); if (n) { n.focus(); n.setSelectionRange(n.value.length, n.value.length); } }, 250); };
}
function roleSummary(e) { const by = {}; e.recs.filter(r => !r.ne && r.sev).forEach(r => by[r.role || '—'] = r.sev); return Object.entries(by).map(([k, v]) => `${k}:${v}`).join(' · '); }

/* ---- kural çekmecesi ---- */
document.addEventListener('click', e => {
  const b = e.target.closest('[data-rule]'); if (!b) return; const id = b.dataset.rule;
  const pm = STATE.view === 'cift' ? getPair(STATE.pair) : STATE.view === 'dokumantasyon' && STATE.docSection === 'cift' ? getPair(STATE.docPair) : null; const E = pm ? (pm._E || evidenceFor(pm)) : null; openRule(id, E ? E[id] : null, b);
});
function openRule(id, e, opener) {
  const r = CATI[id];
  if (!r) return;
  const st = STATUS[r.status] || {t:'Bekleyen', c:'NA'};
  const dv = document.createElement('dialog');
  dv.className = 'qc-dialog';
  dv.setAttribute('aria-labelledby', 'qc-rule-title');
  dv.innerHTML = `<div class="qc-dialog-body">
    <header class="qc-dialog-head">
      <div><div class="eyebrow">${esc(FAMTR[r.fam])}</div><span class="mono qc-code">${esc(id)}</span>
      <h2 id="qc-rule-title">${esc(r.name)}</h2></div>
      <button class="iconbtn" type="button" data-close aria-label="Kapat" autofocus>${ic('x')}</button>
    </header>
    <div class="row qc-dialog-meta"><span class="small muted">Katalog kapsamı</span>${sev(st.c, st.t)}</div>
    <div class="qc-explain">
      <section><div class="qc-explain-mark" aria-hidden="true">01</div><div><h3>Neyi kontrol eder?</h3><p>${esc(r.check)}</p></div></section>
      <section><div class="qc-explain-mark" aria-hidden="true">02</div><div><h3>Neden önemli?</h3><p>${esc(r.why)}</p></div></section>
    </div>
    ${QC_SCOPE_NOTES[id] ? `<div class="note">${ic('info')}<div><b>Bu kurulumun kapsamı</b><p>${esc(QC_SCOPE_NOTES[id])}</p></div></div>` : ''}
    ${e && e.sev ? `<div class="row qc-dialog-meta"><span class="small muted">Seçili çiftin sonucu</span>${sev(e.sev)}</div>` : ''}
    ${e ? (e.recs.length ? e.recs.map(x => `<div style="border:1px solid var(--line);border-radius:10px;padding:12px;margin-bottom:10px">
        <div class="row" style="gap:6px;margin-bottom:8px">${x.role ? `<span class="tag">${esc(x.role)}</span>` : ''}<span class="tag ${x.src === 'tarayıcı' ? 'br' : 'eng'}">${esc(x.src)}</span>${x.id !== id ? `<span class="tag">${esc(x.id)}</span>` : ''}<span class="sp" style="flex:1"></span>${x.sev ? sev(x.sev) : sev('NA', 'değerlendirilmedi')}</div>
        ${x.value && x.value.alignment ? dimerHTML(x.value.alignment) : ''}
        ${x.value != null ? kvHTML(flatten(x.value)) : ''}
        ${x.thr && Object.keys(x.thr).length ? `<div class="small" style="margin-top:8px"><span class="faint">Eşik</span> ${Object.entries(x.thr).map(([k, v]) => `<span class="chip mono">${esc(k)}: ${esc(v)}</span>`).join(' ')}</div>` : ''}
        ${x.note ? `<div class="small muted" style="margin-top:8px">${esc(x.note)}</div>` : ''}
        ${x.model ? `<details style="margin-top:8px"><summary class="tiny" style="cursor:pointer">model künyesi</summary>${kvHTML(x.model)}</details>` : ''}</div>`).join('')
        : `<div class="note">${ic('info')}<div>${esc(e.reason || 'Bu çift için kayıt yok.')}</div></div>`) : ''}
    ${e && e.recs.length ? `<details style="margin-top:6px"><summary class="small" style="cursor:pointer;font-weight:600">Ham kanıt kaydı (JSON)</summary><pre class="json" style="margin-top:8px">${esc(evidenceJson(e.recs, 2))}</pre></details>` : ''}
    </div>`;
  document.body.appendChild(dv);
  dv.addEventListener('click', ev => { if (ev.target === dv || ev.target.closest('[data-close]')) dv.close(); });
  dv.addEventListener('close', () => {
    dv.remove();
    if (opener?.isConnected) opener.focus({preventScroll:true});
  }, {once:true});
  dv.showModal();
}
function flatten(v, pre='', out={}) { if (v == null || typeof v !== 'object' || Array.isArray(v)) { out[pre || 'değer'] = Array.isArray(v) ? (v.length > 12 ? `[${v.slice(0, 12).map(fmtv).join(', ')} … ${v.length}]` : JSON.stringify(v.map(fmtv))) : fmtv(v); return out; }
  for (const [k, x] of Object.entries(v)) { if (k === 'alignment' || technicalField(k)) continue; if (x && typeof x === 'object' && !Array.isArray(x) && Object.keys(x).length < 12) flatten(x, pre ? pre + '.' + k : k, out); else out[pre ? pre + '.' + k : k] = Array.isArray(x) ? (x.length > 10 ? `[${x.slice(0, 10).map(y => typeof y === 'object' ? JSON.stringify(y) : fmtv(y)).join(', ')} … ${x.length}]` : JSON.stringify(x)) : fmtv(x); }
  return out; }
