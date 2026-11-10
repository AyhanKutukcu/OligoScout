/* ================= GENEL BAKIŞ ================= */
RENDER.genel = el => {
  const nPass = PROJ.filter(p => p.verdict === 'PASS').length, nWarn = PROJ.filter(p => p.verdict === 'WARN').length, nFail = PROJ.filter(p => p.verdict === 'FAIL').length;
  const ct = D.meta.contract; const offt = PROJ.reduce((s, p) => s + (+p.raw.offtarget_total || 0), 0);
  const underPrimers = new Set(); PROJ.forEach(p => { varsIn(p.fStart, p.fEnd).forEach(v => underPrimers.add(v.rsid + v.pos)); varsIn(p.rStart, p.rEnd).forEach(v => underPrimers.add(v.rsid + v.pos)); });
  const fams = FAMORDER.map(f => ({f, rs: CAT.filter(r => r.fam === f)}));
  const ring = (() => { const R = 46, C = 2 * Math.PI * R; const full = CAT_FULL / CAT.length, part = (CAT_IMPL - CAT_FULL) / CAT.length;
    return `<svg width="132" height="132" viewBox="0 0 132 132" role="img" aria-label="${CAT_IMPL}/${CAT.length} kural uygulanmış"><circle cx="66" cy="66" r="${R}" fill="none" stroke="${css('--sunk')}" stroke-width="13"/>
      <circle cx="66" cy="66" r="${R}" fill="none" stroke="${css('--accent')}" stroke-width="13" stroke-dasharray="${C*full} ${C}" transform="rotate(-90 66 66)" stroke-linecap="butt"/>
      <circle cx="66" cy="66" r="${R}" fill="none" stroke="${css('--warn')}" stroke-width="13" stroke-dasharray="${C*part} ${C}" stroke-dashoffset="${-C*full}" transform="rotate(-90 66 66)"/>
      <text x="66" y="64" text-anchor="middle" style="font:700 26px var(--display);fill:${css('--ink')}">${CAT_IMPL}</text><text x="66" y="82" text-anchor="middle" class="ax">/ ${CAT.length} kural</text></svg>`; })();
  el.innerHTML = `
  <div class="stack">
    <div class="card"><div class="bd hero">
      <div style="min-width:0">
        <div class="eyebrow">Aktif proje · gerçek motor çıktıları</div>
        <h2 style="font-size:24px;margin:6px 0 8px">chr22 ctDNA paneli — 20 aday primer çifti</h2>
        <p class="muted" style="margin:0 0 12px;max-width:68ch">Hedef bölge <span class="mono">${esc(D.meta.header)}</span>. Çiftler <code>design_primers.py</code> ile üretildi, C++ çekirdeğiyle tüm GRCh38'de tarandı (${ni(ct.binding_records)} bağlanma kaydı, ${ct.shards} shard, ${nf(+ct.search_seconds,1)} s), gnomAD v4 varyantlarıyla ve multipleks çapraz-reaktiviteyle sıralandı.</p>
        <div class="row" style="gap:6px">
          <span class="chip">profil <b>ctdna</b> · 60–100 bp</span>
          <span class="chip">Na⁺ <b>50</b> · Mg²⁺ <b>1,5</b> · dNTP <b>0,6</b> mM · primer <b>200</b> nM · Ta <b>60</b> °C</span>
          <span class="chip mono">param ${esc(D.engine_params.registry_digest.slice(0,12))}…</span>
          <span class="chip">çapa <b>${ct.anchor_length}</b> nt · maks <b>${ct.max_mismatches}</b> uyumsuzluk · ürün ${ct.min_amplicon_length}–${ct.max_amplicon_length} bp</span>
        </div>
      </div>
      <div class="ringbox">${ring}<div class="small" style="line-height:1.7"><div><i class="dot" style="background:var(--accent)"></i> ${CAT_FULL} tam uygulandı</div><div><i class="dot" style="background:var(--warn)"></i> ${CAT_IMPL-CAT_FULL} kısmi (SPC-08)</div><div><i class="dot" style="background:var(--sunk);border:1px solid var(--line2)"></i> ${CAT.length-CAT_IMPL} altyapı/veri bekliyor</div><a href="#katalog" class="small">Kataloğu aç →</a></div></div>
    </div></div>

    <div class="g4">
      ${kpi('Aday çift · son karar', PROJ.length, '', `${nPass} PASS · ${nWarn} WARN · ${nFail} FAIL (varyant katmanı dahil)`, [[nPass,'--ok'],[nWarn,'--warn'],[nFail,'--fail']])}
      ${kpi('Genom-geneli bağlanma', ni(ct.binding_records), '', `${ni(ct.unique_primer_searches)} arama · 3′ ${ct.anchor_length} nt çapa`)}
      ${kpi('Hedef dışı ürün', offt, '', `${PROJ.filter(p => +p.raw.offtarget_total > 0).length} çiftte · hiçbiri tam eşleşme değil`)}
      ${kpi('gnomAD varyantı', VARS.length, '', `${underPrimers.size} tanesi en az bir primerin altında`)}
    </div>

    <div class="card"><div class="hd"><h2>Şablon üzerindeki çiftler</h2><div class="seg" id="ovSeg"><button data-w="z" aria-pressed="${STATE.ovWin !== 'all'}">Çift bölgesi</button><button data-w="all" aria-pressed="${STATE.ovWin === 'all'}">Tüm şablon</button></div><span class="sp"></span><div class="legend"><span><i class="dot" style="background:var(--fwd)"></i>Forward</span><span><i class="dot" style="background:var(--rev)"></i>Reverse</span><span><i class="dot" style="background:var(--fail)"></i>AF ≥ %1</span><span><i class="dot" style="background:var(--warn)"></i>AF ≥ %0,1</span><span><i class="dot" style="background:var(--ink3)"></i>nadir</span></div></div>
      <div class="bd scrollx" id="ovTrack"></div></div>

    <div class="card"><div class="hd"><h2>Panel sıralaması</h2><span class="faint small">final_ranked.tsv · özgüllük ve varyant katmanı sonrası</span><span class="sp"></span><span class="tag cal">skor kalibre değil · REP-03</span></div>
      <div class="bd flush scrollx"><table class="t" style="min-width:880px"><thead><tr><th>#</th><th>Çift</th><th>Primerler (5′→3′)</th><th class="n">Ürün</th><th class="n">Tm F / R</th><th class="n">θ(Ta) F / R</th><th>Özgüllük</th><th>Varyant</th><th>Karar</th><th>Termo skoru</th></tr></thead>
      <tbody>${PROJ.map(p => `<tr class="click" data-pair="${p.key}"><td class="mono">${p.rank}</td><td class="mono"><b>${p.id}</b></td>
        <td class="seq small"><span style="color:var(--fwd)">▸</span> ${p.F}<br><span style="color:var(--rev)">◂</span> ${p.R}</td>
        <td class="n">${p.prod}</td><td class="n">${nf(p.tmF,2)}<br>${nf(p.tmR,2)}</td><td class="n">${nf(p.thetaF,3)}<br>${nf(p.thetaR,3)}</td>
        <td>${sev(p.spec)} <span class="tiny faint">${+p.raw.offtarget_total ? p.raw.offtarget_total + ' hedef dışı' : ''}</span></td><td>${sev(p.varv)}${p.raw.variant_driver ? `<div class="tiny faint mono" style="margin-top:3px;max-width:190px;white-space:normal">${esc(p.raw.variant_driver)}</div>` : ''}</td><td>${sev(p.verdict)}</td>
        <td><div style="display:flex;align-items:center;gap:8px"><div style="width:70px;height:6px;background:var(--sunk);border-radius:4px;overflow:hidden"><i style="display:block;height:100%;width:${Math.min(100,p.score/0.35*100)}%;background:var(--accent)"></i></div><span class="mono tiny">${nf(p.score,4)}</span></div></td></tr>`).join('')}</tbody></table></div></div>

    <div class="g2">
      <div class="card"><div class="hd"><h2>Boru hattı</h2><span class="sp"></span><span class="faint small">bu panelin geçtiği adımlar</span></div><div class="bd">${pipeline()}</div></div>
      <div class="card"><div class="hd"><h2>Aile başına kapsam</h2><span class="sp"></span><a class="small" href="#katalog">tüm kurallar →</a></div><div class="bd" style="display:grid;gap:7px">
        ${fams.map(({f, rs}) => { const full = rs.filter(r => r.status === 'etiketli_uygulandi').length, part = rs.filter(r => r.status === 'etiketli_kismi').length, pend = rs.length - full - part;
          return `<div style="display:grid;grid-template-columns:130px 1fr 44px;gap:10px;align-items:center;font-size:12.5px"><span>${FAMTR[f]}</span><div style="display:flex;height:9px;border-radius:5px;overflow:hidden;background:var(--sunk)">${'<i style="background:var(--accent);width:'+full/rs.length*100+'%"></i>'}${'<i style="background:var(--warn);width:'+part/rs.length*100+'%"></i>'}${'<i style="background:var(--probe);opacity:.5;width:'+pend/rs.length*100+'%"></i>'}</div><span class="mono tiny faint" style="text-align:right">${full+part}/${rs.length}</span></div>`; }).join('')}
      </div></div>
    </div>

    <div class="card"><div class="hd"><h2>Açık maddeler</h2><span class="faint small">acik_maddeler.tsv · kodla kapatılamayan, karar/altyapı bekleyen işler</span></div>
      <div class="bd flush scrollx"><table class="t" style="min-width:720px"><thead><tr><th>Madde</th><th>Alan</th><th>Tür</th><th>Bu planda</th><th>Gerekçe</th></tr></thead><tbody>
      ${D.open_items.map(o => `<tr><td class="mono"><b>${esc(o.madde)}</b></td><td class="mono small">${esc(o.alan)}</td><td><span class="tag">${esc(o.tur)}</span></td><td>${sev(o.bu_planda_kapatilabilir === 'evet' ? 'PASS' : o.bu_planda_kapatilabilir === 'kismen' ? 'WARN' : 'NA', o.bu_planda_kapatilabilir)}</td><td class="small muted">${esc(o.gerekce)}</td></tr>`).join('')}
      </tbody></table></div></div>
  </div>`;
  const lo = Math.min(...PROJ.map(p => p.fStart)) - 25, hi = Math.max(...PROJ.map(p => p.rEnd)) + 25;
  $('#ovTrack', el).innerHTML = overviewTrack(Math.max(720, W_OF($('#ovTrack', el)) - 36), {pairs: PROJ.map(p => ({...p, short: p.id.replace('DESIGN_', '#')})), sel: null, win: STATE.ovWin === 'all' ? null : [lo, hi]});
  $$('#ovSeg button', el).forEach(b => b.onclick = () => { STATE.ovWin = b.dataset.w; rerender(); });
};
function kpi(k, v, unit, d, bar) {
  const tot = bar ? bar.reduce((s, b) => s + b[0], 0) : 0;
  return `<div class="kpi"><div class="k">${k}</div><div class="v">${v}${unit ? `<small>${unit}</small>` : ''}</div><div class="d">${d}</div>${bar ? `<div class="bar">${bar.map(b => `<i style="width:${b[0]/tot*100}%;background:var(${b[1]})"></i>`).join('')}</div>` : ''}</div>`;
}
function pipeline() {
  const ct = D.meta.contract;
  const steps = [
    ['Tasarım', 'design_primers.py · ctdna', '2000 aday çift · sert filtreler + Pareto', 'PASS'],
    ['Dizi / çift QC', 'sequence_qc.py', 'SEQ · CMP · SS-02/04 · PAIR · AMP', 'PASS'],
    ['Genom taraması', 'oligoscout_panel (C++20, FM-index)', `${ni(ct.binding_records)} bağlanma · ${ct.product_records} ürün · ${nf(+ct.search_seconds,1)} s`, 'PASS'],
    ['Termo + özgüllük QC', 'panel_qc.py', 'TH · DIM · SS-01/03 · SPC-01…09 · COND', 'PASS'],
    ['Varyant katmanı', 'variant_qc.py · gnomAD v4', `${VARS.length} varyant · VAR-01…08`, 'PASS'],
    ['Multipleks', 'multiplex_qc.py', 'MPX-01…03 · DIM-06 · çapraz ürün', 'WARN'],
    ['Yeniden sıralama', 'rerank.py', 'final_ranked.tsv', 'PASS']];
  return `<ol style="list-style:none;margin:0;padding:0;display:grid;gap:0">${steps.map((s, i) => `<li style="display:grid;grid-template-columns:28px 1fr auto;gap:12px;align-items:start;padding:9px 0;${i < steps.length-1 ? 'border-bottom:1px dashed var(--line)' : ''}">
    <span style="width:24px;height:24px;border-radius:7px;background:var(--accent-soft);color:var(--accent2);display:grid;place-items:center;font:700 11px var(--display)">${i+1}</span>
    <div style="min-width:0"><div style="font-weight:600;font-size:13px">${s[0]} <span class="mono tiny faint">${s[1]}</span></div><div class="small muted">${s[2]}</div></div>${sev(s[3], s[3] === 'WARN' ? 'çapraz ürün' : 'tamam')}</li>`).join('')}</ol>`;
}

/* ================= PRIMER TASARLA ================= */
const DSG = { src: 'proje', fasta: '', fname: '', region: 'chr22:20000550-20000870' };
RENDER.tasarla = el => {
  const pr = PROFILES.ctdna;
  el.innerHTML = `
  <div style="display:grid;grid-template-columns:minmax(0,410px) minmax(0,1fr);gap:18px;align-items:start" id="dgrid">
    <form class="card" id="dform" autocomplete="off">
      <div class="fsec"><h3><span class="n">1</span>PCR şablonu</h3>
        <div class="seg" id="srcSeg" style="margin-bottom:12px">
          <button type="button" data-src="proje">Proje şablonu</button><button type="button" data-src="bolge">Genomik bölge</button><button type="button" data-src="fasta">FASTA / dizi</button><button type="button" data-src="dosya">Dosya</button></div>
        <div data-pane="proje"><div class="note">${ic('info')}<div><b>${esc(D.meta.header)}</b><br>${TPL.length} bp · ${VARS.length} gnomAD v4 varyantı · varyant farkındalığı açılabilir</div></div></div>
        <div data-pane="bolge"><label class="f"><span>Bölge <em>chrom:başlangıç-bitiş (1-tabanlı)</em></span><input type="text" id="gregion" value="${esc(DSG.region)}"></label>
          <p class="tiny faint" style="margin:6px 0 0">Bu sürümde dizisi pakette bulunan bölge: <span class="mono">chr22:20.000.001–20.001.000</span>. Başka bir bölge için FASTA yapıştırın; genom-geneli dizi çekimi arka uç (FastAPI) ister.</p></div>
        <div data-pane="fasta"><label class="f"><span>Dizi <em>FASTA ya da ham nükleotit</em></span><textarea id="fasta" spellcheck="false" placeholder=">hedef&#10;ACGT...">${esc(DSG.fasta)}</textarea></label>
          <div class="row" style="gap:6px;margin-top:8px"><button type="button" class="btn sm" id="exRandom">Örnek: hedef.fa (1 kb sentetik)</button><button type="button" class="btn sm" id="exChr">Örnek: chr22 amplikon bölgesi</button></div></div>
        <div data-pane="dosya"><label class="f"><span>FASTA dosyası <em>.fa .fasta .txt</em></span><input type="file" id="ffile" accept=".fa,.fasta,.fna,.txt,text/plain"></label><div class="tiny faint" id="fname" style="margin-top:6px">${esc(DSG.fname)}</div></div>
        <div class="fgrid c3" style="margin-top:12px"><label class="f"><span>Aralık başı</span><input type="number" id="tFrom" placeholder="1"></label><label class="f"><span>Aralık sonu</span><input type="number" id="tTo" placeholder="son"></label><label class="f"><span>Dönen çift</span><select id="nRet"><option>10</option><option selected>20</option><option>50</option></select></label></div>
        <div class="tiny faint" id="tplInfo" style="margin-top:8px"></div>
      </div>
      <div class="fsec"><h3><span class="n">2</span>Ürün ve primer konumu</h3>
        <div class="fgrid" style="grid-template-columns:1fr 1fr"><label class="f"><span>Ürün min <em>bp</em></span><input type="number" id="pmin" value="${pr.product_bp[0]}"></label><label class="f"><span>Ürün max <em>bp</em></span><input type="number" id="pmax" value="${pr.product_bp[1]}"></label></div>
        <div class="fgrid" style="grid-template-columns:repeat(4,minmax(0,1fr));margin-top:12px">
          <label class="f"><span>F başı</span><input type="number" id="fFrom"></label><label class="f"><span>F sonu</span><input type="number" id="fTo"></label>
          <label class="f"><span>R başı</span><input type="number" id="rFrom"></label><label class="f"><span>R sonu</span><input type="number" id="rTo"></label></div>
        <div class="fgrid" style="grid-template-columns:1fr 1fr;margin-top:12px"><label class="f"><span>Hedef bölge başı <em>ürün kapsamalı</em></span><input type="number" id="cFrom"></label><label class="f"><span>Hedef bölge sonu</span><input type="number" id="cTo"></label></div>
        <p class="tiny faint" style="margin:8px 0 0">Konumlar şablon içinde 1-tabanlıdır; proje şablonunda <span class="mono">1 = chr22:20.000.001</span>. Boş alan kısıt yok demektir.</p>
      </div>
      <div class="fsec"><h3><span class="n">3</span>Primer parametreleri</h3>
        <label class="f" style="margin-bottom:12px"><span>Uygulama profili</span><select id="prof">${(D.app_profiles||[]).map(p => `<option value="${p.profile}" ${p.profile==='ctdna'?'selected':''}>${p.profile} — ${esc(p.description)}</option>`).join('')}</select></label>
        <div class="fgrid" style="grid-template-columns:repeat(3,minmax(0,1fr))">
          <label class="f"><span>Tm min</span><input type="number" step="0.5" id="tmin" value="${pr.tm_window_c[0]}"></label><label class="f"><span>Tm opt</span><input type="number" step="0.5" id="topt" value="60"></label><label class="f"><span>Tm max</span><input type="number" step="0.5" id="tmax" value="${pr.tm_window_c[1]}"></label>
          <label class="f"><span>Uzunluk min</span><input type="number" id="lmin" value="${pr.length_nt[0]}"></label><label class="f"><span>Uzunluk max</span><input type="number" id="lmax" value="${pr.length_nt[1]}"></label><label class="f"><span>maks ΔTm</span><input type="number" step="0.5" id="dtm" value="2"></label>
          <label class="f"><span>GC % min</span><input type="number" id="gmin" value="${pr.gc_fraction[0]*100}"></label><label class="f"><span>GC % max</span><input type="number" id="gmax" value="${pr.gc_fraction[1]*100}"></label><label class="f"><span>maks homopolimer</span><input type="number" id="homo" value="4"></label></div>
        <details style="margin-top:12px"><summary class="small" style="cursor:pointer;font-weight:600">Dimer ve clamp eşikleri (design.py)</summary>
          <div class="fgrid" style="grid-template-columns:repeat(3,minmax(0,1fr));margin-top:10px">
            <label class="f"><span>self ΔG</span><input type="number" step="0.5" id="sd" value="-6"></label><label class="f"><span>self 3′ ΔG</span><input type="number" step="0.5" id="sd3" value="-5"></label><label class="f"><span>θ(Ta) min</span><input type="number" step="0.05" id="thmin" value="0"></label>
            <label class="f"><span>hetero ΔG</span><input type="number" step="0.5" id="hd" value="-6"></label><label class="f"><span>hetero 3′ ΔG</span><input type="number" step="0.5" id="hd3" value="-5"></label><label class="f"><span>clamp min/max</span><input type="text" id="clamp" value="1-3"></label></div></details>
      </div>
      <div class="fsec" id="varSec"><h3><span class="n">4</span>Varyant farkındalığı <span class="tag">SNPCheck</span></h3>
        <label class="switch"><input type="checkbox" id="vOn" checked> gnomAD varyantı taşıyan adayları ele</label>
        <div class="fgrid" style="grid-template-columns:1fr 1fr;margin-top:12px"><label class="f"><span>AF eşiği</span><select id="vAf"><option value="0.01">≥ %1 (yaygın)</option><option value="0.001" selected>≥ %0,1</option><option value="0.0001">≥ %0,01</option></select></label>
          <label class="f"><span>Kapsam</span><select id="vScope"><option value="3p">yalnız 3′ son 5 nt</option><option value="all" selected>tüm bağlanma yeri</option></select></label></div>
        <p class="tiny faint" style="margin:8px 0 0" id="varNote">Elenen aday VAR-02 (3′ penceresi) ya da VAR-03 (gövde) koduyla kaydedilir.</p></div>
      <div class="fsec"><h3><span class="n">5</span>Reaksiyon koşulları <span class="tag">COND</span></h3>
        <div class="fgrid" style="grid-template-columns:repeat(3,minmax(0,1fr))"><label class="f"><span>Na⁺/K⁺ mM</span><input type="number" step="1" id="na" value="50"></label><label class="f"><span>Mg²⁺ mM</span><input type="number" step="0.1" id="mg" value="1.5"></label><label class="f"><span>dNTP mM</span><input type="number" step="0.1" id="dntp" value="0.6"></label>
          <label class="f"><span>Primer nM</span><input type="number" id="pnm" value="200"></label><label class="f"><span>Ta °C</span><input type="number" step="0.5" id="ta" value="60"></label><label class="f"><span>Tuz modeli</span><select id="salt"><option value="owczarzy_2008">Owczarzy 2008</option><option value="owczarzy_2004">Owczarzy 2004</option><option value="von_ahsen_2001">von Ahsen 2001</option></select></label></div>
        <div class="note" style="margin-top:12px" id="freeMg"></div>
      </div>
      <div class="fsec" style="display:flex;gap:10px;align-items:center;position:sticky;bottom:0;background:var(--panel);border-radius:0 0 var(--r) var(--r)"><button class="btn pri" type="submit" id="goBtn" style="flex:1;justify-content:center;padding:11px">${ic('design')} Primer tasarla</button><button class="btn" type="button" id="resetBtn">Sıfırla</button></div>
    </form>
    <div id="dres" style="min-width:0"></div>
  </div>`;
  if (innerWidth < 1100) $('#dgrid', el).style.gridTemplateColumns = 'minmax(0,1fr)';
  const setSrc = s => { DSG.src = s; $$('#srcSeg button', el).forEach(b => b.setAttribute('aria-pressed', b.dataset.src === s)); $$('[data-pane]', el).forEach(p => p.hidden = p.dataset.pane !== s);
    const proj = s === 'proje' || s === 'bolge'; $('#vOn', el).disabled = !proj; $('#varNote', el).textContent = proj ? 'Elenen aday VAR-02 (3′ penceresi) ya da VAR-03 (gövde) koduyla kaydedilir.' : 'Varyant katmanı yalnız proje şablonunda (gnomAD önbelleği) kullanılabilir.'; tplInfo(); };
  $('#srcSeg', el).addEventListener('click', e => { const b = e.target.closest('button'); if (b) setSrc(b.dataset.src); });
  const tplInfo = () => { const t = currentTemplate(true); $('#tplInfo', el).innerHTML = t.err ? `<span style="color:var(--fail)">${esc(t.err)}</span>` : `<b class="mono">${t.seq.length.toLocaleString('tr-TR')} bp</b> · GC %${nf(OS.gcf(t.seq)*100,1)} · ${(t.seq.match(/[^ACGT]/g)||[]).length} belirsiz baz${t.range ? ` · aralık <span class="mono">${t.range[0]+1}–${t.range[1]}</span>` : ''}`; };
  ['fasta','gregion','tFrom','tTo'].forEach(id => $('#'+id, el).addEventListener('input', () => { if (id === 'fasta') DSG.fasta = $('#fasta', el).value; if (id === 'gregion') DSG.region = $('#gregion', el).value; tplInfo(); }));
  $('#ffile', el).addEventListener('change', e => { const f = e.target.files[0]; if (!f) return; const rd = new FileReader(); rd.onload = () => { DSG.fasta = String(rd.result); DSG.fname = f.name + ' · ' + f.size.toLocaleString('tr-TR') + ' bayt'; $('#fname', el).textContent = DSG.fname; $('#fasta', el).value = DSG.fasta; tplInfo(); }; rd.readAsText(f); });
  $('#exRandom', el).onclick = () => { $('#fasta', el).value = DSG.fasta = EX_RANDOM; tplInfo(); };
  $('#exChr', el).onclick = () => { $('#fasta', el).value = DSG.fasta = '>chr22:20000501-20000900 GRCh38.p14\n' + TPL.slice(500, 900).match(/.{1,60}/g).join('\n'); tplInfo(); };
  $('#prof', el).addEventListener('change', () => { const p = PROFILES[$('#prof', el).value]; if (!p) return;
    $('#pmin', el).value = p.product_bp[0]; $('#pmax', el).value = p.product_bp[1]; $('#tmin', el).value = p.tm_window_c[0]; $('#tmax', el).value = p.tm_window_c[1]; $('#topt', el).value = (p.tm_window_c[0] + p.tm_window_c[1]) / 2;
    $('#lmin', el).value = p.length_nt[0]; $('#lmax', el).value = p.length_nt[1]; $('#gmin', el).value = p.gc_fraction[0]*100; $('#gmax', el).value = p.gc_fraction[1]*100; });
  const fm = () => { const mg = +$('#mg', el).value, dn = +$('#dntp', el).value, f = mg - dn; const n = $('#freeMg', el);
    n.className = 'note' + (f <= 0 ? ' warn' : ''); n.innerHTML = `${ic(f <= 0 ? 'warn' : 'info')}<div>Serbest Mg²⁺ <b class="mono">${nf(f,2)} mM</b> = [Mg] − [dNTP] <span class="faint">(COND-02)</span>${f <= 0 ? ' — dNTP tüm Mg²⁺’u bağlıyor; tuz düzeltmesi tek değerlikliye düşer.' : ''}</div>`; };
  ['mg','dntp'].forEach(id => $('#'+id, el).addEventListener('input', fm)); fm();
  $('#resetBtn', el).onclick = () => rerender();
  $('#dform', el).addEventListener('submit', e => { e.preventDefault(); runDesign(el); });
  setSrc(DSG.src);
  if (STATE.design) showDesign(el, STATE.design); else $('#dres', el).innerHTML = designIntro();
};
const EX_RANDOM = `>Random_DNA_1000bp (hedef.fa)
GGCCCCATTACCAATTGCTCCTCGTACTTGTGTGGTGCTTTGCAGCCGCTCACCCGAGATGTGACACACCTCGGCCGGGAAGAGTTCGCGGCCCCCTATT
GGACCGTGGCAGTGGTAATGCAGCGGTTAAATCGGTTACTGCCGAGACCCTAATGATGTAGAATTTTAAGAAACGACTTTAACACGATTAAGGAGCGCTC
ATGCGGCCAAATAGAAATGTGTCGCAAGTTGTCTATTGCCGCAGCTTGGGCCCAACCCAAAAAAGACCCTCTAAGCCGAATAGTCCACAGTACGTGCGTA
AACTTCTTCCGCTTCCGGTGCTCCCGGCTAGGTGATACGTAGTAAGGTGTACATGATAGCCTCGTAGATACTAGCTGAAATATAAGACAGATGGACCGGA
AGGTTGAATGAGGTAAAGGCTCCCCTAAAGGTTGGAGTTAATGAGGAGACTTGTTTTCGTCCTTCATAGAATGCCGTAACCTTCCTCTTGCCATCCCCTT
CAGATCACGATTGTCCCCGATTTAGACATAATCCACCCTTCTTCTCAACTGTCGGCCCGTTTTAATCGGACGGCTTTATGCCGAAGAAAGGCGGTTACGC
TTTCGCTGAATTTTCAACAACCTGCACCAAAGACGTTCTACACGCCACTCAGGGATCGGGTTGTCCTGTCCCCTAAATTCATATCCTTTCCAGCTGATCC
GGCCAGAGTTGTTCTGAATGGCCATGTACTGAGTACAACATCTGCCTCGCGTGGCTTCCCGCGTGACGAAATGCACCCTTACGGTAATAGAGCCGCTCAC
GCCTTCGGACCATGTCAGGGGTAGATACTGTGGATGTCGAACACACGTCTGGACTGCAGTATTGAAGTGCTGAAACCTTAATTAAGTAAAGATTGTGGGC
GTTGCACCTACAGTCCAGTGGAGAGCCCGTGGCTCCCGCGGAACCCAACTGTTTTCGGGGGCTACGGCAATCTCTAACTCCTGACACCTACCTGTGTACA`;
function parseFasta(txt) { const lines = String(txt).split(/\r?\n/); let name = ''; const seq = [];
  for (const l of lines) { if (l.startsWith('>')) { if (!name) name = l.slice(1).trim(); else break; } else seq.push(l); }
  return {name, seq: seq.join('').replace(/[\s\d]/g, '').toUpperCase().replace(/U/g, 'T')}; }
function currentTemplate(quiet) {
  const v = id => { const e = document.getElementById(id); return e && e.value !== '' ? +e.value : null; };
  let seq, name, off = 0, hasVar = false, range = null;
  if (DSG.src === 'proje') { seq = TPL; name = D.meta.header; off = OFF; hasVar = true; }
  else if (DSG.src === 'bolge') {
    const m = String(DSG.region).replace(/[,.\s]/g, '').match(/^(chr)?(\w+):(\d+)-(\d+)$/i);
    if (!m) return {err:'Bölge biçimi: chr22:20000550-20000870'};
    const ch = 'chr' + m[2].replace(/^chr/i, ''), a = +m[3], b = +m[4];
    if (ch.toLowerCase() !== 'chr22' || a < OFF + 1 || b > OFF + TPL.length || b <= a) return {err:`Bu bölgenin dizisi pakette yok (yalnız chr22:${(OFF+1).toLocaleString('tr-TR')}–${(OFF+TPL.length).toLocaleString('tr-TR')}). FASTA yapıştırın.`};
    seq = TPL; name = D.meta.header; off = OFF; hasVar = true; range = [a - 1 - OFF, b - OFF];
  } else { const p = parseFasta(DSG.fasta); seq = p.seq; name = p.name; if (!seq) return {err:'Dizi yapıştırın ya da dosya seçin.'};
    const bad = seq.match(/[^ACGTN]/g); if (bad) return {err:`Tanınmayan karakter: ${[...new Set(bad)].slice(0,5).join(', ')} (SEQ-01)`};
    if (seq.length > 60000) return {err:'Tarayıcı tasarımı için en fazla 60 kb; daha uzun hedefler design_primers.py ile.'}; }
  const f = v('tFrom'), t = v('tTo'); if (f || t) { const a = Math.max(0, (f || 1) - 1), b = Math.min(seq.length, t || seq.length); range = range ? [Math.max(range[0], a), Math.min(range[1], b)] : [a, b]; }
  return {seq, name, off, hasVar, range};
}
function designIntro() {
  return `<div class="card"><div class="bd" style="padding:28px">
    <div class="eyebrow">Nasıl çalışır</div><h2 style="font-size:20px;margin:6px 0 10px">Referans motorun kendisi, tarayıcıda</h2>
    <p class="muted" style="max-width:62ch;margin:0 0 16px">Bu form <code>design_primers.py</code> akışını birebir çalıştırır: aday üretimi → sert filtreler (CMP-01…04, TH-01, TH-03, DIM-01/02) → çift kurgulama (PAIR-01, DIM-03/04) → Pareto katmanları ve kalibre edilmemiş skaler skor. Tm, θ(Ta) ve dimer hesapları projedeki parametre tablolarıyla yapılır; aynı şablonda Python motoruyla 2000/2000 çift aynı sırada çıkar (Doğrulama sayfasında canlı sınanır).</p>
    <div class="g3">${[['1','Şablonu seçin','Proje bölgesi, genomik koordinat ya da yapıştırılan FASTA'],['2','Ürün boyunu verin','Primer-BLAST gibi F/R konum aralığı ve kapsanacak hedef bölge'],['3','Sonuçları inceleyin','Her çift harita, varyant, dimer, termodinamik ve 91 kuralla açılır']].map(s => `<div style="border:1px solid var(--line);border-radius:10px;padding:12px"><div class="mono tiny faint">${s[0]}</div><div style="font-weight:600;margin:2px 0">${s[1]}</div><div class="small muted">${s[2]}</div></div>`).join('')}</div>
    <div class="note" style="margin-top:16px">${ic('info')}<div>Genom-geneli özgüllük (SPC) tarayıcıda yapılamaz; sonuç sayfası, seçtiğiniz çiftleri C++ taramasına göndermek için hazır panel TSV'sini ve komutu verir.</div></div></div></div>`;
}
function readCons(el) {
  const v = id => +$('#'+id, el).value; const cl = ($('#clamp', el).value || '1-3').split(/[-–,]/).map(Number);
  return {length_min:v('lmin'), length_max:v('lmax'), tm_min:v('tmin'), tm_max:v('tmax'), tm_opt:v('topt'), gc_min:v('gmin')/100, gc_max:v('gmax')/100,
    max_homopolymer:v('homo'), max_repeat_units:4, gc_clamp_min:cl[0] ?? 1, gc_clamp_max:cl[1] ?? 3, self_dimer_dg:v('sd'), self_dimer3p_dg:v('sd3'),
    hetero_dimer_dg:v('hd'), hetero_dimer3p_dg:v('hd3'), product_min:v('pmin'), product_max:v('pmax'), max_delta_tm:v('dtm'), max_pairs:2000, theta_min:v('thmin'), salt_model:$('#salt', el).value};
}
function runDesign(el) {
  const t = currentTemplate(); const res = $('#dres', el);
  if (t.err) { res.innerHTML = `<div class="note warn">${ic('warn')}<div>${esc(t.err)}</div></div>`; return; }
  const cons = readCons(el), cond = {na:+$('#na', el).value, mg:+$('#mg', el).value, dntp:+$('#dntp', el).value, pnm:+$('#pnm', el).value, ta:+$('#ta', el).value};
  if (cons.product_min > cons.product_max || cons.length_min > cons.length_max || cons.tm_min > cons.tm_max) { res.innerHTML = `<div class="note warn">${ic('warn')}<div>Min değerler max değerlerden büyük olamaz.</div></div>`; return; }
  const pos = id => { const e = $('#'+id, el); return e.value === '' ? null : +e.value - 1; };
  const r0 = t.range || [0, t.seq.length];
  const clip = (a, b) => [Math.max(r0[0], a ?? r0[0]), Math.min(r0[1], b != null ? b + 1 : r0[1])];
  const rel = x => x == null ? null : x + (DSG.src === 'bolge' ? 0 : 0);
  const opts = {leftRegion: clip(rel(pos('fFrom')), rel(pos('fTo'))), rightRegion: clip(rel(pos('rFrom')), rel(pos('rTo')))};
  const cf = pos('cFrom'), ctt = pos('cTo'); if (cf != null && ctt != null) opts.mustCover = [cf, ctt + 1];
  const vOn = $('#vOn', el).checked && t.hasVar, vAf = +$('#vAf', el).value, vSc = $('#vScope', el).value;
  if (vOn) opts.variantPenalty = (st, en, o) => { const vs = varsIn(st, en - 1).filter(x => (x.af || 0) >= vAf); if (!vs.length) return null;
    const near = vs.filter(x => (o === 'F' ? (en - 1) - x.idx : x.idx - st) < 5);
    if (near.length) { const w = near.reduce((a, b) => (b.af > a.af ? b : a)); return {rule:'VAR-02', why:"3' penceresinde varyant", value:`${w.rsid} AF=${w.af.toExponential(1)}`}; }
    if (vSc === 'all') { const w = vs.reduce((a, b) => (b.af > a.af ? b : a)); return {rule:'VAR-03', why:'bağlanma yerinde varyant', value:`${w.rsid} AF=${w.af.toExponential(1)}`}; }
    return null; };
  const btn = $('#goBtn', el); btn.disabled = true; btn.innerHTML = '<span class="spin"></span> Hesaplanıyor…';
  res.innerHTML = `<div class="card"><div class="bd empty"><span class="spin" style="border-color:var(--line);border-top-color:var(--accent);display:inline-block;width:22px;height:22px"></span><div style="margin-top:10px">Adaylar üretiliyor, sert filtreler ve dimer taraması uygulanıyor…</div></div></div>`;
  setTimeout(() => {
    const t0 = performance.now(); let out;
    try { out = OS.design(t.seq, cons, cond, opts); } catch (err) { res.innerHTML = `<div class="note warn">${ic('warn')}<div>Hesaplama hatası: ${esc(err.message)}</div></div>`; btn.disabled = false; btn.innerHTML = ic('design') + ' Primer tasarla'; return; }
    const ms = performance.now() - t0; const nRet = +$('#nRet', el).value; const prof = $('#prof', el).value;
    const pms = out.pairs.slice(0, nRet).map((p, i) => ({
      key: 'L:' + Date.now().toString(36) + ':' + i, kind:'local', id: 'TASARIM-' + String(i + 1).padStart(2, '0'), label: 'Tasarım #' + (i + 1), short: '#' + (i + 1), rank: i + 1,
      F: p.left.sequence, R: p.right.sequence, fStart: p.left.start, fEnd: p.left.end - 1, rStart: p.right.start, rEnd: p.right.end - 1, prod: p.product_length,
      tmF: p.left.tm_c, tmR: p.right.tm_c, dtm: p.delta_tm, thetaF: p.left.theta_ta, thetaR: p.right.theta_ta, score: p.score, pareto: p.pareto_rank, d: p,
      tpl: t.seq, tplOff: t.off, tplName: t.name, chrom: t.off ? CHR : (t.name || 'şablon').split(/\s/)[0], hasVar: t.hasVar, profile: prof, cond, salt: cons.salt_model,
      match: PROJ.find(q => q.F === p.left.sequence && q.R === p.right.sequence) }));
    STATE.local = pms; STATE.design = {out, pms, ms, cons, cond, t, nRet, prof, vOn};
    btn.disabled = false; btn.innerHTML = ic('design') + ' Primer tasarla'; showDesign(el, STATE.design);
  }, 30);
}
function showDesign(el, S) {
  const res = $('#dres', el); const {out, pms, ms, cons, t} = S;
  const rc = {}; out.rejections.forEach(r => rc[r.rule] = (rc[r.rule] || 0) + 1); const rcs = Object.entries(rc).sort((a, b) => b[1] - a[1]); const rmax = rcs.length ? rcs[0][1] : 1;
  if (!out.pairs.length) { res.innerHTML = `<div class="card"><div class="bd"><div class="verdict fail"><div class="ic">${ic('x')}</div><div><h3>Uygun çift bulunamadı</h3><div class="small">${out.left.length} forward, ${out.right.length} reverse aday kaldı. Ürün aralığını genişletin ya da Tm/GC pencerelerini gevşetin.</div></div></div>${rejBars(rcs, rmax)}</div></div>`; return; }
  res.innerHTML = `<div class="stack">
    <div class="g4">${kpi('Bulunan çift', ni(out.pairs.length), '', `ilk ${pms.length} gösteriliyor · ${ni(ms)} ms`)}${kpi('Forward aday', ni(out.left.length), '', 'sert filtreler + dimer sonrası')}${kpi('Reverse aday', ni(out.right.length), '', '')}${kpi('Elenen', ni(out.rejections.length), '', 'her elenme kural koduyla')}</div>
    <div class="card"><div class="hd"><h2>Grafik özet</h2><span class="faint small">${esc(t.name || 'şablon')} · ${t.seq.length.toLocaleString('tr-TR')} bp</span><span class="sp"></span><span class="tag cal">Pareto + kalibre edilmemiş skor</span></div><div class="bd scrollx" id="dTrack"></div></div>
    <div class="card"><div class="hd"><h2>Elenme gerekçeleri</h2><span class="faint small">design_rejections.tsv ile aynı kodlar</span></div><div class="bd">${rejBars(rcs, rmax)}</div></div>
    ${pms.map(blastBlock).join('')}
    <div class="card"><div class="hd"><h2>Genom-geneli özgüllük için</h2><span class="sp"></span>${copyBtn(panelTSV(pms), 'Panel TSV kopyala')}</div><div class="bd small muted">Bu çiftler henüz GRCh38'de taranmadı (SPC-01…09). Panel TSV'sini kopyalayıp projede bir dosyaya kaydedin ve <code>RUN_PANEL.cmd panel.tsv 166</code> çalıştırın; çıkan <span class="mono">panel_qc_records.jsonl</span> <code>web/build_web_data.py</code> ile bu siteye aktarılır.<pre class="json" style="margin-top:10px">${esc(panelTSV(pms).split('\n').slice(0,6).join('\n'))}${pms.length > 5 ? '\n…' : ''}</pre></div></div>
  </div>`;
  const lo2 = Math.min(...pms.map(p => p.fStart)), hi2 = Math.max(...pms.map(p => p.rEnd)); const pad2 = Math.max(30, (hi2 - lo2) * .15);
  $('#dTrack', res).innerHTML = overviewTrack(Math.max(640, W_OF($('#dTrack', res)) - 36), {tpl: t.seq, off: t.off, pairs: pms, showVars: t.hasVar, region: t.range, win: (hi2 - lo2) < t.seq.length * .7 ? [lo2 - pad2, hi2 + pad2] : null});
}
function rejBars(rcs, rmax) { return `<div style="display:grid;gap:6px">${rcs.map(([k, n]) => `<div style="display:grid;grid-template-columns:76px minmax(0,1fr) 64px;gap:10px;align-items:center;font-size:12.5px"><span class="mono"><b>${k}</b></span><div style="position:relative;height:18px;background:var(--sunk);border-radius:5px;overflow:hidden"><i style="position:absolute;inset:0 auto 0 0;width:${n/rmax*100}%;background:color-mix(in srgb,var(--accent) 55%,transparent)"></i><span class="tiny" style="position:relative;padding-left:8px;line-height:18px;color:var(--ink)">${esc((CATI[k] && CATI[k].name) || (k === 'LIMIT' ? 'dimer taraması sınırı' : k))}</span></div><span class="mono tiny" style="text-align:right">${ni(n)}</span></div>`).join('')}</div>`; }
function panelTSV(pms) { return 'pair_id\tprimer1\tprimer2\n' + pms.map(p => `${p.id}\t${p.F}\t${p.R}`).join('\n'); }
function blastBlock(p) {
  const row = (lab, s, strand, st, en, tm, g, self, self3) => `<tr><td><b>${lab}</b></td><td class="seq">${seqC(s)} ${copyBtn(s)}</td><td>${strand}</td><td class="n">${s.length}</td><td class="n">${st}</td><td class="n">${en}</td><td class="n">${nf(tm,2)}</td><td class="n">${nf(g*100,2)}</td><td class="n">${nf(self,2)}</td><td class="n">${nf(self3,2)}</td></tr>`;
  const d = p.d;
  return `<div class="card"><div class="hd"><h2>Primer çifti ${p.rank}</h2><span class="chip">Pareto katmanı <b>${p.pareto}</b></span><span class="chip">skor <b class="mono">${nf(p.score,4)}</b></span>${p.match ? `<span class="chip" style="border-color:var(--accent);color:var(--accent2)">${ic('target')} proje panelinde: <b>${p.match.id}</b> · genom taraması var</span>` : ''}<span class="sp"></span><button class="btn sm" data-open="${p.key}">${ic('arrow')} Çift incelemesinde aç</button></div>
    <div class="bd flush scrollx"><table class="t" style="min-width:860px"><thead><tr><th></th><th>Dizi (5′→3′)</th><th>İplik</th><th class="n">Uzunluk</th><th class="n">Başlangıç</th><th class="n">Bitiş</th><th class="n">Tm °C</th><th class="n">GC %</th><th class="n">Self ΔG</th><th class="n">Self 3′ ΔG</th></tr></thead><tbody>
      ${row('Forward', p.F, 'Plus', p.fStart + 1, p.fEnd + 1, p.tmF, d.left.gc, d.left.self_dimer_dg, d.left.self_dimer3p_dg)}
      ${row('Reverse', p.R, 'Minus', p.rEnd + 1, p.rStart + 1, p.tmR, d.right.gc, d.right.self_dimer_dg, d.right.self_dimer3p_dg)}
      <tr><td colspan="10" class="small"><span class="muted">Ürün uzunluğu</span> <b class="mono">${p.prod} bp</b> · <span class="muted">ΔTm</span> <b class="mono">${nf(p.dtm,2)} °C</b> · <span class="muted">hetero ΔG</span> <b class="mono">${nf(d.hetero_dg,2)}</b> / 3′ <b class="mono">${nf(d.hetero3p_dg,2)}</b> kcal/mol · <span class="muted">θ(Ta)</span> <b class="mono">${nf(p.thetaF,3)} / ${nf(p.thetaR,3)}</b>${p.tplOff ? ` · <span class="mono faint">${CHR}:${gposAny(p,p.fStart)}–${gposAny(p,p.rEnd)}</span>` : ''}</td></tr></tbody></table></div></div>`;
}
document.addEventListener('click', e => {
  const o = e.target.closest('[data-open]'); if (o) { STATE.pair = o.dataset.open; STATE.ptab = 'harita'; go('cift'); return; }
  const r = e.target.closest('[data-pair]'); if (r && !e.target.closest('.copy')) { STATE.pair = r.dataset.pair; if (STATE.view !== 'cift') { STATE.ptab = 'harita'; go('cift'); } else rerender(); }
});
