/* ================= GENEL BAKIŞ ================= */
RENDER.genel = el => {
  const ready = LIVE.status;
  const capabilities = [
    ['Referans genom', ready?.reference.ready, 'GRCh38.p14 · 1–22, X ve Y'],
    ['Genom özgüllüğü', ready?.cpp && ready?.indexes.primary.ready, 'GRCh38 üzerinde 24 kromozom'],
    ['Bisülfit', ready?.cpp && ready?.indexes.CT.ready && ready?.indexes.GA.ready, 'CT ve GA indeksleri · BIS-04'],
    ['Anotasyonlar', ready?.annotations.isoform && ready?.annotations.pseudogene && ready?.annotations.segdup,
      'GENCODE v47 eksonları, psödogen ve segmental duplikasyon'],
    ['gnomAD v4', ready?.variants.ready, ready?.variants.api ? 'GraphQL API · genom ve ekzom bölgeleri' : ready?.variants.regions.length
      ? ready.variants.regions.map(r => r.chrom + ':' + ni(r.start) + '–' + ni(r.end)).join(' · ')
      : 'Yerel varyant önbelleğinin bildirdiği bölgeler'],
  ];
  el.innerHTML = `<div class="stack">
    <section class="card"><div class="bd intro-hero overview-hero">
      <div class="overview-copy">
        <h2>Primer tasarımından genom kontrolüne</h2>
        <p>İnsan referans genomu üzerinde PCR primerleri tasarlayın. Özgüllük, erime sıcaklığı, dimer riski ve gnomAD varyantlarını birlikte inceleyin.</p>
        <div class="row live-actions"><a class="btn pri" href="/design">Primer tasarla</a>
          <a class="btn" href="/multiplex">Multipleks tasarla</a><a class="btn" href="/docs">Dokümantasyon</a></div>
      </div>
      <div class="overview-visual" aria-hidden="true"><img src="__OVERVIEW_ART__" alt="" width="1254" height="1254"></div>
    </div></section>
    <div class="g3">
      ${[['design','Tasarım','GRCh38 bölgesine ve PCR koşullarına göre aday çiftler bulun.','tasarla'],
        ['target','Özgüllük','Primerlerin 24 kromozomda oluşturabileceği hedef dışı ürünleri inceleyin.','yerel'],
        ['list','Kanıt','Termodinamik, varyant ve dimer kayıtlarını ölçüm ve eşikleriyle okuyun.','cift']].map(([icon,title,desc,view]) =>
        `<a class="card intro-step" href="${routeUrl(view)}"><span class="overview-step-icon">${ic(icon)}</span><h3>${title}</h3><p>${desc}</p></a>`).join('')}
    </div>
    <section class="card"><div class="hd"><h2>Analiz ve veri kaynakları</h2></div><div class="bd">
      ${!ready ? liveNote(localAvailability()) : ''}
      <div class="live-capabilities">${capabilities.map(([label,ok,detail]) => `<div>
        <b>${label}</b><span class="small muted">${esc(detail)}</span>${ready ? sev(ok ? 'INFO' : 'NA', ok ? 'Hazır' : 'Eksik') : ''}</div>`).join('')}</div>
      <p class="tiny faint" style="margin:10px 0 0">gnomAD hedef bölgeyi sorgular ve başarılı yanıtı önbelleğe alır. Veri alınamayan bölgeler değerlendirilmemiş olarak gösterilir.</p>
    </div></section>
  </div>`;
};

RENDER.panel = el => {
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
        <div class="eyebrow">Kayıtlı çalışma · arşivlenmiş motor çıktıları</div>
        <h2 style="font-size:24px;margin:6px 0 8px">chr22 ctDNA paneli — 20 aday primer çifti</h2>
        <p class="muted" style="margin:0 0 12px;max-width:68ch">Hedef bölge <span class="mono">${esc(D.meta.header)}</span>. Çiftler tasarlandı ve GRCh38 genelinde tarandı (${ni(ct.binding_records)} bağlanma kaydı, ${ct.shards} shard, ${nf(+ct.search_seconds,1)} s), gnomAD v4 varyantlarıyla ve multipleks çapraz-reaktiviteyle sıralandı.</p>
        <div class="row" style="gap:6px">
          <span class="chip">profil <b>ctdna</b> · 60–100 bp</span>
          <span class="chip">Na⁺ <b>50</b> · Mg²⁺ <b>1,5</b> · dNTP <b>0,6</b> mM · primer <b>200</b> nM · Ta <b>60</b> °C</span>
          <span class="chip">çapa <b>${ct.anchor_length}</b> nt · maks <b>${ct.max_mismatches}</b> uyumsuzluk · ürün ${ct.min_amplicon_length}–${ct.max_amplicon_length} bp</span>
        </div>
      </div>
      <div class="ringbox">${ring}<div class="small" style="line-height:1.7"><div><i class="dot" style="background:var(--accent)"></i> ${CAT_FULL} tam uygulandı</div><div><i class="dot" style="background:var(--warn)"></i> ${CAT_IMPL-CAT_FULL} kısmi (SPC-08)</div><div><i class="dot" style="background:var(--sunk);border:1px solid var(--line2)"></i> ${CAT.length-CAT_IMPL} altyapı/veri bekliyor</div><a href="/catalog" class="small">Kataloğu aç →</a></div></div>
    </div></div>

    <div class="g4">
      ${kpi('Aday çift · son karar', PROJ.length, '', `${nPass} PASS · ${nWarn} WARN · ${nFail} FAIL (varyant katmanı dahil)`, [[nPass,'--ok'],[nWarn,'--warn'],[nFail,'--fail']])}
      ${kpi('Genom-geneli bağlanma', ni(ct.binding_records), '', `${ni(ct.unique_primer_searches)} arama · 3′ ${ct.anchor_length} nt çapa`)}
      ${kpi('Hedef dışı ürün', offt, '', `${PROJ.filter(p => +p.raw.offtarget_total > 0).length} çiftte · ${PROJ.reduce((n,p) => n + (+p.raw.offtarget_exact || 0), 0)} tam eşleşmeli`)}
      ${kpi('gnomAD varyantı', VARS.length, '', `${underPrimers.size} tanesi en az bir primerin altında`)}
    </div>

    <div class="card"><div class="hd"><h2>Şablon üzerindeki çiftler</h2><div class="seg" id="ovSeg"><button data-w="z" aria-pressed="${STATE.ovWin !== 'all'}">Çift bölgesi</button><button data-w="all" aria-pressed="${STATE.ovWin === 'all'}">Tüm şablon</button></div><span class="sp"></span><div class="legend"><span><i class="dot" style="background:var(--fwd)"></i>Forward</span><span><i class="dot" style="background:var(--rev)"></i>Reverse</span><span><i class="dot" style="background:var(--fail)"></i>AF ≥ %1</span><span><i class="dot" style="background:var(--warn)"></i>AF ≥ %0,1</span><span><i class="dot" style="background:var(--ink3)"></i>nadir</span></div></div>
      <div class="bd scrollx" id="ovTrack"></div></div>

    <div class="card"><div class="hd"><h2>Panel sıralaması</h2><span class="faint small">final_ranked.tsv · özgüllük ve varyant katmanı sonrası</span></div>
      <div class="bd flush scrollx"><table class="t" style="min-width:880px"><thead><tr><th>#</th><th>Çift</th><th>Primerler (5′→3′)</th><th class="n">Ürün</th><th class="n">Tm F / R</th><th class="n">θ(Ta) F / R</th><th>Özgüllük</th><th>Varyant</th><th>Karar</th><th>Termo skoru</th></tr></thead>
      <tbody>${PROJ.map(p => `<tr class="click" data-pair="${p.key}"><td class="mono">${p.rank}</td><td class="mono"><b>${p.id}</b></td>
        <td class="seq small"><span style="color:var(--fwd)">▸</span> ${p.F}<br><span style="color:var(--rev)">◂</span> ${p.R}</td>
        <td class="n">${p.prod}</td><td class="n">${nf(p.tmF,2)}<br>${nf(p.tmR,2)}</td><td class="n">${nf(p.thetaF,3)}<br>${nf(p.thetaR,3)}</td>
        <td>${sev(p.spec)} <span class="tiny faint">${+p.raw.offtarget_total ? p.raw.offtarget_total + ' hedef dışı' : ''}</span></td><td>${sev(p.varv)}${p.raw.variant_driver ? `<div class="tiny faint mono" style="margin-top:3px;max-width:190px;white-space:normal">${esc(p.raw.variant_driver)}</div>` : ''}</td><td>${sev(p.verdict)}</td>
        <td><div style="display:flex;align-items:center;gap:8px"><div style="width:70px;height:6px;background:var(--sunk);border-radius:4px;overflow:hidden"><i style="display:block;height:100%;width:${Math.min(100,p.score/0.35*100)}%;background:var(--accent)"></i></div><span class="mono tiny">${nf(p.score,4)}</span></div></td></tr>`).join('')}</tbody></table></div></div>

    <div class="g2">
      <div class="card"><div class="hd"><h2>Boru hattı</h2><span class="sp"></span><span class="faint small">bu panelin geçtiği adımlar</span></div><div class="bd">${pipeline()}</div></div>
      <div class="card"><div class="hd"><h2>Aile başına kapsam</h2><span class="sp"></span><a class="small" href="/catalog">tüm kurallar →</a></div><div class="bd" style="display:grid;gap:7px">
        ${fams.map(({f, rs}) => { const full = rs.filter(r => r.status === 'etiketli_uygulandi').length, part = rs.filter(r => r.status === 'etiketli_kismi').length, pend = rs.length - full - part;
          return `<div style="display:grid;grid-template-columns:130px 1fr 44px;gap:10px;align-items:center;font-size:12.5px"><span>${FAMTR[f]}</span><div style="display:flex;height:9px;border-radius:5px;overflow:hidden;background:var(--sunk)">${'<i style="background:var(--accent);width:'+full/rs.length*100+'%"></i>'}${'<i style="background:var(--warn);width:'+part/rs.length*100+'%"></i>'}${'<i style="background:var(--probe);opacity:.5;width:'+pend/rs.length*100+'%"></i>'}</div><span class="mono tiny faint" style="text-align:right">${full+part}/${rs.length}</span></div>`; }).join('')}
      </div></div>
    </div>

    <div class="card"><div class="hd"><h2>Açık maddeler</h2><span class="faint small">acik_maddeler.tsv · kodla kapatılamayan, karar/altyapı bekleyen işler</span></div>
      <div class="bd flush scrollx"><table class="t" style="min-width:720px"><thead><tr><th>Madde</th><th>Alan</th><th>Tür</th><th>Bu planda</th><th>Gerekçe</th></tr></thead><tbody>
      ${D.open_items.filter(o => !o.madde.startsWith('REP-')).map(o => `<tr><td class="mono"><b>${esc(o.madde)}</b></td><td class="mono small">${esc(o.alan)}</td><td><span class="tag">${esc(o.tur)}</span></td><td>${sev(o.bu_planda_kapatilabilir === 'evet' ? 'PASS' : o.bu_planda_kapatilabilir === 'kismen' ? 'WARN' : 'NA', o.bu_planda_kapatilabilir)}</td><td class="small muted">${esc(o.gerekce)}</td></tr>`).join('')}
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
    ['Genom taraması', 'GRCh38 özgüllük taraması', `${ni(ct.binding_records)} bağlanma · ${ct.product_records} ürün · ${nf(+ct.search_seconds,1)} s`, 'PASS'],
    ['Termo + özgüllük QC', 'panel_qc.py', 'TH · DIM · SS-01/03 · SPC-01…09 · COND', 'PASS'],
    ['Varyant katmanı', 'variant_qc.py · gnomAD v4', `${VARS.length} varyant · VAR-01…08`, 'PASS'],
    ['Multipleks', 'multiplex_qc.py', 'MPX-01…03 · DIM-06 · çapraz ürün', 'WARN'],
    ['Yeniden sıralama', 'rerank.py', 'final_ranked.tsv', 'PASS']];
  return `<ol style="list-style:none;margin:0;padding:0;display:grid;gap:0">${steps.map((s, i) => `<li style="display:grid;grid-template-columns:28px 1fr auto;gap:12px;align-items:start;padding:9px 0;${i < steps.length-1 ? 'border-bottom:1px dashed var(--line)' : ''}">
    <span style="width:24px;height:24px;border-radius:7px;background:var(--accent-soft);color:var(--accent2);display:grid;place-items:center;font:700 11px var(--display)">${i+1}</span>
    <div style="min-width:0"><div style="font-weight:600;font-size:13px">${s[0]} <span class="mono tiny faint">${s[1]}</span></div><div class="small muted">${s[2]}</div></div>${sev(s[3], s[3] === 'WARN' ? 'çapraz ürün' : 'tamam')}</li>`).join('')}</ol>`;
}

/* ================= PRIMER TASARLA ================= */
const DSG = { src: 'bolge', fasta: '', fname: '', region: '', profile:'sanger' };
RENDER.tasarla = el => {
  const pr = PROFILES[DSG.profile];
  el.innerHTML = `
  <div style="display:grid;grid-template-columns:minmax(0,410px) minmax(0,1fr);gap:18px;align-items:start" id="dgrid">
    <form class="card" id="dform" autocomplete="off">
      <div class="fsec"><h3><span class="n">1</span>PCR şablonu</h3>
        <div class="seg" id="srcSeg" style="margin-bottom:12px">
          <button type="button" data-src="bolge">Genomik bölge</button><button type="button" data-src="fasta">FASTA / dizi</button><button type="button" data-src="dosya">Dosya</button></div>
        <div data-pane="bolge"><label class="f"><span>Bölge <em>chrom:başlangıç-bitiş (1-tabanlı)</em></span><input type="text" id="gregion" value="${esc(DSG.region)}"></label>
          <p class="tiny faint" style="margin:6px 0 0">Yerel modda GRCh38.p14 üzerindeki 1–22, X veya Y kromozomundan en fazla 5.000 baz seçin.</p></div>
        <div data-pane="fasta"><label class="f"><span>Dizi <em>Canlı modda GRCh38 koordinatlı FASTA</em></span><textarea id="fasta" spellcheck="false" placeholder="Koordinat başlıklı FASTA girin">${esc(DSG.fasta)}</textarea></label>
          </div>
        <div data-pane="dosya"><label class="f"><span>FASTA dosyası <em>.fa .fasta .txt</em></span><input type="file" id="ffile" accept=".fa,.fasta,.fna,.txt,text/plain"></label><div class="tiny faint" id="fname" style="margin-top:6px">${esc(DSG.fname)}</div></div>
        <div class="fgrid c3" style="margin-top:12px"><label class="f"><span>Aralık başı</span><input type="number" id="tFrom" placeholder="1"></label><label class="f"><span>Aralık sonu</span><input type="number" id="tTo" placeholder="son"></label><label class="f"><span>En fazla çift</span><select id="nRet"><option>10</option><option selected>20</option><option>50</option></select></label></div>
        <div class="tiny faint" id="tplInfo" style="margin-top:8px"></div>
      </div>
      <div class="fsec"><h3><span class="n">2</span>Ürün ve primer konumu</h3>
        <div class="fgrid" style="grid-template-columns:1fr 1fr"><label class="f"><span>Ürün min <em>bp</em></span><input type="number" id="pmin" value="${pr.product_bp[0]}"></label><label class="f"><span>Ürün max <em>bp</em></span><input type="number" id="pmax" value="${pr.product_bp[1]}"></label></div>
        <div class="fgrid" style="grid-template-columns:repeat(4,minmax(0,1fr));margin-top:12px">
          <label class="f"><span>F başı</span><input type="number" id="fFrom"></label><label class="f"><span>F sonu</span><input type="number" id="fTo"></label>
          <label class="f"><span>R başı</span><input type="number" id="rFrom"></label><label class="f"><span>R sonu</span><input type="number" id="rTo"></label></div>
        <div class="fgrid" style="grid-template-columns:1fr 1fr;margin-top:12px"><label class="f"><span>Hedef bölge başı <em>ürün kapsamalı</em></span><input type="number" id="cFrom"></label><label class="f"><span>Hedef bölge sonu</span><input type="number" id="cTo"></label></div>
        <p class="tiny faint" style="margin:8px 0 0">Konumlar seçilen şablon içinde 1 tabanlıdır. Genomik bölge seçildiğinde 1, o bölgenin ilk bazıdır. Boş alan kısıt koymaz.</p>
      </div>
      <div class="fsec"><h3><span class="n">3</span>Primer parametreleri</h3>
        <label class="f" style="margin-bottom:12px"><span>Uygulama profili</span><select id="prof">${(D.app_profiles||[]).map(p => `<option value="${p.profile}" ${p.profile===DSG.profile?'selected':''}>${p.profile} — ${esc(p.description)}</option>`).join('')}</select></label>
        <label class="f" id="convField" hidden><span>Bisülfit dönüşümü</span><select id="conversion"><option value="CT">CT</option><option value="GA">GA</option></select></label>
        <div class="fgrid" style="grid-template-columns:repeat(3,minmax(0,1fr))">
          <label class="f"><span>Tm min</span><input type="number" step="0.5" id="tmin" value="${pr.tm_window_c[0]}"></label><label class="f"><span>Tm opt</span><input type="number" step="0.5" id="topt" value="${(pr.tm_window_c[0] + pr.tm_window_c[1]) / 2}"></label><label class="f"><span>Tm max</span><input type="number" step="0.5" id="tmax" value="${pr.tm_window_c[1]}"></label>
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
        <p class="tiny faint" style="margin:8px 0 0" id="varNote">Elenen aday VAR-02 (3′ penceresi) ya da VAR-03 (gövde) koduyla kaydedilir.</p>
        <div style="margin-top:12px"><label class="f"><span>Genom özgüllük taraması</span><select id="sMode"><option value="anchored" selected>Hızlı · tam 12-nt 3′ çapa</option><option value="sensitive">Hassas · Primer-BLAST kuralı (son 5 nt içinde en fazla 1 uyumsuzluk)</option></select></label></div></div>
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
    $('#vOn', el).disabled = LIVE.mode === 'file' && s !== 'bolge';
    $('#varNote', el).textContent = 'Filtre, seçilen bölge için gnomAD API veya önbellek verisi bulunduğunda uygulanır.';
    tplInfo();
  };
  $('#srcSeg', el).addEventListener('click', e => { const b = e.target.closest('button'); if (b) setSrc(b.dataset.src); });
  const tplInfo = () => {
    if (LIVE.mode !== 'file') {
      const info = $('#tplInfo', el);
      try {
        const region = liveDesignTemplate();
        info.textContent = region.region + ' · ' + ni(region.size) + ' bp · GRCh38.p14';
      } catch (error) { info.textContent = error.message; }
      return;
    }
    const t = currentTemplate(true); $('#tplInfo', el).innerHTML = t.err ? `<span style="color:var(--fail)">${esc(t.err)}</span>` : `<b class="mono">${t.seq.length.toLocaleString(numLocale())} bp</b> · GC ${pct(OS.gcf(t.seq)*100,1)} · ${(t.seq.match(/[^ACGT]/g)||[]).length} belirsiz baz${t.range ? ` · aralık <span class="mono">${t.range[0]+1}–${t.range[1]}</span>` : ''}`; };
  ['fasta','gregion','tFrom','tTo'].forEach(id => $('#'+id, el).addEventListener('input', () => { if (id === 'fasta') DSG.fasta = $('#fasta', el).value; if (id === 'gregion') DSG.region = $('#gregion', el).value; tplInfo(); }));
  $('#ffile', el).addEventListener('change', e => { const f = e.target.files[0]; if (!f) return; if (f.size > 100000) { $('#fname', el).textContent = 'FASTA dosyası 100 KB sınırını aşıyor.'; return; } const rd = new FileReader(); rd.onload = () => { DSG.fasta = String(rd.result); DSG.fname = f.name + ' · ' + f.size.toLocaleString(numLocale()) + ' bayt'; $('#fname', el).textContent = DSG.fname; $('#fasta', el).value = DSG.fasta; tplInfo(); }; rd.readAsText(f); });
  $('#prof', el).addEventListener('change', () => { const p = PROFILES[$('#prof', el).value]; if (!p) return; DSG.profile = p.profile;
    $('#convField', el).hidden = p.profile !== 'bisulfite' || LIVE.mode === 'file';
    $('#pmin', el).value = p.product_bp[0]; $('#pmax', el).value = p.product_bp[1]; $('#tmin', el).value = p.tm_window_c[0]; $('#tmax', el).value = p.tm_window_c[1]; $('#topt', el).value = (p.tm_window_c[0] + p.tm_window_c[1]) / 2;
    $('#lmin', el).value = p.length_nt[0]; $('#lmax', el).value = p.length_nt[1]; $('#gmin', el).value = p.gc_fraction[0]*100; $('#gmax', el).value = p.gc_fraction[1]*100; });
  const fm = () => { const mg = +$('#mg', el).value, dn = +$('#dntp', el).value, f = mg - dn; const n = $('#freeMg', el);
    n.className = 'note' + (f <= 0 ? ' warn' : ''); n.innerHTML = `${ic(f <= 0 ? 'warn' : 'info')}<div>Serbest Mg²⁺ <b class="mono">${nf(Math.max(0,f),2)} mM</b> = max(0, [Mg] − [dNTP]) <span class="faint">(COND-02)</span>${f <= 0 ? ' — dNTP tüm Mg²⁺’u bağlıyor; tuz düzeltmesi tek değerlikliye düşer.' : ''}</div>`; };
  ['mg','dntp'].forEach(id => $('#'+id, el).addEventListener('input', fm)); fm();
  $('#resetBtn', el).onclick = () => {
    Object.assign(DSG, {src:'bolge', region:'', fasta:'', fname:'', profile:'sanger'});
    LIVE.errors.delete('design');
    rerender();
  };
  const form = $('#dform', el);
  $('#goBtn', el).disabled = LIVE.pending.has('design') || !liveReady();
  form.insertAdjacentHTML('beforeend', jobNotice('design') + (!liveReady() ? liveNote(localAvailability()) : ''));
  form.addEventListener('submit', e => {
    e.preventDefault();
    if (liveReady()) submitLiveDesign(el);
  });
  setSrc(DSG.src);
  if (STATE.design) showDesign(el, STATE.design);
  else if (!LIVE.pending.has('design')) $('#dres', el).innerHTML = designIntro();
};
function parseFasta(txt) { const lines = String(txt).split(/\r?\n/); let name = ''; const seq = [];
  for (const l of lines) { if (l.startsWith('>')) { if (!name) name = l.slice(1).trim(); else break; } else seq.push(l); }
  return {name, seq: seq.join('').replace(/[\s\d]/g, '').toUpperCase().replace(/U/g, 'T')}; }
function currentTemplate(quiet) {
  const v = id => { const e = document.getElementById(id); return e && e.value !== '' ? +e.value : null; };
  let seq, name, off = 0, hasVar = false, range = null;
  if (DSG.src === 'bolge') {
    const m = String(DSG.region).replace(/[,.\s]/g, '').match(/^(chr)?(\w+):(\d+)-(\d+)$/i);
    if (!m) return {err:'Bölgeyi kromozom:başlangıç-bitiş biçiminde girin.'};
    const ch = 'chr' + m[2].replace(/^chr/i, ''), a = +m[3], b = +m[4];
    if (ch.toLowerCase() !== 'chr22' || a < OFF + 1 || b > OFF + TPL.length || b <= a) return {err:`Bu bölgenin dizisi pakette yok (yalnız chr22:${(OFF+1).toLocaleString(numLocale())}–${(OFF+TPL.length).toLocaleString(numLocale())}). FASTA yapıştırın.`};
    seq = TPL; name = D.meta.header; off = OFF; hasVar = true; range = [a - 1 - OFF, b - OFF];
  } else { const p = parseFasta(DSG.fasta); seq = p.seq; name = p.name; if (!seq) return {err:'Dizi yapıştırın ya da dosya seçin.'};
    const bad = seq.match(/[^ACGTN]/g); if (bad) return {err:`Tanınmayan karakter: ${[...new Set(bad)].slice(0,5).join(', ')} (SEQ-01)`};
    if (seq.length > 60000) return {err:'Tarayıcı tasarımı için en fazla 60 kb; daha uzun hedefler design_primers.py ile.'}; }
  const f = v('tFrom'), t = v('tTo'); if (f || t) { const a = Math.max(0, (f || 1) - 1), b = Math.min(seq.length, t || seq.length); range = range ? [Math.max(range[0], a), Math.min(range[1], b)] : [a, b]; }
  return {seq, name, off, hasVar, range};
}
function designIntro() {
  return `<div class="card"><div class="bd" style="padding:28px">
    <div class="eyebrow">TASARIM AKIŞI</div><h2 style="margin:6px 0 10px">Bir hedef bölge seçin</h2>
    <p class="muted" style="max-width:68ch">Primer adayları sıralanır ve 24 kromozomda taranır. Sonuç listesi, beklenen tek ürün ve her primer için tek tam bağlanma yeri bulunan çiftleri gösterir.</p>
    <p class="small faint">FASTA dosyasının başlığı ve dizisi GRCh38.p14 aralığıyla doğrulanır.</p>
    <a class="btn sm" href="/docs">Giriş ve çıktı örneklerini gör</a>
    ${LIVE.mode === 'file' ? liveNote('Özgüllük taramasıyla primer tasarlamak için RUN_WEB.cmd dosyasını çalıştırın.') : ''}
    </div></div>`;
}
function readCons(el) {
  const v = id => +$('#'+id, el).value; const cl = ($('#clamp', el).value || '1-3').split(/[-–,]/).map(Number);
  return {length_min:v('lmin'), length_max:v('lmax'), tm_min:v('tmin'), tm_max:v('tmax'), tm_opt:v('topt'), gc_min:v('gmin')/100, gc_max:v('gmax')/100,
    max_homopolymer:v('homo'), max_repeat_units:4, gc_clamp_min:cl[0] ?? 1, gc_clamp_max:cl[1] ?? 3, self_dimer_dg:v('sd'), self_dimer3p_dg:v('sd3'),
    hetero_dimer_dg:v('hd'), hetero_dimer3p_dg:v('hd3'), product_min:v('pmin'), product_max:v('pmax'), max_delta_tm:v('dtm'), max_pairs:2000, theta_min:v('thmin'), salt_model:$('#salt', el).value};
}
function showDesign(el, S) {
  const res = $('#dres', el);
  const {pms, ms, t, live} = S;
  const result = live.result;
  const {screening, counts} = result;
  const rejections = Object.entries(result.rejections).sort((a, b) => b[1] - a[1]);
  const status = screening.scanned
    ? `${ni(screening.scanned)} aday tarandı · ${ni(screening.rejected)} aday elendi`
    : 'Ön filtreleri geçen aday yok.';
  const coverage = result.variant_coverage
    ? result.variant_source + (result.variant_filter ? ' · AF filtresi uygulandı' : ' · varyant kanıtları kaydedildi')
    : result.variant_reason || 'Dönüştürülmüş referans için varyant filtresi uygulanmadı.';
  const limits = (screening.unscanned || result.candidate_cap_reached)
    ? '<p class="tiny faint">Sıralama taranan adaylarla sınırlıdır. Aday listesi tasarım motorunun üst sınırına tabidir.</p>' : '';
  const rejectionPanel = rejections.length ? `<details class="card"><summary class="hd">Tasarımda elenme gerekçeleri</summary><div class="bd">${rejBars(rejections, rejections[0][1])}</div></details>` : '';
  if (!pms.length) {
    res.innerHTML = `<div class="stack">
      <div class="card"><div class="bd"><h2>Özgüllük filtresinden geçen çift bulunamadı</h2>
      <p class="small muted">${esc(status)} · ${esc(coverage)}</p>${limits}</div></div>${rejectionPanel}</div>`;
    return;
  }
  res.innerHTML = `<div class="stack">
    <div class="card"><div class="hd"><h2>Özgüllüğü doğrulanan çiftler</h2><span class="sp"></span>${copyBtn(panelTSV(pms), 'Panel TSV kopyala')}</div>
      <div class="bd"><b>${pms.length} çift</b><p class="small muted">${esc(status)} · ${ni(counts.pairs)} tasarım adayı · ${nf(ms / 1000, 1)} sn</p>
        <p class="small muted">${esc(coverage)}</p>${limits}
        <div class="row live-actions"><button class="btn sm" data-open="${esc(pms[0].key)}" data-open-tab="varyant">İlk çiftin grafiklerini aç</button>
        <span class="tiny muted">Varyant, Yapı & Dimer ve Termodinamik · Çift İnceleme</span></div>
        ${live.result.specificity.multiplex ? localToolLink('scan', 'Multipleks sonuçlarını aç', pms[0]) : ''}</div></div>
    <div class="card"><div class="hd"><h2>Grafik özet</h2><span class="faint small">${esc(t.name)}</span></div><div class="bd scrollx" id="dTrack"></div></div>
    ${rejectionPanel}
    ${pms.map(blastBlock).join('')}
  </div>`;
  const lo = Math.min(...pms.map(p => p.fStart)), hi = Math.max(...pms.map(p => p.rEnd));
  const pad = Math.max(30, (hi - lo) * .15);
  $('#dTrack', res).innerHTML = overviewTrack(Math.max(640, W_OF($('#dTrack', res)) - 36),
    {tpl:t.seq, off:t.off, pairs:pms, showVars:false, win:[lo - pad, hi + pad]});
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
      <tr><td colspan="10" class="small"><span class="muted">Ürün uzunluğu</span> <b class="mono">${p.prod} bp</b> · <span class="muted">ΔTm</span> <b class="mono">${nf(p.dtm,2)} °C</b> · <span class="muted">hetero ΔG</span> <b class="mono">${nf(d.hetero_dg,2)}</b> / 3′ <b class="mono">${nf(d.hetero3p_dg,2)}</b> kcal/mol · <span class="muted">θ(Ta)</span> <b class="mono">${nf(p.thetaF,3)} / ${nf(p.thetaR,3)}</b>${p.tplOff ? ` · <span class="mono faint">${esc(p.chrom)}:${gposAny(p,p.fStart)}–${gposAny(p,p.rEnd)}</span>` : ''}</td></tr></tbody></table></div></div>`;
}
document.addEventListener('click', e => {
  const o = e.target.closest('[data-open]'); if (o) { go('cift', {pair:o.dataset.open, ptab:o.dataset.openTab || 'harita'}); return; }
  const r = e.target.closest('[data-pair]');
  if (r && !e.target.closest('.copy')) {
    const pair = getPair(r.dataset.pair);
    if (pair?.kind === 'project') go('dokumantasyon/cift/' + pair.key, {ptab:'harita'});
    else if (pair) go('cift', {pair:pair.key, ptab:'harita'});
  }
});
