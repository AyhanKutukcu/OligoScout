
/* Multi-target primer design and compatible multiplex pool selection. */
const MPD = {
  targets:[{name:'', region:''}, {name:'', region:''}], profile:'sanger',
  candidates:3, tm_spread:3, product_gap:10, dimer_any:-6, dimer_3p:-5,
  product_min:PROFILES.sanger.product_bp[0], product_max:PROFILES.sanger.product_bp[1],
  variants:true, conditions:{na:50, mg:1.5, dntp:.6, pnm:200, ta:60},
  result:null, source:null, mode:'dG37_any', pms:[],
};

RENDER.multipleks = el => {
  const busy = LIVE.pending.has('multiplex_design');
  el.innerHTML = `<div class="stack">
    <section class="card"><div class="bd mx-intro">
      <div><div class="eyebrow">ÇOK HEDEFLİ PRIMER TASARIMI</div><h2>Birlikte çalışacak primerleri seçin</h2>
      <p class="muted">Hedef bölgelerinizi girin. Her hedef için özgüllüğü kontrol edilen adaylar üretilir; çapraz ürün, dimer, Tm ve ürün boyu koşullarına uyan çiftler aynı panelde birleştirilir.</p></div>
      <a class="btn sm" href="/docs/multiplex">Örnek sonuçları gör</a>
    </div></section>
    <div class="mx-layout">
      <form class="card" id="mxForm" autocomplete="off">
        <div class="fsec"><div class="row"><h3><span class="n">1</span>Hedef bölgeler</h3><span class="sp"></span><span class="small faint">2–8 hedef · GRCh38.p14</span></div>
          <div class="mx-targets">${MPD.targets.map((target,i) => `<div class="mx-target">
            <div class="mx-target-num">${i+1}</div>
            <label class="f"><span>Hedef adı</span><input type="text" name="target_name_${i}" value="${esc(target.name)}" maxlength="60" required aria-label="${i+1}. hedef adı"></label>
            <label class="f"><span>Genomik bölge</span><input type="text" name="target_region_${i}" value="${esc(target.region)}" required aria-label="${i+1}. hedef bölgesi"></label>
            <button class="iconbtn" type="button" data-mx-remove="${i}" aria-label="${i+1}. hedefi kaldır" ${MPD.targets.length<=2||busy?'disabled':''}>${ic('close')}</button>
          </div>`).join('')}</div>
          <button class="btn sm" type="button" id="mxAdd" ${MPD.targets.length>=8||busy?'disabled':''}>+ Hedef ekle</button>
          <p class="tiny faint">Koordinat biçimi: kromozom:başlangıç-bitiş. Her bölge 1 tabanlı ve en fazla 5.000 baz olmalı.</p>
          <details class="mx-batch"><summary class="small">Hedef listesini topluca yapıştır</summary>
            <label class="f"><span>Her satıra hedef adı ve genomik bölge</span><textarea id="mxBatch" aria-label="Toplu hedef listesi"></textarea></label>
            <p class="tiny faint">Ad ile koordinat arasına sekme veya boşluk koyun.</p><button class="btn sm" type="button" id="mxImport">Listeyi uygula</button>
          </details>
        </div>
        <div class="fsec"><h3><span class="n">2</span>Ortak tasarım koşulları</h3>
          <label class="f"><span>Uygulama profili</span><select name="profile" id="mxProfile">${(D.app_profiles||[]).filter(p=>p.profile!=='bisulfite').map(p=>`<option value="${p.profile}" ${MPD.profile===p.profile?'selected':''}>${esc(p.profile)} — ${esc(p.description)}</option>`).join('')}</select></label>
          <div class="fgrid c3" style="margin-top:12px">
            <label class="f"><span>Ürün min (bp)</span><input name="product_min" type="number" min="20" max="5000" value="${MPD.product_min}" required></label>
            <label class="f"><span>Ürün max (bp)</span><input name="product_max" type="number" min="20" max="5000" value="${MPD.product_max}" required></label>
            <label class="f"><span>Hedef başına aday</span><select name="candidates">${[1,2,3,4].map(n=>`<option value="${n}" ${n===+MPD.candidates?'selected':''}>${n}</option>`).join('')}</select></label>
          </div>
          <label class="switch" style="margin-top:14px"><input name="variants" type="checkbox" ${MPD.variants?'checked':''}>gnomAD verisi varsa varyantlı adayları ele</label>
          <p class="tiny faint">Her hedef ayrı genom taramasından geçer. Hedef ve aday sayısı arttıkça işlem süresi uzar.</p>
        </div>
        <div class="fsec"><h3><span class="n">3</span>Aynı havuza kabul koşulları</h3>
          <div class="fgrid" style="grid-template-columns:1fr 1fr">
            <label class="f"><span>En fazla Tm aralığı (°C)</span><input name="tm_spread" type="number" step=".1" min=".1" max="15" value="${MPD.tm_spread}" required></label>
            <label class="f"><span>En az ürün boyu farkı (bp)</span><input name="product_gap" type="number" min="0" max="500" value="${MPD.product_gap}" required></label>
            <label class="f"><span>Dimer alt sınırı (kcal/mol)</span><input name="dimer_any" type="number" step=".5" min="-30" max="0" value="${MPD.dimer_any}" required></label>
            <label class="f"><span>3′ dimer alt sınırı (kcal/mol)</span><input name="dimer_3p" type="number" step=".5" min="-30" max="0" value="${MPD.dimer_3p}" required></label>
          </div>
          <p class="tiny faint">Ürün boyu farkını 0 yaparsanız boy ayrımı koşulu uygulanmaz. Dimer enerjisi belirtilen sınırdan daha düşük olan çiftler aynı havuza alınmaz.</p>
          <details><summary class="small">Reaksiyon koşulları</summary><div style="margin-top:12px">${localConditionsHTML()}</div></details>
        </div>
        <div class="fsec">${jobNotice('multiplex_design')}${!liveReady()?liveNote(localAvailability()):''}
          <button class="btn pri" id="mxRun" type="submit" ${busy||!liveReady()?'disabled':''}>${ic('mpx')} Multipleks panel tasarla</button>
        </div>
      </form>
      <div id="mxResult" class="stack" aria-live="polite"></div>
    </div></div>`;
  const form=$('#mxForm',el);
  for(const [name,value] of Object.entries(MPD.conditions)) if(form.elements[name]) form.elements[name].value=value;
  const capture=()=>{
    MPD.targets=MPD.targets.map((target,i)=>({name:form.elements['target_name_'+i].value,region:form.elements['target_region_'+i].value}));
    for(const name of ['profile','product_min','product_max','candidates','tm_spread','product_gap','dimer_any','dimer_3p']) MPD[name]=form.elements[name].value;
    MPD.variants=form.elements.variants.checked;
    MPD.conditions=liveConditions(form);
  };
  form.addEventListener('input',capture);
  form.addEventListener('change',capture);
  $('#mxProfile',form).addEventListener('change',()=>{
    const p=PROFILES[MPD.profile]; MPD.product_min=p.product_bp[0]; MPD.product_max=p.product_bp[1]; rerender();
  });
  $('#mxAdd',form).onclick=()=>{capture(); if(MPD.targets.length<8)MPD.targets.push({name:'',region:''});rerender();};
  $$('[data-mx-remove]',form).forEach(b=>b.onclick=()=>{capture();MPD.targets.splice(+b.dataset.mxRemove,1);rerender();});
  $('#mxImport',form).onclick=()=>{
    try {
      const rows=$('#mxBatch',form).value.trim().split(/\r?\n/).filter(s=>s.trim());
      if(rows.length<2||rows.length>8)throw Error('2–8 hedef satırı girin.');
      MPD.targets=rows.map((row,i)=>{
        const m=row.trim().match(/^(.+?)\s+((?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y):[\d.,]+-[\d.,]+)$/i);
        if(!m)throw Error((i+1)+'. satırda hedef adı ve genom bölgesini kontrol edin.');
        return {name:m[1],region:readLiveRegion(m[2]).region};
      });
      LIVE.errors.delete('multiplex_design');rerender();
    } catch(error){LIVE.errors.set('multiplex_design',error.message);updateJobNotice('multiplex_design');}
  };
  form.onsubmit=e=>{
    e.preventDefault();capture();
    if(!form.reportValidity()||busy)return;
    try {
      const names=new Set(),regions=new Set();
      const targets=MPD.targets.map((t,i)=>{
        const name=t.name.trim(),region=readLiveRegion(t.region).region;
        if(names.has(name.toLocaleLowerCase())||regions.has(region))throw Error('Hedef adları ve bölgeleri birbirinden farklı olmalı.');
        names.add(name.toLocaleLowerCase());regions.add(region);
        return {name,region};
      });
      const payload={tool:'multiplex_design',targets,profile:MPD.profile,candidates:+MPD.candidates,
        constraints:{product_min:+MPD.product_min,product_max:+MPD.product_max},
        conditions:MPD.conditions,variant_filter:MPD.variants,variant_af:.001,variant_scope:'all',
        pool:Object.fromEntries(['tm_spread','product_gap','dimer_any','dimer_3p'].map(k=>[k,+MPD[k]]))};
      MPD.result=null;MPD.pms=[];$('#mxRun',form).disabled=true;
      renderMultiplexResult($('#mxResult',el));
      runLocalJob('multiplex_design',payload,(result,id)=>{
        acceptMultiplexResult(result,id);
        if(STATE.view==='multipleks'||STATE.view==='cift'||STATE.view==='yerel')rerender();
      }).finally(()=>{if(STATE.view==='multipleks')rerender();});
    } catch(error){LIVE.errors.set('multiplex_design',error.message);updateJobNotice('multiplex_design');}
  };
  renderMultiplexResult($('#mxResult',el));
};

function acceptMultiplexResult(result,id){
  MPD.result=result;MPD.source=id;
  const pms=result.targets.flatMap(target=>{
    const pairs=makeLivePairs(target.result,id);
    for(const p of pairs){
      p.targetName=target.name;p.targetId=target.id;
      p.label=target.name+' · Aday '+p.rank;p.short=target.id;
      p.sourceLabel='Multipleks · '+result.targets.length+' hedef';
    }
    return pairs;
  });
  const recent=[...new Set(STATE.local.map(p=>p.source).filter(Boolean))].slice(-2);
  STATE.local=[...STATE.local.filter(p=>recent.includes(p.source)),...pms];
  MPD.pms=pms;
  for(const target of result.targets){
    const scan=target.result.specificity;
    if(scan)applyScan(id,scan);
    for(const p of pms.filter(p=>p.targetId===target.id)){
      const qc=target.result.qc?.[p.id];if(qc)resultForPair(p,'qc',qc);
    }
  }
  const chosen=pms.filter(p=>result.selected.includes(p.id));
  LIVE.last=null;
  LIVE.source=id;LIVE.selected=chosen[0]?.key;LIVE.selectedKeys=chosen.map(p=>p.key);
  if(chosen[0])STATE.pair=chosen[0].key;
}

function multiplexPanelTSV(){
  const header=['pair_id','target','chrom','start_1based','end_1based','primer_F_5to3','primer_R_5to3','product_bp','tm_F_c','tm_R_c'];
  const selected=MPD.pms.filter(p=>MPD.result.selected.includes(p.id));
  return [header,...selected.map(p=>[p.id,p.targetName,p.chrom,p.tplOff+p.fStart+1,p.tplOff+p.rEnd+1,p.F,p.R,p.prod,p.tmF,p.tmR])].map(row=>row.join('\t')).join('\n');
}

function renderMultiplexResult(host){
  const result=MPD.result;
  if(!result){
    host.innerHTML=`<div class="card"><div class="bd mx-empty"><span class="mx-empty-icon">${ic('mpx')}</span>
      <h2>Paneliniz burada oluşacak</h2><ol class="doc-steps">
      <li>Hedef bölgeleri ve ortak koşulları girin.</li><li>Adayları üretin ve genom özgüllüğünü tarayın.</li>
      <li>Birlikte kullanılabilecek çiftleri ve etkileşim matrisini inceleyin.</li></ol>
      <a class="btn sm" href="/docs/multiplex">Multipleks örneğini incele</a></div></div>`;return;
  }
  const pairs=MPD.pms.filter(p=>result.selected.includes(p.id));
  const summary=result.summary;
  const hasPool=pairs.length>=2;
  const failed=pairs.filter(p=>p.liveQC?.records.some(r=>r.severity==='FAIL'));
  const pairQC=p=>worst((p.liveQC?.records||[]).map(r=>r.severity).filter(Boolean));
  host.innerHTML=`
    <section class="card"><div class="hd"><h2>${result.complete?'Tüm hedefleri kapsayan panel':hasPool?'Kısmi panel önerisi':'Uyumlu multipleks panel bulunamadı'}</h2></div>
      <div class="bd stack"><div class="row">${sev(result.complete?'PASS':'WARN',result.complete?'Tam panel':hasPool?'Hedeflerin bir kısmı kapsandı':'Hedefleri veya koşulları gözden geçirin')}</div>
      <div class="mini-stats"><div><div class="k">Kapsanan hedef</div><div class="v">${summary.selected} / ${summary.targets}</div></div>
      <div><div class="k">Kontrol edilen aday</div><div class="v">${summary.candidates}</div></div>
      <div><div class="k">Panel Tm aralığı</div><div class="v">${nf(summary.tm_spread,2)} °C</div></div>
      <div><div class="k">En yakın ürün boyu</div><div class="v">${nf(summary.product_gap,0)} bp</div></div></div>
      <p class="small muted">${hasPool?'Seçilen havuzda arama kapsamı içinde genomik çapraz ürün bulunmadı.':'Multipleks değerlendirme için en az iki hedefin birlikte seçilebilmesi gerekir.'} ${nf(result.elapsed_ms/1000,1)} sn · ${esc(result.profile)}</p>
      ${failed.length?liveNote(failed.length+' seçilen çiftte FAIL düzeyinde QC kaydı var. Çift İnceleme’de ilgili kuralları ve reaksiyon koşullarını kontrol edin.',true):''}
      <p class="tiny faint">Çapraz ürün kapsamı: ${result.settings.min_product}–${ni(result.settings.max_product)} bp · referans GRCh38.p14. “Tam panel”, hedeflerin tamamının havuz koşulları içinde kapsandığını belirtir.</p>
      ${result.search.limit_reached?liveNote('Seçim için hesaplama sınırına ulaşıldı; bulunan en geniş kısmi panel gösteriliyor. Tüm kombinasyonlar denenmedi.',true):''}
      ${pairs.length?`<div class="row live-actions"><button class="btn pri sm" id="mxDownload">Panel TSV indir</button>${copyBtn(multiplexPanelTSV(),'Panel TSV kopyala')}<button class="btn sm" id="mxJSON">Tüm kanıtları JSON indir</button></div>`:''}
      </div></section>
    ${result.missing.length?liveResult('Panele alınamayan hedefler',`<ul>${result.targets.filter(t=>result.missing.includes(t.id)).map(t=>`<li><b translate="no">${esc(t.name)}</b>: ${t.result.pairs.length?'Diğer seçilen hedeflerle tüm havuz koşullarını karşılayan kombinasyon bulunamadı.':'Özgüllük filtresinden geçen primer çifti bulunamadı.'}</li>`).join('')}</ul>`):''}
    ${pairs.length?liveResult('Seçilen primer çiftleri',`<div class="stack">${pairs.map(p=>`<article class="mx-pair">
      <div class="row"><span class="chip mono">${p.targetId}</span><b translate="no">${esc(p.targetName)}</b><span class="chip">Primer QC ${p.liveQC?sev(pairQC(p)):sev('NA','Değerlendirilmedi')}</span><span class="sp"></span><button class="btn sm" data-open="${esc(p.key)}">Çift incele</button></div>
      <div class="tiny faint">${esc(p.chrom)}:${gposAny(p,p.fStart)}–${gposAny(p,p.rEnd)} · ${p.prod} bp · Tm ${nf(p.tmF,1)} / ${nf(p.tmR,1)} °C</div>
      <div class="primerline"><span class="r F">F</span><span class="mx-seq">${seqC(p.F)}</span>${copyBtn(p.F)}</div>
      <div class="primerline"><span class="r R">R</span><span class="mx-seq">${seqC(p.R)}</span>${copyBtn(p.R)}</div>
      </article>`).join('')}</div>`):''}
    ${hasPool?liveResult('Primerler arası dimer matrisi',`<div class="seg" role="group" aria-label="Dimer enerjisi türü">
      <button data-mx-energy="dG37_any" aria-pressed="${MPD.mode==='dG37_any'}">Genel ΔG₃₇</button><button data-mx-energy="dG37_3prime" aria-pressed="${MPD.mode==='dG37_3prime'}">3′ uç ΔG₃₇</button></div>
      <p class="tiny muted">Hücredeki sayı kcal/mol cinsindedir. Daha düşük değer, daha güçlü dimer anlamına gelir. Çizgili hücre, seçtiğiniz dimer sınırının altındadır. “—” kararlı yapı bulunmadığını belirtir.</p>
      <div class="scrollx" id="mxMatrix"></div>`):''}
    <details class="card"><summary class="hd">Adayların birlikte seçilememe gerekçeleri · ${result.conflicts.length} kombinasyon</summary><div class="bd">${localRows(result.conflicts.map(c=>({pair:c.a+' × '+c.b,reasons:c.reasons.join(' · '),cross:c.cross_products,any:c.dG37_any,terminal:c.dG37_3prime})),[['pair','Adaylar'],['reasons','Neden'],['cross','Çapraz ürün'],['any','Genel ΔG'],['terminal','3′ ΔG']])}</div></details>
    <p class="tiny faint">Panel önerisi, hedef başına taranan adaylar ve seçilen koşullarla sınırlıdır. Her hedef için en fazla bir çift seçilir.</p>`;
  const down=$('#mxDownload',host);
  if(down)down.onclick=()=>downloadMultiplex('oligoscout-multipleks.tsv',multiplexPanelTSV(),'text/tab-separated-values;charset=utf-8');
  const json=$('#mxJSON',host);if(json)json.onclick=()=>saveLocalResult(result);
  const matrix=$('#mxMatrix',host);
  if(matrix)renderMultiplexMatrix(matrix,result);
  $$('[data-mx-energy]',host).forEach(button=>button.onclick=()=>{
    MPD.mode=button.dataset.mxEnergy;renderMultiplexResult(host);
  });
}

function renderMultiplexMatrix(host,result){
  const labels=[...new Set(result.dimer.flatMap(r=>[r.primer_a,r.primer_b]))].sort();
  const lookup=new Map();
  for(const row of result.dimer){lookup.set(row.primer_a+'|'+row.primer_b,row);lookup.set(row.primer_b+'|'+row.primer_a,row);}
  const threshold=MPD.mode==='dG37_any'?result.settings.dimer_any:result.settings.dimer_3p;
  host.innerHTML=`<table class="t mx-matrix"><caption>${MPD.mode==='dG37_any'?'Genel':'3′ uç'} dimer enerjisi, kcal/mol · F: ileri, R: geri primer</caption>
    <thead><tr><th scope="col">Primer</th>${labels.map(name=>`<th scope="col">${esc(name)}</th>`).join('')}</tr></thead>
    <tbody>${labels.map(a=>`<tr><th scope="row">${esc(a)}</th>${labels.map(b=>{
      const value=lookup.get(a+'|'+b)?.[MPD.mode];const warning=value!=null&&value<threshold;
      const intensity=value==null?0:Math.min(.35,Math.max(0,-value)/30);
      return `<td class="n ${warning?'mx-cell-flag':''}" style="background-color:color-mix(in srgb,var(--accent) ${intensity*100}%,var(--panel))" title="${esc(a+' × '+b)}${warning?' · dimer sınırının altında':''}">${warning?'! ':''}${nf(value,1)}</td>`;
    }).join('')}</tr>`).join('')}</tbody></table>`;
}
function downloadMultiplex(name,content,type){
  const url=URL.createObjectURL(new Blob([content],{type}));const link=document.createElement('a');
  link.href=url;link.download=name;link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
