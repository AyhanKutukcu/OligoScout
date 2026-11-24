const QC_SCOPE_NOTES = {
  'SPC-08': 'Yerel iş akışı birincil 24 kromozom dışındaki tüm GRCh38.p14 dizilerini tarar: chrM, 260 ALT lokus, 254 yama ve 166 iskele (rDNA kopyaları dahil). Ayrı 45S rDNA referansı (KY962518) derlemede yoktur; ALT/yama bağlanmaları alternatif haplotip kopyası olarak bilgi düzeyinde raporlanır.',
  'RT-01': 'Bu kurulumda ekson birleşimi bir profil önerisi olarak kayıtlı. Birleşimin iki yanını ölçen ayrı bir sonuç doğrulanamadı.',
  'RT-02': 'Profil ayrı bir genomik DNA taraması önerir. RT profiline özgü bir sonuç kaydı doğrulanamadı.',
  'RT-04': 'RNA katlanması için gereken parametreler profil kaydında eksik olarak belirtiliyor.',
  'SAN-01': 'Sanger profilinde hedefe uzaklık önerisi bulunur. Otomatik bir uzaklık kontrolü doğrulanamadı.',
  'TH-06': 'Katkı maddelerine ait düzeltme katsayıları bu kurulumda kullanılmıyor.',
  'TH-07': 'Bu tahmin ancak PCR ürününün erime sıcaklığı mevcutsa gösterilir.',
};

const catalogKey = value => value.normalize('NFKD').toLowerCase().replace(/\p{Mark}/gu, '').replace(/ı/g, 'i');

RENDER.katalog = el => {
  const partial = CAT.filter(r => r.status.startsWith('etiketli_kismi')).length;
  const pending = CAT.length - CAT_FULL - partial;
  const filters = [['all', 'Tümü', CAT.length], ['full', 'Tam kapsam', CAT_FULL], ['part', 'Kısmi kapsam', partial]];
  if (pending) filters.push(['pend', 'Bekleyen', pending]);
  el.innerHTML = `<div class="stack qc-catalog">
    <div class="g3">${kpi('Toplam kural', CAT.length, '', FAMORDER.length + ' aile')}${kpi('Tam kapsam', CAT_FULL, '', 'Proje kaydına göre')}${kpi('Kısmi kapsam', partial, '', 'Proje kaydına göre')}</div>
    <div class="card"><div class="bd qc-controls">
      <p class="small muted">Bir kuralı açarak neyi kontrol ettiğini ve neden önemli olduğunu görün. Kapsam, bir primerin test sonucunu göstermez.</p>
      <div class="row"><div class="seg" id="cSeg" aria-label="Kapsama göre filtrele">${filters.map(([key, label, count]) => `<button type="button" data-f="${key}"><span>${label}</span> <span class="mono">${count}</span></button>`).join('')}</div>
        <input type="search" id="cQ" aria-label="Katalogda ara" placeholder="Kural kodu veya konu ara…" value="${esc(STATE.catQ || '')}"></div>
      <p id="cCount" class="small muted" aria-live="polite"></p>
    </div></div>
    <div id="cList" class="stack"></div>
  </div>`;
  const list = $('#cList', el);
  const input = $('#cQ', el);
  const buttons = $$('#cSeg button', el);
  const refresh = () => {
    const filter = STATE.catF || 'all';
    const query = catalogKey(input.value.trim());
    const matches = CAT.filter(r => {
      const kind = r.status === 'etiketli_uygulandi' ? 'full' : r.status.startsWith('etiketli_kismi') ? 'part' : 'pend';
      if (filter !== 'all' && filter !== kind) return false;
      const terms = [r.id, FAMTR[r.fam], r.name, r.check, r.why];
      return !query || catalogKey(terms.flatMap(t => [t, englishMessages.get(t) || t]).join(' ')).includes(query);
    });
    buttons.forEach(b => b.setAttribute('aria-pressed', b.dataset.f === filter));
    $('#cCount', el).textContent = `${matches.length} / ${CAT.length} kural gösteriliyor`;
    list.innerHTML = FAMORDER.map(fam => {
      const rules = matches.filter(r => r.fam === fam);
      if (!rules.length) return '';
      return `<section class="card qc-family"><div class="hd"><h2>${esc(FAMTR[fam])}</h2><span class="mono tiny muted">${rules.length}</span></div>
        <ul class="qc-rules">${rules.map(r => {
          const status = STATUS[r.status] || {c:'NA', t:'Bekleyen'};
          return `<li><button type="button" class="qc-rule-row" data-rule="${r.id}" aria-haspopup="dialog">
            <span class="mono qc-code">${r.id}</span><span class="qc-rule-copy"><strong>${esc(r.name)}</strong><span>${esc(r.check)}</span></span>
            <span class="qc-rule-status">${sev(status.c, status.t)}</span><span class="qc-chevron" aria-hidden="true">›</span>
          </button></li>`;
        }).join('')}</ul></section>`;
    }).join('') || '<p class="note">Aramanıza uyan kural bulunamadı.</p>';
  };
  input.oninput = () => { STATE.catQ = input.value; refresh(); };
  buttons.forEach(button => button.onclick = () => { STATE.catF = button.dataset.f; refresh(); });
  refresh();
};
