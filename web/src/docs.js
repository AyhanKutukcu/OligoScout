
/* Usage guide and the archived, read-only result examples. */
const DOC_PAGES = [
  ['rehber', 'Kullanım rehberi'], ['panel', 'Kayıtlı Panel'], ['cift', 'Çift İnceleme örneği'],
  ['multipleks', 'Multipleks'], ['prob', 'Prob & Blocker'], ['erime', 'qPCR Erime'], ['dogrulama', 'Motor doğrulaması'],
];
const DOC_REGION = 'chr22:20000551-20000870';
const DOC_PROJECT_REGION = CHR + ':' + (OFF + 1) + '-' + (OFF + TPL.length);
const DOC_FASTA = '>chr22:20000551-20000870 GRCh38.p14\n' +
  TPL.slice(20000551 - 1 - OFF, 20000870 - OFF).match(/.{1,60}/g).join('\n');
const DOC_TOOL_STEPS = {
  multipleks:['Multipleks Tasarla’da en az iki hedef adı ve GRCh38 bölgesi girin. Varsayılan Sanger profili 300–800 bp ürün aralığı kullanır.', 'Hedef başına aday sayısını ve aynı havuza kabul koşullarını seçip tasarımı başlatın.', 'Her hedef genomda taranır. Birlikte kullanılabilecek çiftler seçilir; kapsanmayan hedefler varsa belirtilir. Paneli TSV olarak indirebilirsiniz. Aşağıdaki arşiv, eski çalışmanın dimer ve çapraz ürün raporunu okumak için bir örnektir.'],
  prob:['Yerel Analiz’de Hidroliz probu veya Prob ve blocker QC’yi seçin.', 'Hidroliz probu için bir çift seçin ve varsa prob dizisini girin. Blocker QC için mutant ve WT prob dizileri ayrı alanlara girilir.', 'Çıktı Tm, dimer ve ilgili kural kayıtlarını gösterir. Aşağıdaki örnekler kaydedilmiş ölçümleri içerir.'],
  erime:['Yerel Analiz’de Amplikon erime profili’ni seçin.', 'Tasarım kaynağını ve primer çiftini seçip analizi çalıştırın.', 'Sıcaklık eğrisindeki pikleri ve ilgili QC kaydını okuyun. Bu örnek, çiftlerin erime piklerini aynı grafikte karşılaştırır.'],
};
const DOC_RESULT_TABS = [
  ['harita', 'Harita', 'F ve R okları bağlanma yönlerini, iki okun kapsadığı alan amplikonu gösterir.'],
  ['ozgulluk', 'Özgüllük', 'Beklenen ürün ile hedef dışı ürünleri ve kromozomlara göre bağlanma yerlerini karşılaştırın.'],
  ['varyant', 'Varyant', 'İşaretler varyant konumlarını gösterir. AF alel frekansıdır. Canlı sonuç tablosunda genom ve ekzom frekansları ayrı sütunlarda gösterilir.'],
  ['yapi', 'Yapı & Dimer', 'F×F, R×R ve F×R hizalamalarını, 3′ uç eşleşmelerini ve saç tokası çizimlerini inceleyin. ΔG₃₇ kcal/mol cinsindendir.'],
  ['termo', 'Termodinamik', 'İleri ve geri primerin sıcaklığa göre bağlanma eğrilerini karşılaştırın. Tm erime sıcaklığı, Ta reaksiyonun bağlanma sıcaklığıdır.'],
  ['kanit', 'Kanıt', 'Kural kimliğini açarak ölçümü, eşiği ve sonucu üreten motoru görün. NA, kuralın değerlendirilmediğini belirtir.'],
];

RENDER.dokumantasyon = el => {
  const section = STATE.docSection;
  const [, title] = DOC_PAGES.find(([id]) => id === section);
  el.innerHTML = `<div class="stack">
    <nav class="card doc-nav" aria-label="Dokümantasyon bölümleri">${DOC_PAGES.map(([id, label]) => {
      const href = routeUrl(id === 'rehber' ? 'dokumantasyon' : 'dokumantasyon/' + id);
      return `<a href="${href}" ${id === section ? 'aria-current="page"' : ''}>${label}</a>`;
    }).join('')}</nav>
    <div id="docContent"></div>
  </div>`;
  const content = $('#docContent', el);
  if (section === 'rehber') { renderGuide(content); return; }
  content.className = 'stack';
  content.innerHTML = `<div class="note">${ic('info')}<div><b>Arşiv örneği · ${esc(title)}</b><br><span class="small">Bu sayfa daha önce kaydedilmiş bir çalışmayı gösterir. Kendi sonuçlarınız için Primer Tasarla veya Yerel Analiz’i kullanın.</span></div></div>
    ${section === 'multipleks' ? '<a class="btn pri" href="/multiplex">Kendi multipleks panelini tasarla</a>' : ''}
    ${DOC_TOOL_STEPS[section] ? liveResult('Girişten çıktıya', `<ol class="doc-steps">${DOC_TOOL_STEPS[section].map(step => `<li>${esc(step)}</li>`).join('')}</ol>`) : ''}<div></div>`;
  RENDER[section === 'multipleks' ? 'multipleksOrnek' : section](content.lastElementChild);
};

function renderGuide(el) {
  const pair = PROJ[0];
  const examplePairLink = routeUrl('dokumantasyon/cift/' + pair.key);
  el.innerHTML = `<div class="stack">
    <section class="card"><div class="bd intro-hero">
      <div class="eyebrow">OLIGOSCOUT KULLANIM REHBERİ</div>
      <h2>Girişten sonuçlara primer tasarımı</h2>
      <p>Bir insan genom bölgesi seçin, tasarımı çalıştırın ve otomatik özgüllük taramasından geçen çiftleri inceleyin. Aşağıdaki örnekleri formu doldurmak veya bir sonucun nasıl göründüğünü anlamak için kullanabilirsiniz.</p>
      <div class="row live-actions"><a class="btn pri" href="/design">Primer Tasarla</a><a class="btn" href="/docs/panel">Kayıtlı Paneli aç</a></div>
    </div></section>

    <section class="card"><div class="hd"><h2>1 · Bir giriş biçimi seçin</h2></div><div class="bd stack">
      <p class="muted">Genomik bölge, yapıştırılan FASTA ve yüklenen FASTA aynı hedefi tarif etmenin üç yoludur. Yerel motor GRCh38.p14 üzerinde 1–22, X ve Y kromozomlarında 1–5.000 bazlık bölgeleri kabul eder.</p>
      <div class="g2">
        <article class="doc-box">
          <h3>Genomik bölge</h3>
          <p>Koordinatları 1 tabanlı yazın; ilk ve son baz dahildir. Bu örnek 320 bp seçer.</p>
          <div class="row"><code>${DOC_REGION}</code>${copyBtn(DOC_REGION, 'Koordinatı kopyala')}</div>
          <button type="button" class="btn sm" data-doc-input="bolge">Koordinatı formda aç</button>
        </article>
        <article class="doc-box">
          <h3>FASTA metni veya dosyası</h3>
          <p>Tek kayıt kullanın. Başlık genom koordinatıyla başlamalı, dizinin uzunluğu aralıkla, bazları ise referansla uyuşmalıdır. Dosya olarak <code>.fa</code>, <code>.fasta</code>, <code>.fna</code> veya <code>.txt</code> seçebilirsiniz.</p>
          <div class="row"><button type="button" class="btn sm" data-doc-input="fasta">FASTA’yı formda aç</button>${copyBtn(DOC_FASTA, 'FASTA metnini kopyala')}</div>
        </article>
      </div>
      <details class="doc-box"><summary>Tam FASTA örneğini göster</summary><pre>${esc(DOC_FASTA)}</pre></details>
      <details class="doc-box"><summary>Eski “Proje şablonu” hangi bölgeyi kullanıyordu?</summary>
        <p>Kaydedilmiş chr22 panelinin şablonu <code>${esc(D.meta.header)}</code> bölgesidir. İlk örnek, bu şablonun 320 bazlık bir kesitini kullanır.</p>
        <div class="row">
          <button class="btn sm" type="button" data-doc-input="project">Tam şablonu formda aç</button>
          <a class="btn sm" href="/docs/panel">Kayıtlı şablon ve sonuçları incele</a>
        </div>
      </details>
    </div></section>

    <section class="card"><div class="hd"><h2>2 · Tasarımı çalıştırın</h2></div><div class="bd">
      <ol class="doc-steps">
        <li><b>Profilinizi seçin.</b> Bu örnekte ctDNA profilinin varsayılan ayarlarıyla başlayabilirsiniz. Ürün uzunluğu ve primer koşulları, seçtiğiniz profile göre doldurulur.</li>
        <li><b>İsteğe bağlı aralıkları girin.</b> Şablonun ilk bazı <code>1</code> olur. Örneğin hedef bölgenin başı <code>141</code>, sonu <code>160</code> ise amplikonun bu aralığı kapsaması istenir. Konum kısıtı istemiyorsanız alanları boş bırakın.</li>
        <li><b>Primer Tasarla’ya basın.</b> Aday üretimi, gnomAD kontrolü (veri varsa), genom özgüllüğü taraması ve seçilen çiftlerin QC analizi aynı işlemde yürütülür.</li>
        <li><b>Sonuç kartını açın.</b> Uygun çiftin F/R dizileri, ürün boyu ve Tm değerleri kartta görünür. Paneli TSV olarak kopyalayabilir veya Çift İnceleme’ye geçebilirsiniz.</li>
      </ol>
    </div></section>

    <section class="card"><div class="hd"><h2>Birden fazla hedef için multipleks tasarımı</h2></div><div class="bd stack">
      <p class="muted">Multipleks Tasarla’ya birden fazla bölge girin. Her hedef için en fazla bir çift seçilir; tüm primerler arasındaki dimerler ve genomik çapraz ürünler birlikte kontrol edilir.</p>
      <div class="doc-box"><pre>Hedef 1  chr22:20000001-20001000
Hedef 2  chr22:20004001-20005000</pre><button class="btn sm" type="button" id="docMultiplexExample">Bu hedefleri multipleks formunda aç</button></div>
      <p class="small muted">Hedef başına aday sayısı, o bölgeden kaç özgül çiftin panel seçimine gireceğidir. “En fazla Tm aralığı”, paneldeki en yüksek ve en düşük primer Tm farkını sınırlar. Ürün boyu farkı 0 ise bu koşul uygulanmaz. Panel tüm hedefleri kapsamayabilir; sonucu ve gerekçeleri kontrol edin.</p>
      <a class="btn sm" href="/docs/multiplex">Kayıtlı etkileşim raporunu incele</a>
    </div></section>

    <section class="card"><div class="hd"><h2>3 · Bir çıktıyı okuyun</h2><span class="tag">Kayıtlı örnek</span></div><div class="bd stack">
      <div class="doc-box">
        <div class="row"><b class="mono">${esc(pair.id)}</b><span class="chip">Ürün <b>${pair.prod} bp</b></span>
          <span class="chip">Tm F / R <b>${nf(pair.tmF,1)} / ${nf(pair.tmR,1)} °C</b></span>
          <span class="chip">Karar ${sev(pair.verdict)}</span><a class="btn sm" href="${examplePairLink}/harita">Tam örnek çıktıyı aç</a></div>
      </div>
      <div class="doc-results">${DOC_RESULT_TABS.map(([id, label, desc]) => `<a href="${examplePairLink}/${id}"><b>${label}</b><span>${desc}</span></a>`).join('')}</div>
      <p class="small muted">Multipleks, Prob & Blocker ve qPCR Erime örnekleri üstteki sekmelerdedir. Kendi primerlerinizle bu araçları <a href="/analysis">Yerel Analiz</a> sayfasından kullanabilirsiniz.</p>
    </div></section>

    <section class="card"><div class="hd"><h2>Sık karşılaşılan durumlar</h2></div><div class="bd stack">
      <details class="doc-box"><summary>Yerel motor kapalı görünüyor</summary><p>Proje klasöründeki <code>RUN_WEB.cmd</code> dosyasını çalıştırıp <code>http://127.0.0.1:8765</code> adresini açın.</p></details>
      <details class="doc-box"><summary>Özgüllük filtresinden geçen çift bulunamadı</summary><p>Sonuçlardaki elenme gerekçelerini açın. Ön filtrelerde kalan aday yoksa giriş aralığı ve tasarım kısıtlarını kontrol edin. Adaylar genom taramasında elendiyse hedef bölgeyi yeniden değerlendirin.</p></details>
      <details class="doc-box"><summary>Varyant verisi neden yok?</summary><p>“Varyant bulunmadı”, sorgunun tamamlandığını ve bu bölgede kayıt dönmediğini gösterir. “Değerlendirilmedi” ise veri kaynağının bölgeyi kapsamaması veya isteğin tamamlanamaması gibi farklı bir durumu gösterir.</p></details>
    </div></section>
  </div>`;
  $('#docMultiplexExample', el).onclick = () => {
    MPD.targets=[{name:'Hedef 1',region:'chr22:20000001-20001000'},{name:'Hedef 2',region:'chr22:20004001-20005000'}];
    MPD.profile='sanger';MPD.product_min=PROFILES.sanger.product_bp[0];MPD.product_max=PROFILES.sanger.product_bp[1];
    go('multipleks');
  };
  $$('[data-doc-input]', el).forEach(button => button.onclick = () => {
    const source = button.dataset.docInput;
    Object.assign(DSG, {profile:'ctdna', src:source === 'fasta' ? 'fasta' : 'bolge',
      region:source === 'project' ? DOC_PROJECT_REGION : source === 'bolge' ? DOC_REGION : '',
      fasta:source === 'fasta' ? DOC_FASTA : '', fname:''});
    go('tasarla');
  });
}

const CONTACT_EMAIL = 'bioinformaticianayhan@gmail.com';
const CONTACT_TEAM = [
  {name:'Ayhan KÜTÜKÇÜ', role:'Biyoinformatik araştırmacısı', links:[
    ['LinkedIn', 'https://www.linkedin.com/in/ayhan-kutukcu-434a69235'],
    ['E-posta', 'mailto:' + CONTACT_EMAIL],
  ]},
  {name:'Güzin Tunca Alparslan', role:'Dr. Öğr. Üyesi', links:[
    ['LinkedIn', 'https://www.linkedin.com/in/guzin-tunca-5a001848'],
    ['Akademik profil', 'https://personel.trakya.edu.tr/guzintunca/'],
    ['E-posta', 'mailto:guzintunca@trakya.edu.tr'],
  ]},
  {name:'Vuslat Berfin SAKIZCI'},
];

RENDER.iletisim = el => {
  el.innerHTML = `<div class="stack">
    <section class="card"><div class="bd intro-hero">
      <h2>Yardım ve iletişim</h2>
      <p>OligoScout ekibiyle tanışın, kampüsün konumunu görün veya bir soru iletin.</p>
    </div></section>
    <section class="stack" aria-labelledby="contactTeam">
      <h2 id="contactTeam">Proje ekibi</h2>
      <div class="g3 contact-team">
        ${CONTACT_TEAM.map(person => `<article class="card"><div class="bd">
          ${person.role ? `<p class="small muted">${esc(person.role)}</p>` : ''}
          <h3>${esc(person.name)}</h3>
          ${person.links ? `<div class="row contact-links">${person.links.map(([label, url]) =>
            `<a class="btn sm" href="${esc(url)}" aria-label="${esc(person.name + ' · ' + label)}" ${url.startsWith('https:') ? 'target="_blank" rel="noopener noreferrer"' : ''}>${esc(label)}</a>`
          ).join('')}</div>` : ''}
        </div></article>`).join('')}
      </div>
    </section>
    <div class="g2 contact-details">
      <section class="card"><div class="hd"><h2>Kampüs</h2></div><div class="bd stack">
        <address><b>Trakya Üniversitesi</b><br>Ahmet Karadeniz Yerleşkesi<br>22020 Merkez / Edirne, Türkiye</address>
        <iframe class="contact-map" title="Ahmet Karadeniz Yerleşkesi haritası" loading="lazy" referrerpolicy="no-referrer"
          src="https://www.google.com/maps/embed?pb=!1m14!1m8!1m3!1d2981.9862194871325!2d26.6206061!3d41.6344293!3m2!1i1024!2i768!4f13.1!3m3!1m2!1s0x14b33aba49aba687%3A0xeff543895dfd7b2f!2sProfes%C3%B6r%20Dr.%20Ahmet%20Karadeniz%20Yerle%C5%9Fkesi%20M%C3%BChendislik%20Mimarl%C4%B1k%20Fak%C3%BCltesi!5e0!3m2!1str!2str!4v1707029687012!5m2!1str!2str"></iframe>
        <div><a class="btn sm" href="https://www.google.com/maps/search/?api=1&amp;query=41.6344293%2C26.6206061" target="_blank" rel="noopener noreferrer">Haritada aç ↗</a></div>
      </div></section>
      <section class="card"><div class="hd"><h2>Mesaj yazın</h2></div>
        <form class="bd stack" id="contactForm">
          <label class="f"><span>Ad soyad</span><input type="text" name="name" autocomplete="name" maxlength="120" required></label>
          <label class="f"><span>E-posta</span><input type="email" name="email" autocomplete="email" maxlength="254" required></label>
          <label class="f"><span>Mesaj</span><textarea name="message" rows="6" maxlength="1500" required></textarea></label>
          <p class="small muted">Taslak, e-posta uygulamanızda açılır. Mesajı oradan gönderebilirsiniz.</p>
          <button class="btn pri" type="submit">E-posta taslağı aç</button>
        </form>
      </section>
    </div>
  </div>`;
  $('#contactForm', el).addEventListener('submit', event => {
    event.preventDefault();
    const fields = new FormData(event.currentTarget);
    const body = `${fields.get('message')}\r\n\r\n${fields.get('name')}\r\n${fields.get('email')}`;
    location.href = `mailto:${CONTACT_EMAIL}?subject=${encodeURIComponent(translateMessage('OligoScout iletişim'))}&body=${encodeURIComponent(body)}`;
  });
};
