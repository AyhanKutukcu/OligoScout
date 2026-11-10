# Modül 1 — NN termodinamik motoru (Python referans uygulaması)

Tek giriş noktası: **`RUN_THERMO.cmd`** (klasör kökünde, çift tıklayın).
WSL içinde `reference_py/RUN_THERMO.sh` çalışır: parametre defteri → doğrulama
paketi → demo raporu → depo paneli raporu.

Bağımlılık yok — yalnızca Python standart kütüphanesi. C++20 sürümüne
taşındığında ikisi arasında diferansiyel test çalıştırılacak.

## Dosyalar

| Yol | İçerik |
|---|---|
| `oligoscout_thermo/params.py` | Parametre tabloları, künye, `verification_status`, SHA-256 |
| `oligoscout_thermo/nn.py` | Watson-Crick NN duplex termodinamiği, dizi sözleşmesi |
| `oligoscout_thermo/salt.py` | Owczarzy 2004 / 2008, von Ahsen 2001; koşul kaydı |
| `oligoscout_thermo/concentration.py` | Üç konsantrasyon modeli, Tm, θ(T), erime eğrisi |
| `oligoscout_thermo/engine.py` | TH-01…TH-05 kontrolleri, REP-02 kanıt kaydı |
| `tests/test_thermo.py` | 28 test: analitik özdeşlikler, değişmezlikler, oracle |
| `cli/thermo_report.py` | Panel → TSV + JSONL rapor + erime eğrileri |

## Ne doğrulandı

- **28/28 test geçiyor.** Analitik özdeşlikler: her üç konsantrasyon modelinde
  θ(Tm) = 0,5 (1e-9 tolerans); ∂Tm/∂log C > 0; ∂Tm/∂[Na⁺] > 0; ters-tamamlayıcı
  değişmezliği; simetri düzeltmesinin uygulanması; NN toplanabilirliği.
- **Bağımsız oracle.** Watson-Crick tablosu ve Tm formülü, bağımsız bir
  uygulamaya (Biopython `DNA_NN3`) karşı diferansiyel testten geçti: altı dizide
  en büyük fark **0,006 °C**. Bu, transkripsiyon hatası olmadığını gösterir.
  Oracle yalnız geliştirmede kullanılır, üründe yer almaz; kurulu değilse test
  atlanır.
- **Determinizm (REP-01).** Aynı girdi → aynı JSON kaydı.

## Ne doğrulanmadı — bilerek

| Tablo | Durum | Neden |
|---|---|---|
| Watson-Crick NN (SantaLucia 1998) | `transcribed_unverified` | PNAS makalesinin PMC kaydında yalnız özet var; Tablo 1 gövdesi açık erişimde değil. Değerler literatürden aktarıldı, oracle ile sınandı, ama birincil PDF'e karşı ikinci bağımsız aktarım yapılmadı. |
| Uyumsuzluk NN (Allawi/Peyret) | `transcribed_unverified` | 48 komşu yığın parametresi birincil tablolardan aktarıldı. İkinci kişinin bağımsız aktarımı bekleniyor. |
| Sarkan uçlar (Bommarito 2000) | `transcribed_unverified` | Tablo 2’den aktarılan 32 satır OCR ile eşleşti. İkinci kişinin bağımsız aktarımı bekleniyor. |
| DNA döngü uzunluğu (SantaLucia ve Hicks 2004) | `transcribed_unverified` | Tablo 4’teki 47 enerji terimi aktarıldı. İkinci kişinin bağımsız aktarımı bekleniyor. |
| Terminal uyumsuzluk NN (SantaLucia ve Peyret 2001) | `transcribed_unverified` | WO2001094611A2 Tablo 2–3, 48 bağlam. `GG/CG` ΔH değeri aynı satırın ΔG₃₇ ve ΔS değerlerinden türetildi. |

### Primer3 çapraz denetimi

`cli/param_crosscheck.py` yerel Primer3 tablolarını salt okunur karşılaştırma
kaynağı olarak kullanır. Önceki dört tablonun 277 skaler karşılaştırması
tolerans içinde. Terminal tabloda 96 değerin 94'ü eşleşiyor. `AC/TC` ΔS
işaret farkı ve `GG/CG` ΔH türetimi `SOURCE_NOTE` olarak kaydedilir.
Bu kayıtlar bağımsız insan aktarımı bekleyen tablolara işlenir.

SPC-04, tarama kaydında referans dizisi bulunan tekli iç uyumsuzluklar için
θ(Ta) hesaplar. Tuz düzeltmesi primerin GC oranına dayanır. VAR-04, referans
bağlanma dizisi bilinen hedef yerlerde tek baz değişimleri için
ΔG₃₇(alternatif) − ΔG₃₇(referans) değerini 1 M NaCl ve pH 7 koşullarında
raporlar. Son iki bazdaki, ilk iki bazdaki ve çoklu uyumsuzluklar kapsam
dışındadır. Her iki ölçüm de INFO kaydıdır ve kalibre karar eşiği içermez.

Tabloları girdiğinizde yapılacak: `params.py` içindeki ilgili
`ParameterTable.data` alanını doldurun, `verification_status` değerini
`transcribed_unverified` yapın, iki bağımsız aktarımın diff'ini alın, sonra
`verified`e çekin. SHA-256 otomatik değişir ve her rapor satırına yazılır.

## Modelleme kararları (açıkça beyan edilir)

- **Tuz düzeltmesi ve θ(T) tutarlılığı.** Tuz düzeltmesi Tm üzerinde tanımlıdır.
  θ(T) eğrisini tutarlı kılmak için ΔH sabit tutulur ve ΔS, düzeltilmiş Tm'i tam
  verecek şekilde yeniden ölçülür. Böylece θ(Tm) = 0,5 özdeşliği tuz düzeltmesi
  altında da korunur. Kayıtta `ds_rescaled_to_salt_corrected_tm: true`.
- **TH-02 hangi modelleri karşılaştırır.** Mg²⁺ varken Owczarzy 2004 modeli
  Mg'yi hiç görmez; onu yayılıma katmak model belirsizliğini değil *koşul*
  farkını ölçer. Bu yüzden serbest Mg > 0 iken yalnız Owczarzy 2008 ve
  von Ahsen 2001 karşılaştırılır. Kayıtta `models_compared` alanı bunu gösterir.
- **TH-04 tanımı.** 3′ uçtaki son 5 bazın NN yığın ΔG°37 toplamı; başlatma
  terimi dışarıda bırakılır. İhlalin hangi tarafta olduğu (`zayif` /
  `fazla_kararli` / `asiri_kararli`) kayda yazılır.

## İlk bulgular (depo chr22 paneli, 64 primer, Ta 60 °C, 50 mM Na⁺, 1,5 mM Mg²⁺, 0,6 mM dNTP, 200 nM primer)

- Tm aralığı **56,8 – 72,9 °C**, ortanca 63,7 °C. 16 °C'lik bu yayılım paneli
  tek bir Ta ile çalıştırmayı zorlaştırır.
- θ(60 °C) aralığı **0,06 – 1,00**. 25 primer TH-03 eşiğinin (0,80) altında.
  Tm tablosu bunu göstermez; bağlı fraksiyon gösterir. Spesifikasyondaki
  "Tm yerine fiziksel anlamlı metrik" gerekçesinin somut karşılığı budur.
- Tuz modeli yayılımı **2,17 – 2,96 °C**; 64 primerin 64'ü de TH-02'nin 2 °C
  WARN eşiğini aşıyor. Owczarzy 2008 ile von Ahsen 2001 arasındaki bu sistematik
  fark gerçek bir model belirsizliğidir; eşik ya yükseltilmeli ya da tek bir
  model kanonik seçilmelidir. **Kalibrasyon adayı 1.**
- TH-04, 64 primerin 48'inde "3′ uç zayıf" tarafında WARN veriyor. Bant
  (−9…−6 kcal/mol) gerçek primerlerin çoğunu işaretliyor. **Kalibrasyon adayı 2.**

Bu eşikler QC kataloğundan gelen literatür başlangıç değerleridir ve
ıslak-lab kalibrasyonu beklemektedir; rapor bunu `calibration_status` alanında
taşır.

## Modül 2 — dimer (`align.py`) ve Modül 3 — katlanma (`fold.py`)

**Dimer (DIM-01…DIM-04, DIM-06, DIM-08).** İki oligonun tüm göreli
kaymaları taranır. Tam eşleşen helisler ve her iki yanında en az iki
Watson–Crick çifti bulunan tekli iç uyumsuzluk veya tek bazlık iç bulge
37 °C ve 1 M NaCl’de puanlanır. `3prime` modunda en az bir oligonun
terminal 3′ bazı eşleşmelidir. Hizalamada `X` iç uyumsuzluğu, `B` bulge'u
gösterir. Bulge hesabı `BULGE:1`, komşu WC yığını ve A–T kapanış
terimini kullanır. Üretim varsayılanı `nn_single_internal_mismatch_or_bulge`.
`allow_bulge=False` önceki modeli seçer. Çoklu kusur ve diziye bağlı uç
terimleri kapsam dışındadır. Kaba kuvvet ve sıra simetrisi testleri
`tests/test_dimer_fold.py` içinde. `cli/bulge_compare.py` her iki modelle
paneli puanlar ve `out/bulge` altında TSV ile PNG üretir.

**Katlanma (SS-01 ve SS-03).** DNA modeli 37 °C ve 1 M NaCl’de
Watson–Crick yığınları ile SantaLucia ve Hicks’in Tablo 4’teki döngü
uzunluğu terimlerini ve iç döngü asimetrisini kullanır. Ara uzunluklar
komşu noktalar arasında logaritmik olarak enterpole edilir. 30 bazın
üzerinde Denklem 7, 30 bazdaki tablolanan değere uygulanır. Üç bazlı
saç tokasında A–T kapanış cezası, tek bazlık bulge üzerinde komşu yığın
katkısı eklenir. Çok dallı döngüler, 1×1 iç döngüler, diğer uç terimleri
ve kısa saç tokası bonusları kapsam dışındadır.

SS-01 minimum serbest enerjiyi ve yapıyı, SS-03 her bazın marjinal
eşleşmeme olasılığını ve 3′ pencere ortalamasını `INFO` olarak raporlar.
Ölçümler referans koşullarına aittir ve kalibre karar eşiği içermez.

Doğrulama kaba kuvvettir: n ≤ 14 için tüm pseudoknot'suz yapılar sayılıp MFE,
Z ve **her baz-çifti olasılığı** DP sonucuyla 1e-9 toleransla karşılaştırılır;
geri izlemenin ürettiği yapının enerjisi MFE'ye eşit olmalı ve
−RT ln Z ≤ MFE tutmalıdır. Toplam 44 test geçiyor.

**chr22 paneli, çift düzeyi (32 çift).** Hetero-dimer ΔG°37 ortancası
−3,09 kcal/mol (en kötüsü −10,49). 32 çiftin 23'ü PAIR-01'den düşüyor
(13 FAIL, 10 WARN) — yani ileri ve geri primerlerin Tm farkı büyük. Bir çift
hem DIM-03 hem DIM-04'ten FAIL alıyor: 3′ uçtan uzatılabilir, −10,49 kcal/mol'lük
bir hetero-dimer. Bu panel arama kıyaslaması için üretilmişti, PCR başarımı için
değil; sonuçlar bunu yansıtıyor.

## Modül — hidroliz probu (`probe.py`, PRB-01…PRB-05)

TaqMan/ddPCR probu için beş kontrol ve bir prob seçici. Hepsi mevcut
Watson-Crick tablosuyla hesaplanır; bekleyen parametre tablosu gerektirmez.

| Kural | Ne ölçülür | Eşik |
|---|---|---|
| PRB-01 | Prob Tm'i − en yüksek primer Tm'i | 6–10 °C dışı WARN, <3 veya >13 °C FAIL |
| PRB-02 | 5′ uçtaki baz | G ise FAIL (florofor söndürmesi) |
| PRB-03 | C sayısı > G sayısı | değilse WARN (tümleyen ipliği kullan) |
| PRB-04 | Amplikon içinde konum | primer bağlanma yeriyle çakışma FAIL |
| PRB-05 | Prob self / prob × primer dimerleri | ΔG°37 < −6 WARN, < −9 FAIL (3′: −5 / −7) |

PRB-05'in okunuşu primerlerdekinden farklıdır: probun 3′ ucu söndürücüyle
bloke olduğu için risk uzama değil, **hibridizasyon rekabetidir** — dimer,
probu hedeften uzak tutar. PRB-05, aynı tekli iç uyumsuzluk modelini kullanır.

Prob seçici (`pick_probe`) amplikonun primer dışı bölgesinde iki ipliği de
tarar; uzunluk, GC, homopolimer, poli-G, 5′ G ve C>G filtrelerini uygular,
PRB-01 penceresinin ortasına (8 °C) yakınlığa göre sıralar ve kalanlar
arasında dimer cezası uygular. Her elenme gerekçesiyle sayılır.

Çalıştırma — tasarımla aynı komutta:

```
RUN_DESIGN.cmd hedef_chr22.fa ddpcr          (profil probu kendiliğinden seçer)
RUN_DESIGN.cmd hedef_chr22.fa qpcr_sybr prob (prob adımını zorlar)
```

Çıktılar `probe_panel.tsv` (seçilen problar), `probe_qc.tsv` (satır başına bir
kontrol) ve `probe_qc_records.jsonl` (REP-02 kanıt kaydı). Hazır probları
denetlemek için `cli/probe_qc.py --probes pair_id<TAB>probe` dosyası verilir.

chr22 hedefinde ölçülen: qPCR profili (90–160 bp) 20 çiftin 20'sinde prob
buldu, ΔTm 7,87–8,25 °C, en kötü dimer −4,69 kcal/mol, hepsi INFO. ddPCR
profili (60–120 bp) 20 çiftin 18'inde prob buldu; kalan 2 çiftte primerler
arası bölge 13 nt'lik en kısa proba bile yetmiyor ve bu "PROB YOK" olarak
raporlanıyor. 19 test (`tests/test_probe.py`), toplam 92 test geçiyor.

## Modül — varyant katmanı (`variants.py`, VAR-01…VAR-08)

Primer bağlanma yerine düşen bir varyant o alelde bağlanmayı zayıflatır ya da
bitirir: **alel düşmesi**. Etki, varyantın primer üzerindeki konumuna ve
popülasyondaki sıklığına bağlıdır.

| Kural | Ne ölçülür | Eşik |
|---|---|---|
| VAR-01 | bağlanma yerindeki herhangi bir varyant | AF ≥ 10⁻³ WARN |
| VAR-02 | 3′ terminal 5 nt içindeki varyant | AF ≥ 10⁻⁴ WARN, ≥ %1 FAIL |
| VAR-03 | bağlanma yerinde yaygın varyant | AF ≥ %1 WARN, ≥ %5 FAIL |
| VAR-04 | varyantın ΔΔG₃₇ etkisi | Referansı bilinen hedef yerlerde desteklenen tek baz değişimleri |
| VAR-05 | prob bölgesinde varyant | AF ≥ 10⁻³ WARN, ≥ %1 FAIL |
| VAR-06 | bağlanma yeriyle örtüşen indel | AF ≥ 10⁻⁴ WARN, ≥ %1 FAIL |
| VAR-07 | popülasyona özgü yüksek frekans | pop AF ≥ %1 ve ≥ 5× küresel AF |
| VAR-08 | amplikon içinde, primer dışı yaygın indel | AF ≥ %1 WARN |

Eşikler literatür başlangıç değerleridir ve `VariantConstraints` alanlarından
değiştirilir; QC kataloğunun kendi değerleri geldiğinde kural kimlikleri ve
kayıt şeması değişmeden hizalanır.

İki tasarım kararı: (1) `variants.py` **ağa çıkmaz** — varyantları dışarıdan
alır, böylece referans katmanı standart kütüphane dışına çıkmaz ve testler
çevrimdışı çalışır; veri çekme `cli/variant_qc.py` içindedir (gnomAD v4
GraphQL) ve çekilen her bölge `variant_cache.jsonl`'e yazılır, `OFFLINE=1` ile
ağsız çalışılır. (2) 3′ pencere iplik yönüne göre hesaplanır — `+` ipliğinde
sağ, `−` ipliğinde sol uç; yanlış ucu almak VAR-02'yi anlamsız kılardı ve
testlerde ayrıca doğrulanır.

Çalıştırma:

```
RUN_VARIANT.cmd                                               ctDNA paneli
RUN_VARIANT.cmd reference_py\out\panel_design_qpcr_sybr ^
                reference_py\out\design_qpcr_sybr
```

chr22 hedefinde ölçülen: hedef bölge iki yaygın SNP taşıyor — **rs2106140
(AF %37,5)** ve **rs4819852 (AF %22,1)**. ctDNA panelinde 20 çiftin 12'si
VAR FAIL, 5'i WARN; yalnız 3 çift (DESIGN_010/017/018) varyantsız. qPCR
panelinde hiçbir çift varyantsız değil. Daha çarpıcısı: özgüllük taramasında
PASS alan çiftlerin 9'u (ctDNA) ve 6'sı (qPCR) varyant katmanında FAIL alıyor
ve her iki panelde de **1. sıradaki çift yaygın bir SNP'nin üzerinde oturuyor**.
Varyant katmanı eklendiğinde `tests/test_variants.py` 24, referans paket
toplam 116 test içeriyordu. Bulge ve CT-04 kapanışında tam paket **166 test**
ile geçti.

## Sonraki adımlar

CT-04, ctDNA için mevcut. `thermo_report.py --ct04-alleles` girdisi
`primer_id`, `primer`, `strand`, `matched_plus` ve `mismatched_plus`
sütunlarını ister. İki alel penceresini genomun artı ipliğinde 5′→3′ yazın.
Terminal bağlam için WO2001094611A2, sondan ikinci baz için iç NN
parametreleri kullanılır. θ(Ta), θ marjı ve ΔΔG₃₇ `ct04.tsv`,
`ct04_records.jsonl` ve ana JSONL raporunda yazılır. REP-03 kapsamında
sayısal marjlar kalibre edilmemiş proje başlangıç değerleridir.

`cli/threshold_evidence.py --latest <selected_panel.tsv>` TH-02 ve TH-04
için karar desteği üretir. Eşik seçenekleri çevrimdışı raporda hesaplanır.

1. Parametre tablolarını ikinci bir kişinin birincil kaynaklardan bağımsız aktarması.
2. TH-02, TH-04 ve CT-04 başlangıç eşiklerinin deneysel verilerle değerlendirilmesi.
3. ~~Termodinamik katmanın C++20'ye taşınması ve Python ↔ C++ diferansiyel
   testi.~~ **Tamamlandı.**

   `RUN_DIFF_CPP.cmd` `oligoscout_thermo_dump` hedefini derler, 74 dizilik
   listeyi iki uygulamadan geçirir ve `cli/diff_cpp.py` ile sütun sütun
   karşılaştırır. Son koşu: **1.258 karşılaştırmanın 1.258'i eşleşiyor**
   — 12 çekirdek büyüklük (ΔH/ΔS 1e-13 düzeyinde; Tm, θ(Ta), ΔG₃₇, 3′ uç
   kararlılığı, tuz rejimi, konsantrasyon modeli) ve 5 dimer sütunu
   (self/hetero, `any` ve `3prime`). Dimer farkı artık model boşluğu
   sayılmaz; iki taraf aynı yapı sınıfını tarar ve her fark hatadır.

   Uyumsuzluk tablosu ve bulge sabitleri C++'a elle yazılmaz:
   `tools/gen_thermo_tables.py` bunları doğrudan `params.py`'den üretip
   `include/primerpair/generated/thermo_tables.inc` dosyasına yazar,
   `--check` modu da başlığın güncelliğini denetler (regresyon adımı 16b).

   Bulge'u kapatan A/T çifti başına uygulanan 0,5 kcal/mol'lük ceza artık
   `params.MODEL_CONSTANTS['at_bulge_closure_dg37']` içinde kayıtlıdır.
   Bu kayıt parametre tablolarından **ayrıdır**: literatürden aktarılmış bir
   ölçüm değil, beyan edilmiş bir modelleme kararıdır (SantaLucia & Hicks'in
   A·T kapanış cezası saç tokası döngüsü için tanımlıdır; bulge kapanışına
   uygulamak bu projenin kararıdır), bu yüzden `verification_status` değil
   `rationale` ve `status` taşır. Değeri `registry_digest()`e girer — yani
   değiştirilirse her kayıttaki `param_set_sha256` değişir — ve dimer
   `model_metadata()` çıktısında beyan edilir. `align.py`, test paketindeki
   kaba kuvvet sayımı ve C++ tablo üreticisi aynı tek kaynaktan okur.
   Sayısal değer **kalibre değildir** (REP-03).
4. chr3 hedef varyant konumunun belirlenmesi ve ona göre aday seçilmesi.
