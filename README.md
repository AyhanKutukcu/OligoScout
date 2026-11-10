# OligoScout

GRCh38 insan referansi uzerinde **PCR primer/prob tasarimi ve kapsamli kalite
kontrolu (QC)** yapan, sifirdan yazilmis bir hesaplamali termodinamik sistem.
Hicbir dis termodinamik kutuphaneye bagli degildir: en-yakin-komsu (NN) modeli,
dimer hizalama, ikincil yapi katlanmasi, tuz duzeltmeleri ve genom-geneli FM-index
aramasi projenin kendi kodudur.

| | |
|---|---|
| **QC kurallari** | 19 ailede 91 kural; 90'i uygulandi |
| **Birim test** | 264 (hepsi geciyor) |
| **Uctan uca regresyon** | 20 PASS / 1 SKIP / 0 FAIL |
| **Python ↔ C++ diferansiyel** | 74 dizi × 17 sutun, 0 fark |
| **Parametre capraz denetim** | 373 karsilastirma (primer3 oracle) |
| **C++ tarama (yayin paneli)** | 1.024 cift, 8.262.803 baglanma, 111.8 s, 0.65 GiB |

> **Dikkat — arastirma/bilgi amaclidir.** Esiklerin cogu literatur baslangic
> degerleridir ve **islak-lab verisiyle KALIBRE EDILMEMISTIR (REP-03)**. Ciktilar
> deneysel dogrulama yerine gecmez.

---

## Mimari — iki katman

1. **`reference_py/` — saf Python referans motoru.** Tum termodinamik ve QC
   mantigi burada tanimli; calistirilabilir ve dogrulanabilir "referans gercek".
2. **`OligoScout/` — C++20 cekirdek.** Genom-geneli FM-index (PPFM) aramasi ve
   yuksek basarimli tarama; WSL/Linux gerektirir. Python katmani onun
   spesifikasyonu ve dogrulama oracle'idir.

Bir hesap once Python'da tanimlanir ve test edilir, sonra C++'a tasinir; iki
uygulama **diferansiyel test** ile tolerans dahilinde karsilastirilir.

---

## Ozellikler

- **Termodinamik model:** Watson-Crick NN (SantaLucia 1998), ic/terminal
  uyumsuzluk, sarkan uc, DNA dongu/hairpin katlanmasi (Zuker MFE + McCaskill),
  3 tuz modeli, DNA/RNA hibrit NN (Sugimoto 1995).
- **QC kataloğu (19 aile, 91 kural):** dizi butunlugu, kompozisyon, termodinamik,
  ikincil yapi, dimer, cift uyumu, ozgulluk, amplikon, prob, ctDNA, mRNA/RT,
  qPCR, bisulfit, multipleks, varyant, reaksiyon kosulu, NGS, Sanger, rapor.
- **8 uygulama profili:** ctDNA, ddPCR, RT-qPCR, SYBR qPCR, Sanger, NGS amplikon,
  bisulfit, standart PCR.
- **Primer tasarimi:** aday uretimi, sert filtreler, Pareto siralama, ozgullukle
  yeniden siralama.
- **Ozgulluk:** genom-geneli FM-index aramasi, hedef disi urunler, 3'-capa,
  psodogen/segdup ortusmesi (GENCODE/UCSC), arama sozlesmesi + indeks checksum.
- **Ileri katmanlar:** varyant (gnomAD v4), multipleks capraz-reaktivite, bisulfit
  4-iplik ozgullugu (donusturulmus genom indeksi), izoform kapsami (GENCODE ekson).
- **Modifiye oligolar:** LNA ve uc modifikasyon notasyonu ayristirilir (termodinamik
  DNA iskeleti uzerinden; LNA katkisi modellenmez ve acikca raporlanir).
- **Kanit ureten QC:** her kontrol bir kanit kaydi uretir (deger, esik, onem
  derecesi, model, kosullar, parametre SHA-256). Tekrar uretilebilir.

---

## Depo duzeni

```
reference_py/
  oligoscout_thermo/   termodinamik + QC motoru (nn, align, fold, salt, seqpair,
                       probe, variants, multiplex, bisulfite, isoform, params ...)
  cli/                 komut satiri araclari (21 arac)
  tests/               birim testler (264)
OligoScout/
  src/ include/        C++20 FM-index cekirdek
  tests/ tools/ scripts/ docs/
setup/                 WSL kurulum/boru hatti betikleri (*.sh)
RUN_*.cmd              Windows sarmalayicilari (WSL'yi cagirir)
README_KURULUM.md      ayrintili kurulum
```

> **Veri depoda YOKTUR.** Genomlar, indeksler, anotasyonlar ve uretilen ciktilar
> `.gitignore` ile haric tutulur (asagiya bakin).

---

## Veri dahil DEGILDIR — nasil elde edilir

Depo yalniz kaynak kodu icerir. Asagidaki veriler yereldir ve `.gitignore`
tarafindan dislanir:

| Veri | Kaynak | Yol (yerel) |
|---|---|---|
| GRCh38.p14 referans | NCBI RefSeq | `OligoScout/data/genomes/` |
| PPFM FM-index shard'lari | yerel kurulum (`RUN_SETUP.cmd`) | `OligoScout/data/indexes/ppfm/` |
| GENCODE v47 anotasyon (GTF) | GENCODE | `OligoScout/data/annotations/` |
| UCSC genomicSuperDups | UCSC Table Browser | `OligoScout/data/annotations/` |
| Bisulfit donusturulmus indeksler | `RUN_BISULFITE_INDEX.cmd` | `OligoScout/data/indexes/bisulfite/` |
| primer3 (parametre capraz denetimi icin) | github.com/primer3-org/primer3 | `OligoScout/third_party/primer3/` |

Referans ve indeksler `setup/` betikleriyle yerelde kurulur; detay
`README_KURULUM.md` ve `RUN_SETUP.cmd` icindedir.

---

## Kurulum

**Python referans motoru** (dis bagimlilik gerektirmez; standart kutuphane):

```bash
cd reference_py
python -m unittest discover -s tests -p "test_*.py"   # 264 test
```

**C++ cekirdek** (WSL/Linux, CMake + C++20 derleyici):

```bash
cd OligoScout
cmake --preset release
cmake --build --preset release
```

---

## Hizli baslangic

Windows'ta `RUN_*.cmd` sarmalayicilari WSL uzerinden ilgili boru hattini calistirir:

| Script | Islev |
|---|---|
| `RUN_REGRESYON.cmd` | 264 test + 21 adimli uctan uca regresyon |
| `RUN_SEQUENCE_QC.cmd` | Dizi/cift QC (SEQ/CMP/SS/PAIR/AMP; `--profile bisulfite`, `--probes`, `--tails`) |
| `RUN_PANEL.cmd` | Genom-geneli ozgulluk taramasi + termodinamik QC (SPC ailesi) |
| `RUN_VARIANT.cmd` | Varyant katmani (VAR, gnomAD v4) |
| `RUN_BLOCKER.cmd` | CT-05 blocker/prob alel ayrimi |
| `RUN_BISULFITE_INDEX.cmd` | Bisulfit donusturulmus genom indeksi (JOBS ile paralel) |
| `RUN_BISULFITE_PANEL.cmd` | BIS-04 iplik ozgullugu (CT/GA taramasi) |
| `RUN_ISOFORM.cmd` | RT-03 izoform kapsami (GENCODE ekson) |
| `RUN_ANNOTATIONS.cmd` | Psodogen/segdup/ekson anotasyon indeksleri |
| `RUN_DIFF_CPP.cmd` | Python↔C++ diferansiyel testi |

---

## Dogrulama

- **Birim testler:** `reference_py/tests/` — 264 test.
- **Parametre capraz denetim:** `cli/param_crosscheck.py` tablolari primer3'un
  bagimsiz dosyalarina karsi dener (373 karsilastirma, transkripsiyon hatasi 0).
  Bu `verified` durumu VERMEZ; ikinci insan aktarimi hala gerekir.
- **Python↔C++ diferansiyel:** `cli/diff_cpp.py` — 74 dizi × 17 sutun, 0 fark.
- **Uctan uca regresyon:** `cli/regresyon.py` — 21 adim, 20 PASS / 1 SKIP.

---

## Kalibrasyon ve durustluk uyarilari (REP-03)

1. Esiklerin cogu literatur baslangic degerleridir, kalibre EDILMEMISTIR.
2. LNA/modifiye baz katkisi modellenmez; ilgili Tm'ler ALT SINIRDIR.
3. Parametre tablolari "aktarildi, dogrulanmadi" durumundadir (ikinci insan
   aktarimi bekler).
4. A·T bulge kapanis sabiti bir proje modelleme kararidir.
5. Bazi kurallar dis veri/indeks gerektirir (SPC-08 ALT contig'leri C++ FM-index).

Her kanit kaydi kendi model kunyesini ve `param_set_sha256`'yi tasir; bir karar
hangi girdi ve modelden ciktiysa izlenebilir.

---

## Lisans ve atif

Lisans icin depo koku ile `OligoScout/` altindaki ilgili dosyalara bakin.
Termodinamik parametreler orijinal yayinlara aittir (SantaLucia 1998, Allawi &
Peyret, Bommarito 2000, SantaLucia & Hicks 2004, SantaLucia & Peyret 2001,
Sugimoto 1995); atiflar `reference_py/oligoscout_thermo/params.py` kunyelerindedir.
