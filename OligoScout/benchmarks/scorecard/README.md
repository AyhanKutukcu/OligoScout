# Benchmark scorecard

Bu klasör OligoScout'un tüm doğrulama sonuçlarını tek bir skor tablosunda toplar
(`SCORECARD.md`, MCP'de `oligoscout://benchmark` olarak da sunulur).

| Betik | Ne ölçer | Süre (WSL, 2 çekirdek) |
|---|---|---|
| `thermo_vs_primer3.py` | 2.048 Test #81 primeri: Tm (oligotm) ve dimer/hairpin ΔG (ntthal) | ~15 s |
| `exhaustive_scan.cpp` | FM-indeksten bağımsız, çapasız, ≤ 3 uyumsuzluklu kaba kuvvet genom taraması | ~2 s / primer |
| `blindspot_analysis.py` | Çapa kör noktası, Primer-BLAST kuralı (3′ son 5 nt'de ≤ 1 uyumsuzluk), tek primer ürünleri | saniyeler |
| `build_scorecard.py` | Hepsini `mcp_server/benchmark_scores.json` + `SCORECARD.md` içine toplar | <1 s |

## Yeniden üretim (WSL Ubuntu)

```bash
# primer3 (third_party içinden, sistem glibc'sine göre)
mkdir -p ~/p3 && cp -r OligoScout/third_party/primer3/src ~/p3/ && cd ~/p3/src && make clean && make oligotm ntthal
cd -  # proje kökü

cd OligoScout/benchmarks/scorecard
python3 thermo_vs_primer3.py --primer3 ~/p3/src

g++ -std=c++17 -O3 -march=native -pthread exhaustive_scan.cpp -o ~/exhaustive_scan
# primers.tsv: id <TAB> dizi   (id olarak dizinin kendisi kullanılırsa analiz doğrudan eşler)
~/exhaustive_scan ../../data/genomes/GRCh38.p14.fa primers.tsv hits.tsv 3 2
python3 blindspot_analysis.py --hits hits.tsv --panel panel.tsv \
    --oligoscout-bindings <panel.bindings.tsv> --oligoscout-products <panel.products.tsv> --out results/<ad>

python3 build_scorecard.py
```

## Bağımsız doğrulamalar (bu klasörde üretildi)

* **ctDNA paneli (20 çift):** kaba kuvvet taramasının çapalı alt kümesi `oligoscout_panel` ile
  **site düzeyinde birebir** (4.931 / 4.931) ve **ürün düzeyinde birebir** (26 / 26).
* **Test #81 ilk 64 dengeli çift:** çapalı alt küme Test #81 oracle'ının çift başına bağlanma
  (12.141) ve ürün (1.683) sayılarıyla 64/64 çiftte aynı.

## Primer-BLAST ile karşılaştırma (kullanıcı tarafından çalıştırılacak)

NCBI'ye bu oturumdan erişilemediği için Primer-BLAST web karşılaştırması yapılmadı. Önerilen protokol:

1. `results/blindspot_ctdna/per_pair.tsv` içindeki çiftleri Primer-BLAST'a *"Primer pair specificity
   checking"* modunda girin: veritabanı **Genomes for selected organisms → Homo sapiens GRCh38**,
   diğer ayarlar varsayılan (en az 2 uyumsuzluk, bunların ≥ 2'si son 5 bazda; 6+ uyumsuzluk yok sayılır;
   ürün ≤ 4000 bp).
2. Her çift için Primer-BLAST'ın listelediği hedef dışı ürünlerin sayısını ve konumlarını kaydedin.
3. Karşılaştırın: `oligoscout` sütunu (OligoScout sözleşmesi) ↔ `primer_blast_rule` (aynı kuralın
   kaba kuvvet karşılığı) ↔ Primer-BLAST'ın gerçek çıktısı. Primer-BLAST'ın BLAST tohumlarından
   kaynaklanan kaçırmaları üçüncü sütunla ikinci arasındaki fark gösterir.
