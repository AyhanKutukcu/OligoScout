# OligoScout — bu bilgisayarda kurulum

Tek giriş noktası: **`RUN_SETUP.cmd`** (çift tıklayın).
Tüm iş WSL2 içinde, tek script olarak sırayla çalışır: `setup/00_bootstrap.sh`.

## Neden WSL

Deponun `README.md` dosyası gereksinimleri açıkça yazıyor: Linux veya WSL2,
C++20 derleyici, CMake 3.25+, git (alt modül desteğiyle), GNU Make.
Yayındaki ölçümler de WSL2/DrvFS üzerinde alınmış. Native Windows derlemesi
desteklenmiyor, bu yüzden her şey WSL içinde yürütülür.

WSL kurulu değilse, yönetici PowerShell'de bir kez:

```powershell
wsl --install -d Ubuntu
```

Kurulumdan sonra Ubuntu'yu bir kez açıp kullanıcı adı/parola belirleyin.

## Klasör düzeni

| Yol | İçerik |
|---|---|
| `OligoScout/` | Depo (GitHub `main`, commit `45161ea`) |
| `setup/00_bootstrap.sh` | Tüm fazları çalıştıran tek script |
| `setup/normalize_reference.py` | Referans normalizasyonu + SHA-256 doğrulaması |
| `downloads/` | Ham indirmeler (NCBI RefSeq FASTA, gerekirse CMake) |
| `logs/` | Her çalıştırmanın logları; `logs/state/` faz işaretleri |
| `OligoScout/data/genomes/GRCh38.p14.fa` | Normalize edilmiş referans (~3,1 GB) |
| `OligoScout/data/indexes/ppfm/` | 24 kromozom PPFM indeksi (~8,5 GB) |

## Fazlar

| # | İş | Yaklaşık süre |
|---|---|---|
| 0 | Ortam raporu (WSL, CPU, RAM, disk) | saniyeler |
| 1 | `apt` ile derleme bağımlılıkları | 1–5 dk (sudo parolası sorar) |
| 2 | CMake ≥ 3.25 garantisi (gerekirse yerel kurulum) | 0–2 dk |
| 3 | Depo + `primer3` alt modülü, pinlenmiş commit'ler | 1–2 dk |
| 4 | `primer3` kütüphanesi (`make`) | 1–3 dk |
| 5 | OligoScout derlemesi (Ninja, `release` preset) | 10–30 dk |
| 6 | Hızlı regresyon testleri (genom gerektirmeyenler) | 2–10 dk |
| 7 | GRCh38.p14: indirme + normalizasyon + SHA-256 doğrulama | 10–25 dk |
| 8 | 24 kromozom PPFM indeksi (paralel) + manifest | 30–120 dk |
| 9 | Tam test paketi (78 test) + özet | 10–30 dk |

Her faz bittiğinde `logs/state/faz_N.done` yazılır; script yeniden
çalıştırıldığında tamamlanmış fazları atlar. Yeniden çalıştırmak için:

```
RUN_SETUP.cmd --only 8      # yalnız indeks fazı
RUN_SETUP.cmd --from 5      # 5. fazdan itibaren
RUN_SETUP.cmd --force       # tamamlananları da yeniden çalıştır
```

Ayarlar ortam değişkeniyle ezilebilir (WSL kabuğundan):

```bash
PRESET=release-avx2 INDEX_PAR=8 SA_RATE=8 bash setup/00_bootstrap.sh --only 8
```

## Referansın doğrulanması

Depo, 24 birincil kromozomun **normalize edilmiş dizi SHA-256**'larını
`data/benchmarks/reference/GRCh38.p14.primary24.manifest.tsv` içinde taşıyor.
Normalizasyon politikası aynı klasördeki `...normalization.txt` dosyasında:
uzunluk korunur, A/C/G/T/N korunur, IUPAC belirsizlikleri (R Y S W K M B D H V)
N'e çevrilir; manifestte toplam 99 dönüşüm bekleniyor.

`setup/normalize_reference.py` kaynak olarak **NCBI RefSeq**
`GCF_000001405.40_GRCh38.p14_genomic.fna.gz` dosyasını alır, bu politikayı
uygular ve her kromozomun uzunluğunu, SHA-256'sını ve IUPAC→N sayısını
manifestle karşılaştırır. Tek bir kromozom bile tutmazsa çıkış kodu 1 olur ve
hangi kromozomun farklı olduğu listelenir. Bu, QC kataloğundaki SEQ-06
(referans kimliği beyanı) maddesinin uygulanmış halidir.

**Kaynak neden NCBI olmak zorunda.** `build_grch38_ppfm` kromozomu FASTA
başlığından çözümler. `src/grch38_reference.cpp` bir kaydı ancak kayıt adı
`NC_` ile başlıyorsa, açıklamada `Homo sapiens chromosome <N>,` ve
`GRCh38.p14 Primary Assembly` geçiyorsa kabul eder. Ensembl veya UCSC
başlıkları (`>1`, `>chr1`) diziler doğru olsa bile aracı kırar ve
`ERROR: Requested GRCh38 primary chromosome was not found exactly once`
hatasını verir. Bu yüzden normalizasyon **kaynak başlığını aynen korur** ve
yazdıktan sonra 24 kaydın hepsinin bu kurala uyduğunu denetler; uymazsa çıkış
kodu 3 ile durur ve Faz 7 referansı NCBI kaynağından yeniden üretir.

## Performans notları

- **Disk:** Proje klasörü Windows NTFS üzerinde (DrvFS). Derleme ve indeks
  yazımı ext4'e göre belirgin yavaştır. Çok yavaş gelirse `REPO` ve `INDEXDIR`
  değişkenlerini WSL ev dizinine (`~/oligoscout`) yönlendirip oradan çalışın.
- **RAM:** WSL2 öntanımlı olarak ana belleğin yarısını alır. Bu makinede
  126 GiB'ın yarısı ≈ 63 GiB, indeks derlemesi için fazlasıyla yeterli.
  Sınırlamak veya artırmak için `%USERPROFILE%\.wslconfig`:

  ```ini
  [wsl2]
  memory=64GB
  processors=24
  ```

  Ardından `wsl --shutdown`.
- **Paralellik:** `INDEX_PAR` eş zamanlı kromozom indeksi sayısıdır
  (öntanımlı 6). Kromozom indeksleri birbirinden bağımsızdır; paralel üretim
  sonucu değiştirmez, yalnız süreyi kısaltır. Her indeks tamamlandığında
  builder'ın kendi `ALL_CHECKS YES` kontrolü aranır, geçmezse `.part` dosyası
  silinir.

## Sorun giderme

- `cmake --preset release` "preset bulunamadı" derse: CMake 3.25'ten eski.
  Faz 2 bunu yakalar; yine de olursa `source setup/cmake_path.env`.
- `make -C third_party/primer3/src` boş klasör hatası verirse alt modül
  çekilmemiştir: `RUN_SETUP.cmd --only 3`.
- İndirme yarıda kalırsa `curl -C -` ile kaldığı yerden devam eder;
  `RUN_SETUP.cmd --only 7` yeterlidir.
- `ERROR: Requested GRCh38 primary chromosome was not found exactly once`:
  FASTA başlıkları NCBI biçiminde değil. `RUN_SETUP.cmd --only 7 --force`
  referansı NCBI kaynağından yeniden üretir.
- `ppfm_manifest` testi `Unable to open PPFM manifest` diyorsa indeksler henüz
  üretilmemiştir. Bu test `results/benchmarks/grch38_ppfm_batch/manifest.sa8.tsv`
  dosyasını depo kökünden okur ve başlık satırının tam olarak
  `chromosome⇥file_bytes⇥sha256⇥status` olmasını ister; Faz 8 manifesti oraya
  bu biçimde yazar. Testi Faz 8'den önce çalıştırmak anlamsızdır, bu yüzden
  hızlı test fazından çıkarılmıştır.
- Bir faz tamamlandı işaretliyse ama yeniden çalıştırmak istiyorsanız
  `--only N --force` kullanın; `--from N` tamamlanmış fazları yine atlar.
