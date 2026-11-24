# OligoScout MCP sunucusu

OligoScout motorunu (oligoscout.org'un kullandığı `web/local_engine.py` + C++ `oligoscout_panel`)
**Model Context Protocol** üzerinden her sohbet botuna açar. Ek Python paketi gerekmez
(yalnızca standart kütüphane); protokol çekirdeği `protocol.py` içindedir.

| | |
|---|---|
| Protokol | MCP 2025-06-18 (2025-03-26 ve 2024-11-05 ile geriye uyumlu) |
| Taşıma | **stdio** (yerel istemciler) · **Streamable HTTP** `POST /mcp` (uzak istemciler) |
| Ek | REST `POST /api/tools/<ad>` + OpenAPI 3.1 `/openapi.json` (GPT Actions, Gemini function calling) |
| Güvenlik | `Authorization: Bearer` (OLIGOSCOUT_MCP_TOKEN) · Origin denetimi (DNS-rebinding) · 60 istek/60 s sınır |

## Araçlar (15)

| Araç | Ne yapar |
|---|---|
| `oligoscout_status` | Motor sürümü, genom taraması hazır mı, gnomAD önbelleği, anotasyonlar |
| `oligoscout_list_profiles` | qPCR-SYBR, RT-qPCR, ddPCR, ctDNA, Sanger, NGS, bisülfit, standart PCR profilleri |
| `oligoscout_explain_rule` | 91 QC kuralından birinin sade açıklaması, eşikleri ve kanıt düzeyi |
| `oligoscout_resolve_target` | Gen / ekzon → GRCh38 koordinatı (MANE Select, iplik yönüne göre ekzon numarası) |
| `oligoscout_thermo` | Tm, ΔG, tuz modeli yayılımı (LNA → alt sınır) |
| `oligoscout_dimer` | IDT tarzı hizalamalı self/hetero-dimer |
| `oligoscout_design_primers` | Bölge veya yapıştırılan FASTA üzerinde Primer-BLAST tarzı tasarım + genom çapında özgüllük + varyant |
| `oligoscout_check_pair` | Kullanıcının primer çifti için tam QC (90+ kural) |
| `oligoscout_check_specificity` | GRCh38 genelinde bağlanma yerleri ve hedef dışı ürünler |
| `oligoscout_check_variants` | gnomAD v4 varyantları, 3′ uca uzaklık (SNPCheck tarzı) |
| `oligoscout_design_probe` | Hidroliz probu (TaqMan) tasarımı |
| `oligoscout_check_blocker` | Bloklayıcı / clamp oligo denetimi |
| `oligoscout_amplicon_melt` | Amplikon erime eğrisi (uMelt tarzı) |
| `oligoscout_design_multiplex` | 2–8 hedefli multipleks panel tasarımı ve çapraz dimer matrisi |
| `oligoscout_job_status` | Uzun işler için `job_id` sorgusu |

Kaynaklar: `oligoscout://guide` (asistan yorumlama kuralları), `oligoscout://catalog`,
`oligoscout://profiles`, `oligoscout://benchmark` (skor tablosu), `oligoscout://rule/{rule_id}`.
İstem: `design_assay`.

Her sonuç `provenance` (motor sürümü, parametre özeti, referans) ve `not_evaluated` alanı taşır; her
metin yanıtı REP-03 kalibrasyon notuyla biter. Asistanlar `GUIDE.md` kurallarına uymalıdır: sayı
uydurma yok, değerlendirilmeyen katman açıkça söylenir.

## Gereksinimler

* WSL Ubuntu içinde Python ≥ 3.10 (paket gerekmez)
* Genom taraması için `OligoScout/build/release/oligoscout_panel` ve `OligoScout/data/indexes/ppfm`
  (yoksa tasarım yine çalışır; özgüllük `not_evaluated` olarak döner)
* İsteğe bağlı: `reference_py/out/variant_cache.jsonl` (gnomAD çevrimdışı önbellek)

## Çalıştırma

```bash
# stdio (Claude Desktop/Code, Cursor, VS Code, Gemini CLI, LM Studio, Windsurf …)
python3 mcp_server/server.py

# Streamable HTTP + REST + OpenAPI (ChatGPT, claude.ai, uzak istemciler)
OLIGOSCOUT_MCP_TOKEN=uzun-rastgele-dize python3 mcp_server/server.py --http --port 8766
```

Windows'ta `RUN_MCP.cmd` HTTP sunucusunu WSL içinde başlatır.

## İstemci ayarları

Aşağıda `C:\Users\BioComp\Desktop\AYHAN-OligoScout` klasörü WSL'de
`/mnt/c/Users/BioComp/Desktop/AYHAN-OligoScout` olarak görünür.

### Claude Desktop (Windows)
`%APPDATA%\Claude\claude_desktop_config.json`:
```json
{
  "mcpServers": {
    "oligoscout": {
      "command": "wsl.exe",
      "args": ["-d", "Ubuntu", "--exec", "python3",
               "/mnt/c/Users/BioComp/Desktop/AYHAN-OligoScout/mcp_server/server.py"]
    }
  }
}
```

### Claude Code
```bash
claude mcp add oligoscout -- python3 /mnt/c/Users/BioComp/Desktop/AYHAN-OligoScout/mcp_server/server.py
# uzak sunucu:
claude mcp add --transport http oligoscout https://oligoscout.org/mcp --header "Authorization: Bearer <token>"
```

### Cursor (`.cursor/mcp.json`) · Windsurf (`~/.codeium/windsurf/mcp_config.json`) · LM Studio (`mcp.json`)
```json
{ "mcpServers": { "oligoscout": {
    "command": "wsl.exe",
    "args": ["-d", "Ubuntu", "--exec", "python3", "/mnt/c/Users/BioComp/Desktop/AYHAN-OligoScout/mcp_server/server.py"] } } }
```

### VS Code / GitHub Copilot (`.vscode/mcp.json`)
```json
{ "servers": { "oligoscout": { "type": "stdio", "command": "wsl.exe",
    "args": ["-d", "Ubuntu", "--exec", "python3", "/mnt/c/Users/BioComp/Desktop/AYHAN-OligoScout/mcp_server/server.py"] } } }
```

### Gemini CLI (`~/.gemini/settings.json`)
```json
{ "mcpServers": { "oligoscout": {
    "command": "python3", "args": ["/mnt/c/Users/BioComp/Desktop/AYHAN-OligoScout/mcp_server/server.py"] } } }
```
Uzak sunucu için `"httpUrl": "https://oligoscout.org/mcp"` ve `"headers": {"Authorization": "Bearer <token>"}`.

### ChatGPT (geliştirici modu → Connectors) ve claude.ai (Ayarlar → Connectors → Özel bağlayıcı)
Uzak URL: `https://oligoscout.org/mcp` (HTTPS zorunlu). Sunucu `OLIGOSCOUT_OAUTH=1` ile OAuth 2.1
sunar: istemci `401 + WWW-Authenticate: Bearer resource_metadata=...` yanıtından keşif yapar,
kendini dinamik olarak kaydeder (RFC 7591), kullanıcı tek tıklık onay sayfasında izin verir,
PKCE (S256) ile token alınır. Kayıt açıktır (araçlar yalnızca okur, hız sınırı vardır);
`OLIGOSCOUT_OAUTH_INVITE` verilirse onay sayfası davet kodu ister. Uç noktalar:

| | |
|---|---|
| Korunan kaynak (RFC 9728) | `/.well-known/oauth-protected-resource/mcp` |
| Yetkilendirme sunucusu (RFC 8414) | `/.well-known/oauth-authorization-server` |
| Kayıt / yetki / token | `/mcp/oauth/register`, `/mcp/oauth/authorize`, `/mcp/oauth/token` |

### Custom GPT Actions / Gemini function calling / herhangi bir LLM
`https://oligoscout.org/mcp/openapi.json` adresindeki OpenAPI 3.1 şemasını içe aktarın; her araç
`POST /mcp/api/tools/<ad>` olarak çağrılır (yerelde `/api/tools/<ad>`), kimlik doğrulama Bearer token'dır.

## Sunucuya kurulum (oligoscout.org)

```bash
sudo cp deploy/oligoscout-mcp.service /etc/systemd/system/
sudo systemctl daemon-reload && sudo systemctl enable --now oligoscout-mcp
# nginx: deploy/oligoscout-mcp.nginx.conf içindeki location bloklarını mevcut server {} bloğuna ekleyin
sudo nginx -t && sudo systemctl reload nginx
```
Token'ı `/etc/oligoscout/mcp.env` dosyasına koyun (servis bu dosyayı okur), repoya koymayın.

## Doğrulama

* `python3 mcp_server/tests/sdk_conformance.py` (resmî MCP Python SDK gerekir) ve
  `python3 -m unittest discover -s mcp_server/tests` (OAuth akışı, ağ gerekmez).
* Resmî MCP Python SDK 2.3.0 ile stdio ve Streamable HTTP üzerinden **26/26** kontrol geçti:
  initialize (2025-06-18), tools/list (15), Tm motorla birebir (60.1601 °C), hetero-dimer
  `align.best_dimer` ile birebir (uç terimli model), şema ihlali → `isError`, FASTA üzerinde çift QC,
  FASTA üzerinde tasarım (1. çift `design_primers.py` 1. sırasıyla aynı), kaynaklar, kural şablonu, istem.
* Hızlı duman testi: `printf '{"jsonrpc":"2.0","id":1,"method":"tools/list"}\n' | python3 mcp_server/server.py`
* Skor tablosu: `OligoScout/benchmarks/scorecard/SCORECARD.md` (yeniden üretmek için `build_scorecard.py`).

## Arama modları

`oligoscout_design_primers` ve `oligoscout_check_specificity` `search_mode` alır:
`anchored` (varsayılan, tam 12-nt 3′ çapa) ya da `sensitive` (Primer-BLAST kuralı: primer başına
≤ 3 uyumsuzluk, son 5 nt'de ≤ 1; çapa yok). Hassas mod kaba kuvvet genom taramasıyla birebir
doğrulandı. Ürünler F–R, R–F ve tek primer (F–F, R–R) yönelimlerini kapsar. SPC-08, birincil
24 kromozom dışındaki tüm GRCh38.p14 dizilerini (chrM, ALT, yamalar, iskeleler) tarar.

## Docker

```bash
docker build -t oligoscout:0.4.1 .
docker run --rm -p 8766:8766 -v /srv/oligoscout-data:/data -e OLIGOSCOUT_OAUTH=1 \
  -e OLIGOSCOUT_PUBLIC_URL=https://oligoscout.org/mcp oligoscout:0.4.1
```

## Bilinen sınırlar

* Islak laboratuvar kalibrasyonu yok (REP-03). Sonuçlar araştırma amaçlıdır.
* Yalnızca insan GRCh38.p14.
