#!/usr/bin/env python3
"""Parametre aktarim harness'i: sablon uret, denetle, koda cevir.

Amac: kapali erisim makalelerdeki tablolari elle girerken hatayi *girerken*
yakalamak. Hicbir deger burada uretilmez; arac yalnizca sizin girdiginiz
degerleri denetler.

Alt komutlar:
  template <tablo> <cikti.tsv>     bos sablon: gereken tum anahtarlar yazili
  check <tablo> <dosya.tsv> [ikinci.tsv]
                                   butunluk, ters-tamamlayici tutarliligi,
                                   dG37 = dH - 310.15*dS/1000 denetimi ve
                                   (iki dosya verilirse) cift aktarim diff'i
  emit <tablo> <dosya.tsv>         params.py'ye yapistirilacak Python blogu
                                   + icerik SHA-256

Tablolar:
  mismatch   ic tekli uyumsuzluk NN yiginlari (Allawi/SantaLucia, Peyret)
  dangling   sarkan uclar (Bommarito 2000)
  loop       dongu enerjileri (hairpin/bulge/ic dongu/multiloop)

Denetimler:
  1. Butunluk   : gereken anahtar kumesi programatik uretilir; eksik ya da
                  fazla anahtar hata verir.
  2. RC simetri : ayni yigin diger iplikten okundugunda ayni deger olmali.
                  Tek bir rakam hatasini yakalayan en guclu denetim budur.
  3. dG tutarlilik: makalede verilen dG37 sutunu dH ve dS'ten yeniden
                  hesaplananla tabloya ait yuvarlama payi icinde olmali.
  4. Aralik     : fiziksel olarak makul olmayan degerler isaretlenir
                  (reddedilmez, uyarilir).
  5. Cift aktarim: iki bagimsiz aktarimin birebir ayni olmasi.

Kullanim ornegi:
  python3 cli/param_transcribe.py template mismatch mismatch_bos.tsv
  # iki kisi bagimsiz doldurur
  python3 cli/param_transcribe.py check mismatch ayhan.tsv ikinci.tsv
  python3 cli/param_transcribe.py emit mismatch ayhan.tsv > blok.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys

BASES = "ACGT"
WC = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}

DG_TEMPERATURE = 310.15
DG_TOLERANCE = 0.05          # kcal/mol, yuvarlama payi
# Kaynaklar dH/dS'yi bir, dG37'yi iki ondalik basamakla verir.
DG_ROUNDED_TOLERANCE = 0.05 + DG_TEMPERATURE * 0.05 / 1000.0 + 0.005
DH_RANGE = (-20.0, 10.0)     # kcal/mol
DS_RANGE = (-60.0, 30.0)     # cal/(mol K)


# ---------------------------------------------------------------------------
# Anahtar uzaylari
# ---------------------------------------------------------------------------
def mismatch_keys() -> list[str]:
    """'XY/ZW' : ust 5'-XY-3', alt 3'-ZW-5'; tam olarak bir konum uyumsuz.

    Ters-tamamlayici esdegeri 'WZ/YX'; her sinifin temsilcisi alfabetik
    olarak kucuk olanidir.
    """
    seen: dict[str, str] = {}
    for x in BASES:
        for y in BASES:
            for z in BASES:
                for w in BASES:
                    p1 = (x, z) in WC
                    p2 = (y, w) in WC
                    if p1 == p2:           # ikisi de WC ya da ikisi de uyumsuz
                        continue
                    key = f"{x}{y}/{z}{w}"
                    rc = f"{w}{z}/{y}{x}"
                    rep = min(key, rc)
                    seen[rep] = rep
    return sorted(seen)


def rc_key_mismatch(key: str) -> str:
    top, bottom = key.split("/")
    x, y = top
    z, w = bottom
    return f"{w}{z}/{y}{x}"


def dangling_keys() -> list[str]:
    """'D5:<sarkan>:<ust>-<alt>' ve 'D3:...' — 2 x 4 x 4 = 32 anahtar."""
    keys = []
    for end in ("D5", "D3"):
        for dangle in BASES:
            for top in BASES:
                keys.append(f"{end}:{dangle}:{top}-{COMP[top]}")
    return sorted(keys)


def loop_keys() -> list[str]:
    sizes = (*range(3, 11), 12, 14, 16, 18, 20, 25, 30)
    return [f"{kind}:{size}" for kind in ("HAIRPIN", "BULGE", "INTERNAL")
            for size in ((1, 2, *sizes) if kind == "BULGE" else sizes)]


TABLES = {
    "mismatch": {
        "keys": mismatch_keys,
        "dg_tolerance": DG_ROUNDED_TOLERANCE,
        "rc": rc_key_mismatch,
        "columns": ["key", "dH_kcal", "dS_cal", "dG37_kcal_paper"],
        "citation": ("Allawi & SantaLucia 1997/1998; Peyret et al. 1999 "
                     "(G.T, G.A, C.T, A.C, A.A, C.C, G.G, T.T)"),
        "note": "ust iplik 5'->3', alt iplik 3'->5'",
    },
    "dangling": {
        "keys": dangling_keys,
        "dg_tolerance": DG_ROUNDED_TOLERANCE,
        "rc": None,
        "columns": ["key", "dH_kcal", "dS_cal", "dG37_kcal_paper"],
        "citation": "Bommarito, Peyret & SantaLucia, NAR 2000;28:1929",
        "note": "D5 = 5' sarkan uc, D3 = 3' sarkan uc; <sarkan>:<ust>-<alt>",
    },
    "loop": {
        "keys": loop_keys,
        "rc": None,
        "columns": ["key", "dG37_kcal_paper", "dH_kcal"],
        "citation": "SantaLucia & Hicks 2004, Tablo 4 (DNA, 1 M NaCl, 37 C)",
        "note": "Yalniz yayinlanmis uzunluk terimleri. dH = 0 Tablo 4 varsayimi.",
        "required_value": "dG37_kcal_paper",
    },
}


# ---------------------------------------------------------------------------
# Okuma / yazma
# ---------------------------------------------------------------------------
def write_template(table: str, path: str) -> None:
    spec = TABLES[table]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# tablo: {table}\n")
        fh.write(f"# kaynak: {spec['citation']}\n")
        fh.write(f"# not: {spec['note']}\n")
        fh.write("# bos birakilan satir eksik sayilir; birimler: "
                 "dH kcal/mol, dS cal/(mol K), dG37 kcal/mol\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(spec["columns"])
        for key in spec["keys"]():
            w.writerow([key] + [""] * (len(spec["columns"]) - 1))
    n = len(spec["keys"]())
    print(f"{path}: {n} anahtarlik bos sablon yazildi ({table})")


def read_table(table: str, path: str) -> dict[str, dict[str, float]]:
    spec = TABLES[table]
    rows: dict[str, dict[str, float]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    for r in reader:
        key = (r.get("key") or "").strip()
        if not key:
            continue
        values: dict[str, float] = {}
        for col in spec["columns"][1:]:
            raw = (r.get(col) or "").strip()
            if raw == "":
                continue
            try:
                value = float(raw)
                if not math.isfinite(value):
                    raise ValueError("sonlu olmayan deger")
                values[col] = value
            except ValueError:
                raise SystemExit(f"{path}: '{key}' satirinda sayiya "
                                 f"cevrilemeyen deger: {col}={raw!r}")
        rows[key] = values
    return rows


# ---------------------------------------------------------------------------
# Denetimler
# ---------------------------------------------------------------------------
def check(table: str, path: str, second: str | None) -> int:
    spec = TABLES[table]
    required = set(spec["keys"]())
    rows = read_table(table, path)
    problems: list[str] = []
    warnings: list[str] = []

    rc_fn = spec["rc"]

    # Ters-tamamlayici esdegeri satirlar istege baglidir: sablon her sinifin
    # tek temsilcisini ister, ama makalede iki yon de listeleniyorsa ikisini
    # de girebilirsiniz; o zaman RC simetri denetimi gercek bir hata
    # yakalayiciya donusur.
    accepted = set(required)
    if rc_fn:
        accepted |= {rc_fn(k) for k in required}

    missing = sorted(required - set(rows))
    extra = sorted(set(rows) - accepted)
    required_value = spec.get("required_value")
    empty = sorted(k for k, v in rows.items()
                   if not v or (required_value and required_value not in v))

    if missing:
        problems.append(f"eksik anahtar ({len(missing)}): "
                        f"{', '.join(missing[:8])}"
                        f"{' ...' if len(missing) > 8 else ''}")
    if extra:
        problems.append(f"taninmayan anahtar ({len(extra)}): "
                        f"{', '.join(extra[:8])}")
    if empty:
        problems.append(f"doldurulmamis satir ({len(empty)}): "
                        f"{', '.join(empty[:8])}"
                        f"{' ...' if len(empty) > 8 else ''}")

    filled = {k: v for k, v in rows.items() if v and k in accepted}

    # dG tutarliligi
    dg_tolerance = spec.get("dg_tolerance", DG_TOLERANCE)
    dg_bad = []
    for k, v in filled.items():
        if "dH_kcal" in v and "dS_cal" in v and "dG37_kcal_paper" in v:
            calc = v["dH_kcal"] - DG_TEMPERATURE * v["dS_cal"] / 1000.0
            if abs(calc - v["dG37_kcal_paper"]) > dg_tolerance:
                dg_bad.append(
                    f"{k}: makale dG37={v['dG37_kcal_paper']:.2f}, "
                    f"dH/dS'ten {calc:.2f}"
                )
    if dg_bad:
        problems.append(f"dG37 tutarsizligi ({len(dg_bad)}): " +
                        "; ".join(dg_bad[:6]))

    # RC simetri
    if rc_fn:
        rc_bad = []
        for k, v in filled.items():
            partner = rc_fn(k)
            if partner == k or partner not in filled:
                continue
            for col, val in v.items():
                other = filled[partner].get(col)
                if other is not None and abs(other - val) > 1e-9:
                    rc_bad.append(f"{k} vs {partner}: {col} {val} != {other}")
        if rc_bad:
            problems.append(f"ters-tamamlayici tutarsizligi ({len(rc_bad)}): " +
                            "; ".join(rc_bad[:6]))

    # aralik
    for k, v in filled.items():
        if "dH_kcal" in v and not (DH_RANGE[0] <= v["dH_kcal"] <= DH_RANGE[1]):
            warnings.append(f"{k}: dH={v['dH_kcal']} beklenen araligin disinda")
        if "dS_cal" in v and not (DS_RANGE[0] <= v["dS_cal"] <= DS_RANGE[1]):
            warnings.append(f"{k}: dS={v['dS_cal']} beklenen araligin disinda")

    # cift aktarim
    if second:
        rows2 = read_table(table, second)
        diff = []
        for k in sorted(set(rows) | set(rows2)):
            a, b = rows.get(k, {}), rows2.get(k, {})
            if a != b:
                diff.append(f"{k}: {a} != {b}")
        if diff:
            problems.append(f"cift aktarim farki ({len(diff)}): " +
                            "; ".join(diff[:6]))
        else:
            print(f"cift aktarim: {len(rows)} satirda fark yok")

    print(f"tablo          : {table}")
    print(f"gereken anahtar: {len(required)}")
    print(f"dolu satir     : {len(filled)}")
    for w in warnings[:10]:
        print(f"  UYARI {w}")
    if problems:
        print("\nSORUNLAR:")
        for p in problems:
            print(f"  - {p}")
        print("\nTablo HENUZ kullanilabilir degil.")
        return 1
    print("\nTum denetimler gecti. Durum 'transcribed_unverified' yapilabilir; "
          "ikinci bagimsiz aktarim da eslesiyorsa 'verified'.")
    return 0


def emit(table: str, path: str) -> int:
    spec = TABLES[table]
    rows = read_table(table, path)
    required = spec["keys"]()
    data: dict[str, list[float]] = {}
    for k in required:
        v = rows.get(k, {})
        if "dH_kcal" in v and "dS_cal" in v:
            data[k] = [v["dH_kcal"], v["dS_cal"]]
        elif "dG37_kcal_paper" in v:
            data[k] = [v["dG37_kcal_paper"]]
    blob = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    digest = hashlib.sha256(blob).hexdigest()
    print(f"# {table}: {len(data)}/{len(required)} anahtar")
    print(f"# kaynak: {spec['citation']}")
    print(f"# icerik sha256: {digest}")
    print("DATA = {")
    for k in sorted(data):
        print(f"    {k!r}: {data[k]!r},")
    print("}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["template", "check", "emit"])
    ap.add_argument("table", choices=sorted(TABLES))
    ap.add_argument("path")
    ap.add_argument("second", nargs="?")
    args = ap.parse_args()

    if args.command == "template":
        write_template(args.table, args.path)
        return 0
    if args.command == "check":
        return check(args.table, args.path, args.second)
    return emit(args.table, args.path)


if __name__ == "__main__":
    sys.exit(main())
