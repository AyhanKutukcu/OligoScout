#!/usr/bin/env python3
"""
GRCh38.p14 referans normalizasyonu ve dogrulamasi.

Depodaki politikayi birebir uygular
(data/benchmarks/reference/GRCh38.p14.primary24.normalization.txt):

  * kapsam                 : chr1..chr22, chrX, chrY
  * uzunluk                : degistirilmez
  * korunan bazlar         : A C G T N
  * IUPAC belirsizlikleri  : R Y S W K M B D H V  ->  N
  * kucuk harf (soft-mask) : buyuk harfe cevrilir

Her kromozomun normalize edilmis dizisinin SHA-256'si depodaki manifestle
karsilastirilir; uyusmazlik varsa cikis kodu 1 olur.

Kullanim:
  python3 normalize_reference.py --input <ensembl.fa.gz> --output <GRCh38.p14.fa> \
                                 --manifest <manifest.tsv>
  python3 normalize_reference.py --verify-only --output <GRCh38.p14.fa> \
                                 --manifest <manifest.tsv>
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import os
import sys

IUPAC = b"RYSWKMBDHV"
KEEP = b"ACGTN"
LINE_WIDTH = 60

WANTED = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def build_table() -> bytes:
    """256 baytlik ceviri tablosu: buyuk harf + IUPAC -> N."""
    table = bytearray(range(256))
    for i in range(256):
        c = bytes([i]).upper()
        if c in (KEEP + IUPAC):
            table[i] = c[0]
    for c in IUPAC:
        table[c] = ord("N")
        table[ord(chr(c).lower())] = ord("N")
    return bytes(table)


TABLE = build_table()


def ncbi_primary_alias(header: str) -> str | None:
    """src/grch38_reference.cpp icindeki kuralin birebir Python karsiligi.

    Kabul kosullari:
      * kayit adi 'NC_' ile baslar
      * aciklamada 'Homo sapiens chromosome <token>,' gecer
      * aciklamada 'GRCh38.p14 Primary Assembly' gecer
      * token X, Y veya 1..22
    Bu kural `build_grch38_ppfm` araciyla uyumlulugun olcutudur.
    """
    body = header[1:]
    name = body.split()[0] if body.split() else ""
    if not name.startswith("NC_"):
        return None
    marker = "Homo sapiens chromosome "
    begin = body.find(marker)
    if begin == -1:
        return None
    start = begin + len(marker)
    comma = body.find(",", start)
    if comma == -1:
        return None
    if "GRCh38.p14 Primary Assembly" not in body:
        return None
    token = body[start:comma]
    if token in ("X", "Y"):
        return "chr" + token
    if not token or len(token) > 2 or not token.isdigit():
        return None
    n = int(token)
    return f"chr{n}" if 1 <= n <= 22 else None


def canonical_name(header: str) -> str | None:
    """Kaydi chrN'e cozer. Once NCBI kurali, sonra Ensembl/UCSC adlandirmasi."""
    alias = ncbi_primary_alias(header)
    if alias:
        return alias
    name = header[1:].split()[0] if header[1:].split() else ""
    if name.startswith("chr"):
        name = name[3:]
    if name in {str(i) for i in range(1, 23)} or name in {"X", "Y"}:
        return "chr" + name
    return None


def open_maybe_gzip(path: str) -> io.BufferedReader:
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def read_manifest(path: str) -> dict[str, dict[str, int | str]]:
    out: dict[str, dict[str, int | str]] = {}
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            out[row["chromosome"]] = {
                "length": int(row["length"]),
                "sha256": row["normalized_sequence_sha256"],
                "iupac": int(row["iupac_to_n_count"]),
            }
    return out


def iter_records(path: str):
    """(chrN, orijinal baslik, normalize dizi, IUPAC->N sayisi); digerleri atlanir."""
    name: str | None = None
    header: str = ""
    parts: list[bytes] = []
    n_iupac = 0
    with open_maybe_gzip(path) as fh:
        for raw in fh:
            if raw.startswith(b">"):
                if name is not None:
                    yield name, header, b"".join(parts), n_iupac
                parts, n_iupac = [], 0
                header = raw.decode("ascii", "replace").rstrip("\r\n")
                name = canonical_name(header)
                continue
            if name is None:
                continue
            upper = raw.strip().upper()
            n_iupac += sum(upper.count(bytes([c])) for c in IUPAC)
            parts.append(upper.translate(TABLE))
    if name is not None:
        yield name, header, b"".join(parts), n_iupac


def write_fasta(out_path: str, records: dict[str, tuple[str, bytes]]) -> None:
    """Orijinal baslik satiri aynen korunur.

    `build_grch38_ppfm` kromozomu basliktan cozumler (NC_ + 'Homo sapiens
    chromosome N,' + 'GRCh38.p14 Primary Assembly'). Basligi '>chrN' olarak
    yeniden yazmak araci kirar; bu yuzden kaynak basligi degistirilmez.
    """
    tmp = out_path + ".part"
    with open(tmp, "wb") as out:
        for chrom in WANTED:
            header, seq = records[chrom]
            out.write(header.encode("ascii", "replace") + b"\n")
            for i in range(0, len(seq), LINE_WIDTH):
                out.write(seq[i : i + LINE_WIDTH])
                out.write(b"\n")
    os.replace(tmp, out_path)


def report(rows: list[tuple[str, int, str, int]], manifest: dict) -> int:
    """Kabul olcutu uzunluk + SHA-256'dir.

    IUPAC->N sayisi yalnizca bilgilendirme amaclidir: kaynaga gore degisir.
    Ensembl'in DNA dosyalarinda belirsiz bazlar zaten N olarak gelir, bu yuzden
    donusum sayisi 0 cikar; NCBI/UCSC kaynaginda manifestteki sayilar gorulur.
    Her iki durumda da normalize edilmis dizi ayni oldugu icin SHA-256 tutar.
    """
    bad = 0
    print(f"{'kromozom':<9}{'uzunluk':>12}  {'uzunluk':<8}{'sha256':<8}{'iupac->N (bilgi)':<18}")
    for chrom, length, sha, iupac in rows:
        ref = manifest.get(chrom)
        if ref is None:
            print(f"{chrom:<9}{length:>12}  manifestte yok")
            bad += 1
            continue
        len_ok = length == ref["length"]
        sha_ok = sha == ref["sha256"]
        bad += 0 if (len_ok and sha_ok) else 1
        note = f"{iupac}" if iupac == ref["iupac"] else f"{iupac} (manifest: {ref['iupac']})"
        print(
            f"{chrom:<9}{length:>12}  "
            f"{'OK' if len_ok else 'FARKLI':<8}"
            f"{'OK' if sha_ok else 'FARKLI':<8}"
            f"{note:<18}"
        )
    total_iupac = sum(r[3] for r in rows)
    expected_iupac = sum(int(v["iupac"]) for v in manifest.values())
    print(f"\nIUPAC->N donusumu (bilgi): {total_iupac}; manifest kaynaginda: {expected_iupac}")
    if total_iupac != expected_iupac:
        print("  Not: kaynak dosyada belirsiz bazlar zaten N; normalize dizi ayni, "
              "SHA-256 bunu dogruluyor.")
    if bad:
        print(f"\n{bad} kromozom manifestle uyusmuyor (uzunluk/SHA-256).", file=sys.stderr)
    else:
        print("\nTum 24 kromozom manifestle birebir uyusuyor (uzunluk + SHA-256).")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input")
    ap.add_argument("--output", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument(
        "--require-builder-headers",
        action="store_true",
        default=True,
        help="basliklar build_grch38_ppfm ile uyumlu degilse cikis kodu 3",
    )
    ap.add_argument(
        "--allow-any-headers",
        dest="require_builder_headers",
        action="store_false",
        help="baslik uyumsuzlugunda yalniz uyar",
    )
    args = ap.parse_args()

    manifest = read_manifest(args.manifest)

    source = args.output if args.verify_only else args.input
    if not source or not os.path.exists(source):
        print(f"girdi bulunamadi: {source}", file=sys.stderr)
        return 2

    records: dict[str, tuple[str, bytes]] = {}
    rows: list[tuple[str, int, str, int]] = []
    headers: dict[str, str] = {}
    for chrom, header, seq, iupac in iter_records(source):
        rows.append((chrom, len(seq), hashlib.sha256(seq).hexdigest(), iupac))
        headers[chrom] = header
        if not args.verify_only:
            records[chrom] = (header, seq)

    rows.sort(key=lambda r: WANTED.index(r[0]) if r[0] in WANTED else 99)

    missing = [c for c in WANTED if c not in {r[0] for r in rows}]
    if missing:
        print(f"eksik kromozom: {', '.join(missing)}", file=sys.stderr)
        return 2

    bad = report(rows, manifest)

    if not args.verify_only:
        write_fasta(args.output, records)
        size_gb = os.path.getsize(args.output) / 2**30
        print(f"yazildi: {args.output} ({size_gb:.2f} GiB)")

    # --- build_grch38_ppfm baslik uyumlulugu --------------------------------
    incompatible = [c for c in WANTED if ncbi_primary_alias(headers[c]) != c]
    if incompatible:
        print(
            "\nUYARI: basliklar build_grch38_ppfm ile uyumlu degil "
            f"({len(incompatible)}/24 kromozom).\n"
            "  Arac kromozomu basliktan cozumler ve sunlari arar:\n"
            "    kayit adi 'NC_...', aciklamada 'Homo sapiens chromosome N,'\n"
            "    ve 'GRCh38.p14 Primary Assembly'\n"
            f"  Ornek mevcut baslik: {headers[incompatible[0]][:90]}\n"
            "  Cozum: kaynak olarak NCBI RefSeq dosyasini kullanin\n"
            "    GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
            file=sys.stderr,
        )
        if args.require_builder_headers:
            return 3
    else:
        print("basliklar build_grch38_ppfm ile uyumlu (24/24).")

    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
