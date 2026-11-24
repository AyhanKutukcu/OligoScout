"""Read indexed slices of the normalized local GRCh38.p14 FASTA."""
from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path

CHROMOSOMES = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX", "chrY")
REGION = re.compile(r"(?:chr)?([1-9]|1[0-9]|2[0-2]|X|Y):([0-9]+)-([0-9]+)", re.I)
HEADER_CHROM = re.compile(rb"\bchromosome (1[0-9]|2[0-2]|[1-9]|X|Y),")
MAX_REGION = 5000


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_region(value):
    if not isinstance(value, str):
        raise ValueError("GRCh38 bölgesi gerekli.")
    match = REGION.fullmatch(value.strip().replace(".", "").replace(",", ""))
    if not match:
        raise ValueError("Bölgeyi chr22:20000550-20000870 biçiminde girin.")
    chrom, first, last = match.groups()
    start, end = int(first) - 1, int(last)
    if start < 0 or not 1 <= end - start <= MAX_REGION:
        raise ValueError(f"Bölge 1 tabanlı olmalı ve en fazla {MAX_REGION} baz içermeli.")
    return f"chr{chrom.upper()}", start, end


class Reference:
    """Use the normalization manifest to jump over chromosomes, then validate bounds."""

    def __init__(self, fasta, manifest):
        self.fasta = Path(fasta)
        rows = read_tsv(manifest)
        lengths = {r["chromosome"]: int(r["length"]) for r in rows}
        self.entries = {}
        with self.fasta.open("rb") as handle:
            header = handle.readline()
            while header:
                chrom_match = HEADER_CHROM.search(header)
                if not header.startswith(b">") or not chrom_match:
                    raise ValueError("GRCh38 FASTA ve normalizasyon manifesti eşleşmiyor.")
                chrom = "chr" + chrom_match.group(1).decode("ascii")
                if chrom not in lengths or chrom in self.entries:
                    raise ValueError(f"FASTA başlığı beklenen kromozomla eşleşmiyor: {chrom}")
                offset = handle.tell()
                line = handle.readline()
                bases = len(line.rstrip(b"\r\n"))
                if bases == 0 or line.startswith(b">"):
                    raise ValueError(f"FASTA kromozomu boş: {chrom}")
                length = lengths[chrom]
                self.entries[chrom] = (offset, length, bases, len(line))
                last = length - 1
                last_base_offset = offset + (last // bases) * len(line) + last % bases
                handle.seek(last_base_offset)
                if handle.read(1) not in (b"A", b"C", b"G", b"T", b"N"):
                    raise ValueError(f"FASTA satır biçimi veya kromozom uzunluğu uyuşmuyor: {chrom}")
                if handle.readline().strip():
                    raise ValueError(f"FASTA kromozomu beklenen uzunluktan büyük: {chrom}")
                header = handle.readline()
        if set(self.entries) != set(CHROMOSOMES):
            raise ValueError("GRCh38 referansında beklenen 24 kromozom bulunamadı.")

    def read(self, chrom, start, end):
        offset, length, bases, width = self.entries[chrom]
        if not 0 <= start < end <= length:
            raise ValueError(f"{chrom} koordinatları 1–{length:,} aralığında olmalı.")
        byte_start = offset + start // bases * width + start % bases
        last = end - 1
        byte_end = offset + last // bases * width + last % bases + 1
        with self.fasta.open("rb") as handle:
            handle.seek(byte_start)
            sequence = handle.read(byte_end - byte_start).replace(b"\n", b"").replace(b"\r", b"")
        if len(sequence) != end - start or not re.fullmatch(rb"[ACGTN]+", sequence):
            raise ValueError("Referans bölgesi okunamadı. FASTA biçimini kontrol edin.")
        return sequence.decode("ascii")


def inspect_shards(index_dir):
    manifest = index_dir / "manifest.sa8.tsv"
    if not manifest.is_file():
        return {"ready": False, "shards": 0}
    try:
        rows = read_tsv(manifest)
        if sorted(r["chromosome"] for r in rows) != sorted(CHROMOSOMES):
            return {"ready": False, "shards": len(rows)}
        for row in rows:
            path = index_dir / f"{row['chromosome']}.sa8.ppfm"
            if path.stat().st_size != int(row["file_bytes"]) or not re.fullmatch(
                r"[a-f0-9]{64}", row["sha256"]
            ):
                return {"ready": False, "shards": len(rows)}
    except (OSError, KeyError, ValueError):
        return {"ready": False, "shards": 0}
    return {
        "ready": True, "shards": len(rows),
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
    }
