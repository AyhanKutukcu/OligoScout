"""Izoform kapsami (RT-03): bir primer ciftinin hangi transkript izoformlarini
hedefledigini transkript-ekson anotasyonundan belirler.

Transkriptom FM-index GEREKMEZ: GENCODE GTF'den cikarilan ekson araliklari
(cli/build_annotations.py -> transcript_exons.tsv) ve primerlerin GENOM
koordinatlari kullanilir. Bir transkript, her iki primer footprint'i de kendi
EKSONLARI icinde barindiriyorsa (intron cikarilmis) o cift tarafindan
HEDEFLENIR. Hedeflenen listesi ile genin tum izoformlari karsilastirilir;
hedeflenmeyenler rapor edilir (RT-03: hedeflenen/hedeflenmeyen transkript listesi).
"""

from __future__ import annotations

from dataclasses import dataclass, field


def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Ortusmesiz birlesmis araliklar (yari-acik, 0-tabanli)."""
    out: list[tuple[int, int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def covered_by_exons(iv: tuple[int, int], exons: list[tuple[int, int]]) -> bool:
    """`iv` (start,end) tam olarak eksonlarin BIRLESIMI icinde mi (her baz bir
    eksonda)? Ekson-ekson kavsagini asan primerler de kapsanir."""
    s, e = iv
    for es, ee in _merge(exons):
        if s >= es and s < ee:
            # bu eksondan basliyor; bitis bu ekson zincirinde mi?
            if e <= ee:
                return True
            # kavsak: kalan kismi sonraki bitisik ekson(lar) kapsamali
            return covered_by_exons((ee, e), exons) if ee > s else False
    return False


@dataclass
class TranscriptIndex:
    """gene_name -> {transcript_id: (chrom, [ (exon_start, exon_end), ... ])}."""
    by_gene: dict = field(default_factory=dict)

    def genes(self) -> int:
        return len(self.by_gene)

    def transcripts(self, gene: str) -> dict:
        return self.by_gene.get(gene, {})


def load_transcript_exons(path: str) -> TranscriptIndex:
    """transcript_exons.tsv (chrom,start,end,gene_name,transcript_id) yukler."""
    idx = TranscriptIndex()
    import csv
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["gene_name"]
            tid = row["transcript_id"]
            g = idx.by_gene.setdefault(gene, {})
            chrom, exons = g.get(tid, (row["chrom"], []))
            exons.append((int(row["start"]), int(row["end"])))
            g[tid] = (chrom, exons)
    return idx


def _rec(cid, value, threshold, severity, note=""):
    return {"id": cid, "value": value, "threshold": threshold,
            "severity": severity, "note": note}


def check_rt03(gene: str, chrom: str,
               fwd_iv: tuple[int, int], rev_iv: tuple[int, int],
               index: TranscriptIndex) -> dict:
    """RT-03: izoform kapsami. `fwd_iv`/`rev_iv` primerlerin GENOM koordinatlari
    (0-tabanli, yari-acik). Gen izoformlari arasindan her iki primeri de
    eksonlarinda barindiran transkriptler HEDEFLENIR."""
    transcripts = index.transcripts(gene)
    if not transcripts:
        return _rec("RT-03",
                    {"gene": gene, "status": "no_annotation",
                     "targeted": [], "untargeted": [], "n_isoforms": 0},
                    {}, "INFO",
                    "gen transkript anotasyonunda yok (transcript_exons.tsv)")
    targeted, untargeted = [], []
    for tid, (tchrom, exons) in transcripts.items():
        if tchrom != chrom:
            untargeted.append(tid)
            continue
        if covered_by_exons(fwd_iv, exons) and covered_by_exons(rev_iv, exons):
            targeted.append(tid)
        else:
            untargeted.append(tid)
    n = len(transcripts)
    sev = "INFO"
    if not targeted:
        sev = "WARN"             # hicbir izoformu hedeflemiyor
    return _rec("RT-03",
                {"gene": gene, "status": "evaluated",
                 "n_isoforms": n, "n_targeted": len(targeted),
                 "targeted": sorted(targeted), "untargeted": sorted(untargeted)},
                {"warn": "hicbir izoform hedeflenmedi"},
                sev,
                "izoform kapsami GENCODE ekson anotasyonundan; "
                "transkriptom FM-index gerekmez")
