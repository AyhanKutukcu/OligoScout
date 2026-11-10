"""Multipleks capraz reaktivite: MPX-01, MPX-02, MPX-03.

Makalenin 3. sinirliligi: bir panelde yalniz tasarlanmis ciftler degil,
**tum primer kombinasyonlari** urun verebilir. F x F, R x R ve ciftler arasi
kombinasyonlar dahil.

Onemli: bunun icin yeni bir genom taramasi gerekmez. `oligoscout_panel`
zaten her primerin tum genomdaki tum baglanma yerlerini yaziyor. Capraz
urunler bu kayitlardan kurulur:

    ayni kromozomda, biri + biri - iplikte, birbirine bakan (inward) ve
    urun uzunlugu penceresinde olan her baglanma yeri cifti bir urundur.

Bu, aramanin eksiksizlik garantisini korur: tarama zaten eksiksizdi, burada
yalnizca ayni kayit kumesi uzerinde farkli bir eslestirme yapiliyor.
"""

from __future__ import annotations

import collections
import itertools
from dataclasses import dataclass

from . import align, concentration, salt


@dataclass(frozen=True)
class BindingSite:
    primer_id: str
    primer: str
    chrom: str
    start: int
    end: int
    strand: str
    mismatches: int


@dataclass(frozen=True)
class CrossProduct:
    left_id: str
    right_id: str
    chrom: str
    start: int
    end: int
    length: int
    mismatches: int
    designed_pair: bool          # tasarlanmis cift mi, capraz mi

    def as_row(self) -> list:
        return [self.left_id, self.right_id, self.chrom, self.start, self.end,
                self.length, self.mismatches,
                "designed" if self.designed_pair else "cross"]


def load_bindings(path: str) -> tuple[list[BindingSite], dict[str, set[str]]]:
    """oligoscout_panel bindings TSV'sini okur.

    KIMLIK PRIMER DIZISIDIR, pair_id degil. Ayni primer birden cok cift
    tarafindan paylasilabilir; kimligi pair_id'ye baglamak paylasilan
    primerlerin tasarlanmis eslesmesini kaybettirir ve o urunleri yanlislikla
    "capraz" gosterir.

    (siteler, dizi -> onu kullanan pair_id kumesi) dondurur.
    """
    import csv

    seen: set[tuple] = set()
    sites: list[BindingSite] = []
    owners: dict[str, set[str]] = collections.defaultdict(set)
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            owners[r["primer"]].add(f'{r["pair_id"]}_{r["role"]}')
            key = (r["primer"], r["chrom"], r["start"], r["strand"])
            if key in seen:
                continue
            seen.add(key)
            sites.append(BindingSite(
                primer_id=r["primer"],          # kimlik = dizi
                primer=r["primer"],
                chrom=r["chrom"],
                start=int(r["start"]),
                end=int(r["end_exclusive"]),
                strand=r["strand"],
                mismatches=int(r["mismatches"]),
            ))
    return sites, dict(owners)


def iter_products(
    sites: list[BindingSite],
    designed: set[frozenset[str]],
    min_length: int = 50,
    max_length: int = 3000,
) -> "collections.abc.Iterator[CrossProduct]":
    """Tum primer kombinasyonlari icin urunleri uretir (MPX-02).

    Kural `assemble_primer_pair_hits` ile ayni: ayni kromozom, ileri site
    (+ iplik) solda, geri site (- iplik) sagda, amplikon uzunlugu pencerede.
    Uretec olmasinin nedeni: sayim her zaman EKSIKSIZ yapilmali, kirpma
    yalnizca dosyaya yazilan satirlara uygulanmali.
    """
    by_chrom: dict[str, list[BindingSite]] = collections.defaultdict(list)
    for s in sites:
        by_chrom[s.chrom].append(s)

    for chrom, group in by_chrom.items():
        plus = sorted((s for s in group if s.strand == "+"), key=lambda s: s.start)
        minus = sorted((s for s in group if s.strand == "-"), key=lambda s: s.start)
        if not plus or not minus:
            continue

        j0 = 0
        for f in plus:
            # ilk olasi geri site: amplikon >= min_length
            while j0 < len(minus) and minus[j0].end < f.start + min_length:
                j0 += 1
            j = j0
            while j < len(minus):
                r = minus[j]
                length = r.end - f.start
                if length > max_length:
                    break
                if length >= min_length:
                    combo = frozenset({f.primer_id, r.primer_id})
                    yield CrossProduct(
                        left_id=f.primer_id, right_id=r.primer_id,
                        chrom=chrom, start=f.start, end=r.end, length=length,
                        mismatches=f.mismatches + r.mismatches,
                        designed_pair=combo in designed,
                    )
                j += 1


def assemble_cross_products(sites, designed, min_length=50, max_length=3000
                            ) -> list[CrossProduct]:
    """`iter_products`'in liste sarmalayicisi (testler ve kucuk paneller)."""
    return list(iter_products(sites, designed, min_length, max_length))


def dimer_matrix(primers: dict[str, str], cond: salt.Conditions | None = None
                 ) -> list[dict]:
    """MPX-01 ve dimer matrisi: n(n+1)/2 kombinasyon, kendisiyle dahil."""
    cond = cond or salt.Conditions()
    ids = sorted(primers)
    rows = []
    for a, b in itertools.combinations_with_replacement(ids, 2):
        d_any = align.best_dimer(primers[a], primers[b], "any")
        d_3p = align.best_dimer(primers[a], primers[b], "3prime")
        rows.append({
            "primer_a": a, "primer_b": b,
            "dG37_any": round(d_any.dG37_kcal, 3) if d_any.found else None,
            "dG37_3prime": round(d_3p.dG37_kcal, 3) if d_3p.found else None,
            "tm_any_c": round(d_any.tm_c(cond.primer_M), 2) if d_any.found else None,
            "self": a == b,
            "mismatches_any": d_any.mismatch_count if d_any.found else None,
            "mismatches_3prime": d_3p.mismatch_count if d_3p.found else None,
        })
    return rows


def tm_homogeneity(primers: dict[str, str],
                   cond: salt.Conditions | None = None) -> dict:
    """MPX-01: paneldeki Tm araligi."""
    cond = cond or salt.Conditions()
    tms = {pid: concentration.build(seq, cond=cond).tm_c
           for pid, seq in primers.items()}
    lo_id = min(tms, key=tms.get)
    hi_id = max(tms, key=tms.get)
    return {
        "n_primers": len(tms),
        "tm_min_c": round(tms[lo_id], 3), "tm_min_primer": lo_id,
        "tm_max_c": round(tms[hi_id], 3), "tm_max_primer": hi_id,
        "tm_range_c": round(tms[hi_id] - tms[lo_id], 3),
        "per_primer": {k: round(v, 3) for k, v in tms.items()},
    }


def size_separation(lengths: list[int], min_gap_bp: int = 10) -> dict:
    """MPX-03: kapiller elektroforezde ayirt edilebilirlik.

    `lengths` HEDEF urun uzunluklaridir (panelde beyan edilen). Hedef disi
    urunleri buraya katmak olcumu anlamsiz kilar, cunku ayni cift onlarca
    farkli boyda hedef disi urun verebilir.
    """
    lengths = sorted(lengths)
    if len(lengths) < 2:
        return {"n_products": len(lengths), "min_gap_bp": None,
                "conflicts": 0}
    gaps = [b - a for a, b in zip(lengths, lengths[1:])]
    return {
        "n_products": len(lengths),
        "min_gap_bp": min(gaps),
        "conflicts": sum(1 for g in gaps if g < min_gap_bp),
        "threshold_bp": min_gap_bp,
    }
