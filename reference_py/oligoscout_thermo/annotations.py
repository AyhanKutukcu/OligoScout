"""Aralik anotasyonlarini hizli ortusme sorgusuyla yukler (SPC-06, SPC-07).

Kompakt araliklar `cli/build_annotations.py` tarafindan uretilir
(chrom<TAB>start<TAB>end<TAB>name, 0-tabanli, bitis haric). Bu modul onlari
kromozom basina siralanmis listelere alir ve ikili arama ile ortusme
sorgular. Ham GTF/superDups'u ayristirmaz; o is on-islemdedir.

Bu modul AG ERISIMI YAPMAZ ve ham dosya yoksa bos index doner (kurallar
sessizce 'anotasyon yok' raporlar).
"""

from __future__ import annotations

import bisect
import os
from dataclasses import dataclass


@dataclass
class IntervalIndex:
    """Kromozom -> (siralanmis start listesi, end listesi, ad listesi)."""

    by_chrom: dict[str, tuple[list[int], list[int], list[str]]]

    @property
    def total(self) -> int:
        return sum(len(v[0]) for v in self.by_chrom.values())

    def overlaps(self, chrom: str, start: int, end: int) -> list[str]:
        """[start, end) ile ortusen araliklarin adlari (bos olabilir)."""
        entry = self.by_chrom.get(chrom)
        if not entry:
            return []
        starts, ends, names = entry
        # Baslangici `end`den kucuk olan tum araliklar aday; bunlardan bitisi
        # `start`ten buyuk olanlar gercekten ortusur.
        hi = bisect.bisect_left(starts, end)
        hits = []
        for i in range(hi):
            if ends[i] > start:
                hits.append(names[i])
        return hits

    def overlaps_any(self, chrom: str, start: int, end: int) -> bool:
        return bool(self.overlaps(chrom, start, end))


def load_intervals(path: str) -> IntervalIndex:
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    if not path or not os.path.isfile(path):
        return IntervalIndex({})
    with open(path, encoding="utf-8") as fh:
        header = fh.readline()
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            name = cols[3] if len(cols) > 3 else ""
            by_chrom.setdefault(chrom, []).append((start, end, name))
    packed: dict[str, tuple[list[int], list[int], list[str]]] = {}
    for chrom, rows in by_chrom.items():
        rows.sort()
        packed[chrom] = ([r[0] for r in rows], [r[1] for r in rows],
                         [r[2] for r in rows])
    return IntervalIndex(packed)


def default_paths(repo_root: str) -> dict[str, str]:
    base = os.path.join(repo_root, "OligoScout", "data", "annotations")
    return {
        "pseudogene": os.path.join(base, "pseudogene_intervals.tsv"),
        "segdup": os.path.join(base, "segdup_intervals.tsv"),
    }
