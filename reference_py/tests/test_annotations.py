"""Aralik anotasyon indeksinin dogrulamasi (SPC-06/07 altyapisi)."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import annotations  # noqa: E402


def make_index(rows):
    return annotations.IntervalIndex({
        chrom: ([r[0] for r in sorted(v)], [r[1] for r in sorted(v)],
                [r[2] for r in sorted(v)])
        for chrom, v in _group(rows).items()
    })


def _group(rows):
    out: dict = {}
    for chrom, s, e, name in rows:
        out.setdefault(chrom, []).append((s, e, name))
    return out


class TestIntervalIndex(unittest.TestCase):
    def setUp(self):
        self.idx = make_index([
            ("chr1", 100, 200, "A"),
            ("chr1", 150, 250, "B"),
            ("chr1", 400, 500, "C"),
            ("chr2", 10, 20, "D"),
        ])

    def test_half_open_overlap(self):
        # [180,190) chr1 -> A ve B ile ortusur.
        self.assertEqual(set(self.idx.overlaps("chr1", 180, 190)), {"A", "B"})

    def test_touching_end_is_not_overlap(self):
        # A [100,200) sorgu [200,210) ile yalniz DOKUNUR -> ortusmez.
        # B [150,250) ise ayni sorguyla gercekten ortusur.
        self.assertEqual(self.idx.overlaps("chr1", 200, 210), ["B"])
        self.assertNotIn("A", self.idx.overlaps("chr1", 200, 210))

    def test_gap_between_intervals(self):
        self.assertEqual(self.idx.overlaps("chr1", 300, 350), [])

    def test_unknown_chrom_is_empty(self):
        self.assertEqual(self.idx.overlaps("chrX", 100, 200), [])
        self.assertFalse(self.idx.overlaps_any("chrX", 100, 200))

    def test_other_chrom_isolated(self):
        self.assertEqual(self.idx.overlaps("chr2", 15, 18), ["D"])
        self.assertEqual(self.idx.overlaps("chr2", 100, 200), [])

    def test_total_count(self):
        self.assertEqual(self.idx.total, 4)


class TestLoadIntervals(unittest.TestCase):
    def test_missing_file_is_empty_index(self):
        idx = annotations.load_intervals("/yok/olmayan.tsv")
        self.assertEqual(idx.total, 0)

    def test_roundtrip_from_tsv(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "iv.tsv")
            with open(p, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("chrom\tstart\tend\tname\n")
                fh.write("chr3\t1000\t2000\tPG1\n")
                fh.write("chr3\t1500\t2500\tPG2\n")
            idx = annotations.load_intervals(p)
            self.assertEqual(idx.total, 2)
            self.assertEqual(set(idx.overlaps("chr3", 1800, 1900)),
                             {"PG1", "PG2"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
