"""Izoform kapsami (RT-03) dogrulamasi."""
from __future__ import annotations
import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oligoscout_thermo import isoform as iso  # noqa: E402

IDX = iso.TranscriptIndex(by_gene={"G": {
    "T1": ("chr1", [(100, 200), (300, 400)]),
    "T2": ("chr1", [(100, 250)]),
}})


class TestCoverage(unittest.TestCase):
    def test_within_single_exon(self):
        self.assertTrue(iso.covered_by_exons((110, 130), [(100, 200), (300, 400)]))

    def test_intron_not_covered(self):
        self.assertFalse(iso.covered_by_exons((210, 230), [(100, 200), (300, 400)]))

    def test_spanning_intron_gap_not_covered(self):
        self.assertFalse(iso.covered_by_exons((190, 310), [(100, 200), (300, 400)]))

    def test_merge_adjacent(self):
        self.assertTrue(iso.covered_by_exons((150, 350), [(100, 300), (300, 400)]))


class TestRt03(unittest.TestCase):
    def test_targets_junction_isoform(self):
        r = iso.check_rt03("G", "chr1", (110, 130), (350, 370), IDX)
        self.assertEqual(r["value"]["targeted"], ["T1"])
        self.assertEqual(r["value"]["untargeted"], ["T2"])
        self.assertEqual(r["severity"], "INFO")

    def test_targets_single_exon_isoform(self):
        r = iso.check_rt03("G", "chr1", (110, 130), (210, 240), IDX)
        self.assertEqual(r["value"]["targeted"], ["T2"])

    def test_no_annotation(self):
        r = iso.check_rt03("X", "chr1", (1, 2), (3, 4), IDX)
        self.assertEqual(r["value"]["status"], "no_annotation")
        self.assertEqual(r["severity"], "INFO")

    def test_no_isoform_targeted_warns(self):
        # her iki primer de intronda -> hicbir izoform
        r = iso.check_rt03("G", "chr1", (210, 230), (260, 280), IDX)
        self.assertEqual(r["value"]["n_targeted"], 0)
        self.assertEqual(r["severity"], "WARN")

    def test_wrong_chrom_untargeted(self):
        r = iso.check_rt03("G", "chr9", (110, 130), (350, 370), IDX)
        self.assertEqual(r["value"]["n_targeted"], 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
