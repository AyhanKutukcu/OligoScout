"""Bisulfit katmani kontrollerinin dogrulamasi (BIS-01..BIS-04)."""

from __future__ import annotations

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import bisulfite as b  # noqa: E402


class TestHelpers(unittest.TestCase):
    def test_count_cpg(self):
        self.assertEqual(b.count_cpg("ACGCGT"), 2)
        self.assertEqual(b.count_cpg("AAAA"), 0)

    def test_non_cpg_c(self):
        # ACCGT: C@1 (next C) sayilir, C@2 (next G=CpG) sayilmaz -> 1
        self.assertEqual(b.count_non_cpg_c("ACCGT"), 1)
        self.assertEqual(b.count_non_cpg_c("CG"), 0)      # CpG C sayilmaz
        self.assertEqual(b.count_non_cpg_c("CA"), 1)

    def test_nearest_cpg_3p(self):
        self.assertEqual(b.nearest_cpg_distance_from_3p("ACGTT"), 3)
        self.assertIsNone(b.nearest_cpg_distance_from_3p("AAAA"))


class TestBis01(unittest.TestCase):
    def test_no_cpg_info(self):
        self.assertEqual(b.check_bis01("ATATATATAT")["severity"], "INFO")

    def test_cpg_near_3p_fail(self):
        self.assertEqual(b.check_bis01("ATATATATCG")["severity"], "FAIL")

    def test_cpg_mid_warn(self):
        self.assertEqual(b.check_bis01("ATCGATATATAT")["severity"], "WARN")

    def test_msp_exempt(self):
        self.assertEqual(b.check_bis01("ATATATATCG", msp=True)["severity"], "INFO")


class TestBis02(unittest.TestCase):
    def test_low_selectivity_warn(self):
        # CpG disi C yok -> dusuk donusum seciciligi
        self.assertEqual(b.check_bis02("ATGTATGTATGT")["severity"], "WARN")

    def test_good_selectivity_info(self):
        self.assertEqual(b.check_bis02("ACACACTATACA")["severity"], "INFO")
        self.assertGreaterEqual(
            b.check_bis02("ACACACTATACA")["value"]["non_cpg_c_count"], 2)


class TestBis03(unittest.TestCase):
    def test_low_complexity_fail(self):
        self.assertEqual(b.check_bis03("TTTTTTTTTTTT")["severity"], "FAIL")

    def test_normal_converted_info(self):
        # 3-harfli dengeli donusmus dizi normalde INFO
        self.assertEqual(b.check_bis03("ATGTAGTATGATGTAG")["severity"], "INFO")


class TestBis04(unittest.TestCase):
    def test_clean_single_strand(self):
        self.assertEqual(b.check_bis04(1, 1, 0, 0, "CT")["severity"], "INFO")

    def test_cross_strand_fail(self):
        self.assertEqual(b.check_bis04(1, 1, 0, 1, "CT")["severity"], "FAIL")

    def test_extra_intended_warn(self):
        self.assertEqual(b.check_bis04(1, 2, 0, 0, "CT")["severity"], "WARN")

    def test_ga_intended(self):
        # GA tasarimi: CT'de urun iplik-caprazdir -> FAIL
        self.assertEqual(b.check_bis04(1, 1, 1, 1, "GA")["severity"], "FAIL")


if __name__ == "__main__":
    unittest.main(verbosity=2)
