"""Amplikon erime modeli testleri (QPCR-02)."""

from __future__ import annotations

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import amplicon as am  # noqa: E402
from oligoscout_thermo import salt  # noqa: E402

AT = "ATATTAATATTAAATATTAATATTAAATATTAATATTAAATATTAATA"
GC = "GCGCGGCCGCGGCGCCGCGGCCGCGGCGCCGCGGCCGCGGCGCCGCGG"
MIXED = ("CTCCTCAAAGGCACCTACTGGCACAGTCCACACTGCTACGGTGTGAATGCTAGTGT"
         "CCAGACCCTTAGGTACTCCAGACTTCCAAGTGCCATCCAGTTGGCTTACCTCAGG")


class TestProfileShape(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = am.melt_profile(MIXED, t_min=40, t_max=110, step=0.5)

    def test_theta_is_monotone_decreasing(self):
        self.assertTrue(all(b <= a + 1e-9
                            for a, b in zip(self.p.helicity, self.p.helicity[1:])))

    def test_theta_spans_paired_to_melted(self):
        self.assertGreater(self.p.helicity[0], 0.95)
        self.assertLess(self.p.helicity[-1], 0.05)

    def test_tm_is_inside_scan(self):
        self.assertTrue(40 <= self.p.tm_c <= 110)

    def test_tm_brackets_half_helicity(self):
        """Tm, theta'nin 0.5'i gectigi iki izgara noktasi arasinda olmali.

        Gecis dik oldugu icin en yakin izgara noktasinda theta 0.5'ten
        belirgin sapabilir; dogru olcut arayi kusatmasidir.
        """
        below = max(i for i, t in enumerate(self.p.temperatures)
                    if t <= self.p.tm_c)
        self.assertGreaterEqual(self.p.helicity[below], 0.5)
        self.assertLess(self.p.helicity[below + 1], 0.5)

    def test_too_short_raises(self):
        with self.assertRaises(ValueError):
            am.melt_profile("A")


class TestPhysics(unittest.TestCase):
    def test_gc_rich_melts_later(self):
        at = am.melt_profile(AT, t_min=20, t_max=130, step=0.5)
        gc = am.melt_profile(GC, t_min=20, t_max=130, step=0.5)
        self.assertGreater(gc.tm_c, at.tm_c + 20)

    def test_salt_raises_tm(self):
        lo = am.melt_profile(AT, cond=salt.Conditions(monovalent_M=0.01),
                             t_min=20, t_max=130, step=0.5)
        hi = am.melt_profile(AT, cond=salt.Conditions(monovalent_M=0.5),
                             t_min=20, t_max=130, step=0.5)
        self.assertGreater(hi.tm_c, lo.tm_c)

    def test_longer_duplex_melts_later(self):
        short = am.melt_profile(MIXED[:30], t_min=20, t_max=130, step=0.5)
        long = am.melt_profile(MIXED, t_min=20, t_max=130, step=0.5)
        self.assertGreater(long.tm_c, short.tm_c)

    def test_two_domains_give_two_peaks(self):
        """AT ve GC alan adi birlesince iki ayri pik cikar."""
        p = am.melt_profile(AT + GC, t_min=20, t_max=130, step=0.5)
        self.assertEqual(len(p.peaks), 2)
        at = am.melt_profile(AT, t_min=20, t_max=130, step=0.5)
        gc = am.melt_profile(GC, t_min=20, t_max=130, step=0.5)
        lo, hi = p.peaks[0][0], p.peaks[1][0]
        self.assertLess(abs(lo - at.tm_c), 12)
        self.assertLess(abs(hi - gc.tm_c), 12)

    def test_homogeneous_sequence_single_peak(self):
        p = am.melt_profile(GC, t_min=20, t_max=130, step=0.5)
        self.assertEqual(len(p.peaks), 1)

    def test_deterministic(self):
        a = am.melt_profile(MIXED, t_min=60, t_max=95, step=0.5)
        b = am.melt_profile(MIXED, t_min=60, t_max=95, step=0.5)
        self.assertEqual(a.helicity, b.helicity)


class TestQpcr02(unittest.TestCase):
    def test_single_peak_is_info(self):
        p = am.melt_profile(GC, t_min=20, t_max=130, step=0.5)
        self.assertEqual(am.qpcr02(p)["severity"], "INFO")

    def test_two_separated_peaks_warn(self):
        p = am.melt_profile(AT + GC, t_min=20, t_max=130, step=0.5)
        v = am.qpcr02(p)
        self.assertEqual(v["severity"], "WARN")
        self.assertEqual(v["n_peaks"], 2)

    def test_separation_threshold_is_respected(self):
        p = am.melt_profile(AT + GC, t_min=20, t_max=130, step=0.5)
        self.assertEqual(am.qpcr02(p, min_separation_c=200.0)["severity"], "INFO")

    def test_model_constants_declared_uncalibrated(self):
        p = am.melt_profile(MIXED, t_min=60, t_max=95, step=1.0)
        d = p.as_dict()
        self.assertFalse(d["model_constants"]["calibrated"])
        self.assertIn("pending", d["loop_entropy"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
