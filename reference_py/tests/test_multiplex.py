"""Multipleks capraz reaktivite ve uygulama profilleri testleri."""

from __future__ import annotations

import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import multiplex as mpx  # noqa: E402
from oligoscout_thermo import profiles  # noqa: E402

F1 = "GTCACTGGCATTCAAGAGCTACAT"
R1 = "ACCTGTGACATTGCCTGGAGATTA"
F2 = "CAAGTGCCTAAGGTCATCGGATTC"
R2 = "TTGCACCAGTGACTTAGCCATGAA"


def write_bindings(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pair_id", "role", "primer", "chrom", "start",
                    "end_exclusive", "strand", "mismatches", "mask_bits",
                    "mask", "site_plus_strand"])
        for r in rows:
            w.writerow(r)


class TestCrossProducts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, "p.bindings.tsv")
        # P0: F1 @1000(+), R1 @1200(-)   P1: F2 @1100(+), R2 @1400(-)
        write_bindings(self.path, [
            ["P0", "F", F1, "chr1", 1000, 1024, "+", 0, 0, "." * 24, F1],
            ["P0", "R", R1, "chr1", 1200, 1224, "-", 0, 0, "." * 24, R1],
            ["P1", "F", F2, "chr1", 1100, 1124, "+", 0, 0, "." * 24, F2],
            ["P1", "R", R2, "chr1", 1400, 1424, "-", 0, 0, "." * 24, R2],
        ])
        self.sites, self.owners = mpx.load_bindings(self.path)

    def test_identity_is_sequence(self):
        self.assertEqual({s.primer_id for s in self.sites},
                         {F1, R1, F2, R2})

    def test_designed_and_cross_split(self):
        designed = {frozenset({F1, R1}), frozenset({F2, R2})}
        prods = mpx.assemble_cross_products(self.sites, designed, 50, 3000)
        designed_p = [p for p in prods if p.designed_pair]
        cross = [p for p in prods if not p.designed_pair]
        self.assertEqual(len(designed_p), 2)          # P0 ve P1 kendi urunleri
        # F1 x R2 ve F2 x R1 capraz urun verir
        self.assertEqual(len(cross), 2)
        self.assertEqual({(p.left_id, p.right_id) for p in cross},
                         {(F1, R2), (F2, R1)})

    def test_length_window_respected(self):
        designed = {frozenset({F1, R1}), frozenset({F2, R2})}
        prods = mpx.assemble_cross_products(self.sites, designed, 300, 3000)
        for p in prods:
            self.assertGreaterEqual(p.length, 300)
        self.assertTrue(all(p.length >= 300 for p in prods))

    def test_shared_primer_keeps_designed_mapping(self):
        """Ayni primeri paylasan iki cift 'capraz' sayilmamali."""
        path = os.path.join(self.tmp, "s.bindings.tsv")
        write_bindings(path, [
            ["P0", "F", F1, "chr1", 1000, 1024, "+", 0, 0, "." * 24, F1],
            ["P0", "R", R1, "chr1", 1200, 1224, "-", 0, 0, "." * 24, R1],
            ["P1", "F", F1, "chr1", 1000, 1024, "+", 0, 0, "." * 24, F1],
            ["P1", "R", R2, "chr1", 1400, 1424, "-", 0, 0, "." * 24, R2],
        ])
        sites, owners = mpx.load_bindings(path)
        self.assertEqual(owners[F1], {"P0_F", "P1_F"})
        designed = {frozenset({F1, R1}), frozenset({F1, R2})}
        prods = mpx.assemble_cross_products(sites, designed, 50, 3000)
        self.assertTrue(all(p.designed_pair for p in prods), prods)

    def test_no_product_on_different_chromosomes(self):
        path = os.path.join(self.tmp, "c.bindings.tsv")
        write_bindings(path, [
            ["P0", "F", F1, "chr1", 1000, 1024, "+", 0, 0, "." * 24, F1],
            ["P0", "R", R1, "chr2", 1200, 1224, "-", 0, 0, "." * 24, R1],
        ])
        sites, _ = mpx.load_bindings(path)
        self.assertEqual(mpx.assemble_cross_products(sites, set(), 50, 3000), [])


class TestMatrixAndSeparation(unittest.TestCase):
    def test_dimer_matrix_shape(self):
        primers = {F1: F1, R1: R1, F2: F2}
        rows = mpx.dimer_matrix(primers)
        self.assertEqual(len(rows), 6)            # n(n+1)/2 = 6
        self.assertEqual(sum(1 for r in rows if r["self"]), 3)

    def test_tm_homogeneity(self):
        out = mpx.tm_homogeneity({F1: F1, R1: R1})
        self.assertEqual(out["n_primers"], 2)
        self.assertGreaterEqual(out["tm_range_c"], 0.0)

    def test_size_separation_uses_given_lengths(self):
        out = mpx.size_separation([100, 110, 112], min_gap_bp=10)
        self.assertEqual(out["min_gap_bp"], 2)
        self.assertEqual(out["conflicts"], 1)
        self.assertEqual(mpx.size_separation([100])["min_gap_bp"], None)


class TestProfiles(unittest.TestCase):
    def test_all_profiles_well_formed(self):
        for name, p in profiles.PROFILES.items():
            with self.subTest(profile=name):
                self.assertEqual(p.name, name)
                c = p.constraints
                self.assertLess(c.product_min, c.product_max)
                self.assertLess(c.length_min, c.length_max)
                self.assertLess(c.gc_min, c.gc_max)
                self.assertLess(p.qc.tm_window_c[0], p.qc.tm_window_c[1])

    def test_ctdna_is_short_amplicon(self):
        self.assertLessEqual(profiles.get("ctdna").constraints.product_max, 100)

    def test_unknown_profile_raises(self):
        with self.assertRaises(KeyError):
            profiles.get("yok_boyle_bir_profil")

    def test_intact_fragment_probability(self):
        """CT-02: P = (L - A + 1) / L."""
        self.assertAlmostEqual(profiles.intact_fragment_probability(70, 166),
                               (166 - 70 + 1) / 166, places=12)
        self.assertEqual(profiles.intact_fragment_probability(200, 166), 0.0)
        self.assertGreater(profiles.intact_fragment_probability(60),
                           profiles.intact_fragment_probability(100))

    def test_detectability_monotone(self):
        """CT-03: kisa amplikon ve yuksek VAF tespit olasiligini artirir."""
        a = profiles.detectability(70, 0.001)
        b = profiles.detectability(150, 0.001)
        self.assertGreater(a["p_detect"], b["p_detect"])
        c = profiles.detectability(70, 0.01)
        self.assertGreater(c["p_detect"], a["p_detect"])

    def test_fragment_weighting_drops_long_products(self):
        """166 bp fragmandan 200 bp urun olusamaz; agirlik tam sifir."""
        self.assertEqual(profiles.intact_fragment_probability(200, 166), 0.0)
        self.assertAlmostEqual(
            profiles.intact_fragment_probability(166, 166), 1 / 166, places=12)

    def test_genome_equivalents(self):
        self.assertAlmostEqual(profiles.genome_equivalents(10.0), 10 / 3.3e-3,
                               places=6)

    def test_detectability_with_distribution(self):
        flat = profiles.detectability(70, 0.001,
                                      fragment_distribution={166: 1.0})
        single = profiles.detectability(70, 0.001, fragment_bp=166)
        self.assertAlmostEqual(flat["p_detect"], single["p_detect"], places=12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
