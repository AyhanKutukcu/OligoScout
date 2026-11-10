"""Varyant katmani testleri: VAR-01..VAR-08 (cevrimdisi, sentetik varyantlar)."""

from __future__ import annotations

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import variants as var  # noqa: E402

# Primer: chr22 uzerinde 0-tabanli [1000, 1020), 20 nt.
FWD = var.Site("P_F", "chr22", 1000, 1020, "+")
REV = var.Site("P_R", "chr22", 1100, 1120, "-")


def snv(pos, af, rsid="", pops=None):
    return var.Variant(chrom="22", pos=pos, ref="A", alt="G", af=af,
                       an=150000, rsid=rsid, pop_af=pops or {})


def indel(pos, af, ref="AT", alt="A"):
    return var.Variant(chrom="22", pos=pos, ref=ref, alt=alt, af=af, an=150000)


def sev(checks, cid):
    order = {"INFO": 0, "WARN": 1, "FAIL": 2}
    return max((c.severity for c in checks if c.id == cid),
               key=lambda s: order[s])


def value(checks, cid):
    return next(c.value for c in checks if c.id == cid)


class TestGeometry(unittest.TestCase):
    def test_chrom_normalisation(self):
        self.assertEqual(var.normalize_chrom("chr22"), "22")
        self.assertEqual(var.normalize_chrom("22"), "22")

    def test_snv_span_is_one_base(self):
        self.assertEqual(snv(1005, 0.1).span, (1004, 1005))

    def test_indel_span_covers_ref(self):
        self.assertEqual(indel(1005, 0.1, "ATTT", "A").span, (1004, 1008))

    def test_three_prime_window_plus_strand(self):
        self.assertEqual(FWD.three_prime_window(5), (1015, 1020))

    def test_three_prime_window_minus_strand(self):
        self.assertEqual(REV.three_prime_window(5), (1100, 1105))

    def test_distance_from_3p_plus(self):
        self.assertEqual(var.distance_from_3p(FWD, snv(1020, 0.1)), 0)
        self.assertEqual(var.distance_from_3p(FWD, snv(1001, 0.1)), 19)

    def test_distance_from_3p_minus(self):
        self.assertEqual(var.distance_from_3p(REV, snv(1101, 0.1)), 0)
        self.assertEqual(var.distance_from_3p(REV, snv(1120, 0.1)), 19)

    def test_overlap_is_half_open(self):
        self.assertEqual(var.overlapping([snv(1000, 0.1)], "chr22", 1000, 1020),
                         [])          # pos 1000 -> span (999,1000), disarida
        self.assertEqual(len(var.overlapping([snv(1001, 0.1)], "chr22",
                                             1000, 1020)), 1)


class TestChecks(unittest.TestCase):
    def test_clean_site_is_info(self):
        checks = var.evaluate_site(FWD, [])
        self.assertEqual(var.worst_severity(checks), "INFO")
        self.assertEqual(value(checks, "VAR-01")["n_variants"], 0)

    def test_var01_warns_on_frequent_variant(self):
        checks = var.evaluate_site(FWD, [snv(1005, 2e-3)])
        self.assertEqual(sev(checks, "VAR-01"), "WARN")

    def test_var01_info_for_singleton(self):
        checks = var.evaluate_site(FWD, [snv(1005, 7e-6)])
        self.assertEqual(sev(checks, "VAR-01"), "INFO")

    def test_var02_fails_for_common_variant_at_3p(self):
        checks = var.evaluate_site(FWD, [snv(1019, 0.2)])
        self.assertEqual(sev(checks, "VAR-02"), "FAIL")
        self.assertEqual(value(checks, "VAR-02")["distance_from_3p_nt"], 1)

    def test_var02_ignores_5p_variant(self):
        checks = var.evaluate_site(FWD, [snv(1002, 0.2)])
        self.assertEqual(sev(checks, "VAR-02"), "INFO")
        self.assertEqual(sev(checks, "VAR-03"), "FAIL")

    def test_var02_uses_correct_end_on_minus_strand(self):
        """- iplikte 3' uc SOLDADIR; yanlis tarafi almak kontrolu bozar."""
        near = var.evaluate_site(REV, [snv(1102, 0.2)])
        far = var.evaluate_site(REV, [snv(1118, 0.2)])
        self.assertEqual(sev(near, "VAR-02"), "FAIL")
        self.assertEqual(sev(far, "VAR-02"), "INFO")

    def test_var03_thresholds(self):
        self.assertEqual(sev(var.evaluate_site(FWD, [snv(1005, 0.02)]),
                             "VAR-03"), "WARN")
        self.assertEqual(sev(var.evaluate_site(FWD, [snv(1005, 0.2)]),
                             "VAR-03"), "FAIL")

    def test_var04_skips_a_site_without_variants(self):
        v = value(var.evaluate_site(FWD, []), "VAR-04")
        self.assertFalse(v["evaluated"])

    def test_var04_scores_both_strands_and_skips_penultimate_variant(self):
        primer = "CAACG"
        forward = var.Site("F", "22", 100, 105, "+", primer)
        reverse = var.Site("R", "22", 110, 115, "-", primer)
        forward_snv = snv(103, 0.01)
        reverse_snv = var.Variant("22", 113, "T", "C", 0.01)
        for site, hit in ((forward, forward_snv), (reverse, reverse_snv)):
            with self.subTest(strand=site.strand):
                result = value(var.evaluate_site(site, [hit]), "VAR-04")
                self.assertEqual(result["n_scored"], 1)
                self.assertGreater(result["variants"][0]["delta_delta_g37_kcal_mol"], 0)
        result = value(var.evaluate_site(forward, [snv(102, 0.01)]), "VAR-04")
        self.assertEqual(result["n_unscored"], 1)

    def test_var05_probe_path(self):
        probe = var.Site("P_P", "chr22", 1030, 1055, "+")
        checks = var.evaluate_site(probe, [snv(1040, 0.2)], probe=True)
        self.assertEqual({c.id for c in checks}, {"VAR-05"})
        self.assertEqual(sev(checks, "VAR-05"), "FAIL")

    def test_var06_singleton_indel_is_info(self):
        checks = var.evaluate_site(FWD, [indel(1005, 7e-6)])
        self.assertEqual(sev(checks, "VAR-06"), "INFO")
        self.assertEqual(value(checks, "VAR-06")["n_indels"], 1)

    def test_var06_common_indel_fails(self):
        checks = var.evaluate_site(FWD, [indel(1005, 0.03)])
        self.assertEqual(sev(checks, "VAR-06"), "FAIL")
        self.assertEqual(value(checks, "VAR-06")["length_change"], -1)

    def test_var07_population_specific(self):
        v = snv(1005, 0.004, "rsPOP", {"afr": 0.08, "nfe": 0.0005})
        checks = var.evaluate_site(FWD, [v])
        self.assertEqual(sev(checks, "VAR-07"), "WARN")
        self.assertEqual(value(checks, "VAR-07")["population"], "afr")

    def test_var07_quiet_when_uniform(self):
        v = snv(1005, 0.05, "rsUNI", {"afr": 0.05, "nfe": 0.05})
        self.assertEqual(sev(var.evaluate_site(FWD, [v]), "VAR-07"), "INFO")


class TestAmplicon(unittest.TestCase):
    def test_var08_excludes_primer_variants(self):
        inside = indel(1060, 0.2, "ATTT", "A")
        on_primer = indel(1005, 0.2, "ATTT", "A")
        checks = var.evaluate_amplicon("P_amp", "chr22", 1000, 1120,
                                       [FWD, REV], [inside, on_primer])
        v = checks[0].value
        self.assertEqual(v["n_common_indels"], 1)
        self.assertEqual(checks[0].severity, "WARN")

    def test_var08_info_without_common_indel(self):
        checks = var.evaluate_amplicon("P_amp", "chr22", 1000, 1120,
                                       [FWD, REV], [snv(1060, 0.4)])
        self.assertEqual(checks[0].severity, "INFO")
        self.assertEqual(checks[0].value["n_variants_inside"], 1)


class TestRecords(unittest.TestCase):
    def test_record_fields_and_summary(self):
        checks = var.evaluate_site(FWD, [snv(1019, 0.2)],
                                   conditions={"dataset": "gnomad_r4"})
        for c in checks:
            self.assertTrue(c.param_set_sha256)
            self.assertEqual(c.conditions["dataset"], "gnomad_r4")
        s = var.summarize(checks)
        self.assertEqual(s["worst_severity"], "FAIL")
        self.assertEqual(s["by_check"]["VAR-02"], "FAIL")
        self.assertEqual(set(s["by_check"]),
                         {"VAR-01", "VAR-02", "VAR-03", "VAR-04", "VAR-06",
                          "VAR-07"})


class TestVerdictCombination(unittest.TestCase):
    def test_verdict_from_severities(self):
        self.assertEqual(var.verdict_from_severities(["INFO", "INFO"]), "PASS")
        self.assertEqual(var.verdict_from_severities(["INFO", "WARN"]), "WARN")
        self.assertEqual(var.verdict_from_severities(["FAIL", "WARN"]), "FAIL")
        self.assertEqual(var.verdict_from_severities([]), "PASS")

    def test_specificity_pass_variant_fail_is_fail(self):
        self.assertEqual(var.combine_verdicts("PASS", "FAIL"), "FAIL")

    def test_empty_layer_is_ignored(self):
        self.assertEqual(var.combine_verdicts("WARN", ""), "WARN")

    def test_worst_of_three(self):
        self.assertEqual(var.combine_verdicts("PASS", "WARN", "PASS"), "WARN")


if __name__ == "__main__":
    unittest.main(verbosity=2)
