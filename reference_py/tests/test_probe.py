"""Hidroliz probu modulu testleri: PRB-01..PRB-05 ve prob secimi."""

from __future__ import annotations

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import nn, probe as prb, salt  # noqa: E402

FWD = "CTCCTCAAAGGCACCTACT"
REV = "GGGCACTAATCACCTCCTAA"
MIDDLE = ("ATCCACCTGAGTCTGCAATCCAGACCCTTAGGTACTCCAGACTTCCAAGTGCCATCC"
          "AGTTGGCTTACCTCAGGACCTTGCCATAGACCTTAGGACCTCCA")
AMPLICON = FWD + MIDDLE + nn.reverse_complement(REV)


def severity_of(result, check_id):
    """Bir kontrolun en agir ciddiyeti (PRB-05 birden cok satir uretir)."""
    order = {"INFO": 0, "WARN": 1, "FAIL": 2}
    sevs = [c["severity"] for c in result["checks"] if c["id"] == check_id]
    return max(sevs, key=lambda s: order[s])


def value_of(result, check_id):
    for c in result["checks"]:
        if c["id"] == check_id:
            return c["value"]
    raise KeyError(check_id)


class TestProbeSelection(unittest.TestCase):
    def setUp(self):
        self.cand = prb.pick_probe(AMPLICON, FWD, REV)

    def test_probe_found(self):
        self.assertIsNotNone(self.cand)

    def test_probe_inside_inner_region(self):
        """PRB-04: secilen prob primer baglanma yerlerine girmez."""
        self.assertGreaterEqual(self.cand.start, len(FWD))
        self.assertLessEqual(self.cand.end, len(AMPLICON) - len(REV))

    def test_probe_tm_window(self):
        """PRB-01: secilen prob 6-10 C penceresinde."""
        self.assertGreaterEqual(self.cand.delta_tm_c, 6.0)
        self.assertLessEqual(self.cand.delta_tm_c, 10.0)

    def test_composition_rules_applied(self):
        """PRB-02 ve PRB-03 tasarim asamasinda zorunlu."""
        self.assertNotEqual(self.cand.sequence[0], "G")
        self.assertGreater(self.cand.c_count, self.cand.g_count)

    def test_rejections_are_counted(self):
        cands, rejections = prb.enumerate_probes(AMPLICON, FWD, REV)
        self.assertTrue(rejections)
        self.assertEqual(cands[0].sequence, self.cand.sequence)

    def test_no_room_between_primers(self):
        tight = FWD + "ACGT" + nn.reverse_complement(REV)
        cands, rej = prb.enumerate_probes(tight, FWD, REV)
        self.assertEqual(cands, [])
        self.assertEqual(rej[0][0], "PRB-04")


class TestProbeChecks(unittest.TestCase):
    def evaluate(self, probe_seq, amplicon=AMPLICON):
        return prb.evaluate_probe("P", probe_seq, FWD, REV, amplicon=amplicon)

    def test_all_pass_for_designed_probe(self):
        cand = prb.pick_probe(AMPLICON, FWD, REV)
        res = self.evaluate(cand.sequence)
        self.assertEqual(res["worst_severity"], "INFO")

    def test_prb01_probe_too_cold(self):
        """Primerlerle ayni Tm'deki bir prob PRB-01'de FAIL."""
        res = self.evaluate(MIDDLE[10:28])
        v = value_of(res, "PRB-01")
        self.assertLess(v["delta_tm_c"], 3.0)
        self.assertEqual(severity_of(res, "PRB-01"), "FAIL")

    def test_prb01_window_warn(self):
        """Pencerenin hemen disi (3-6 C) WARN, FAIL degil."""
        cand = prb.pick_probe(AMPLICON, FWD, REV)
        shorter = cand.sequence
        while True:
            res = self.evaluate(shorter[:-1])
            if value_of(res, "PRB-01")["delta_tm_c"] < 6.0:
                break
            shorter = shorter[:-1]
        self.assertEqual(severity_of(res, "PRB-01"), "WARN")

    def test_prb02_five_prime_g(self):
        probe_seq = "G" + MIDDLE[1:25]
        res = self.evaluate(probe_seq, amplicon=None)
        self.assertEqual(severity_of(res, "PRB-02"), "FAIL")
        self.assertEqual(value_of(res, "PRB-02")["five_prime_base"], "G")

    def test_prb03_more_g_than_c(self):
        probe_seq = "ATGGATGGATGGATGGATGG"
        res = self.evaluate(probe_seq, amplicon=None)
        v = value_of(res, "PRB-03")
        self.assertGreater(v["g_count"], v["c_count"])
        self.assertEqual(severity_of(res, "PRB-03"), "WARN")
        self.assertEqual(v["g_count_other_strand"], v["c_count"])

    def test_prb04_overlap_with_primer_fails(self):
        """Ileri primerin ustune binen prob FAIL verir."""
        probe_seq = AMPLICON[len(FWD) - 6:len(FWD) + 18]
        res = self.evaluate(probe_seq)
        v = value_of(res, "PRB-04")
        self.assertEqual(v["overlap_forward_nt"], 6)
        self.assertEqual(severity_of(res, "PRB-04"), "FAIL")

    def test_prb04_probe_not_in_amplicon(self):
        res = self.evaluate("ACCTTACCTTACCTTACCTTAC")
        self.assertFalse(value_of(res, "PRB-04")["found_in_amplicon"])
        self.assertEqual(severity_of(res, "PRB-04"), "FAIL")

    def test_prb04_reverse_strand_probe_found(self):
        cand = prb.pick_probe(AMPLICON, FWD, REV)
        rc = nn.reverse_complement(cand.sequence)
        v = value_of(self.evaluate(rc), "PRB-04")
        self.assertTrue(v["found_in_amplicon"])
        self.assertEqual(v["strand"], "-")

    def test_prb04_not_evaluated_without_amplicon(self):
        res = self.evaluate(MIDDLE[20:44], amplicon=None)
        self.assertFalse(value_of(res, "PRB-04")["evaluated"])

    def test_prb05_self_complementary_probe(self):
        """Kendine tumleyen prob PRB-05'te FAIL."""
        half = "CCTTAGCAGCTT"
        probe_seq = half + nn.reverse_complement(half)
        res = self.evaluate(probe_seq, amplicon=None)
        self.assertEqual(severity_of(res, "PRB-05"), "FAIL")

    def test_prb05_reports_all_partners(self):
        cand = prb.pick_probe(AMPLICON, FWD, REV)
        res = self.evaluate(cand.sequence)
        targets = {c["target"].split(":")[1] for c in res["checks"]
                   if c["id"] == "PRB-05"}
        self.assertEqual(targets, {"self", "forward", "reverse"})


class TestEvidenceRecords(unittest.TestCase):
    def test_record_fields(self):
        cand = prb.pick_probe(AMPLICON, FWD, REV)
        res = prb.evaluate_probe("P", cand.sequence, FWD, REV,
                                 amplicon=AMPLICON,
                                 cond=salt.Conditions(anneal_c=58.0))
        for c in res["checks"]:
            self.assertTrue(c["param_set_sha256"])
            self.assertTrue(c["engine_version"])
            self.assertEqual(c["conditions"]["anneal_c"], 58.0)
        self.assertFalse(res["calibration_status"]["thresholds_calibrated"])

    def test_all_five_checks_present(self):
        cand = prb.pick_probe(AMPLICON, FWD, REV)
        res = prb.evaluate_probe("P", cand.sequence, FWD, REV,
                                 amplicon=AMPLICON)
        ids = {c["id"] for c in res["checks"]}
        self.assertEqual(ids, {"PRB-01", "PRB-02", "PRB-03", "PRB-04",
                               "PRB-05"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
