"""Tasarim motorunun dogrulamasi."""

from __future__ import annotations

import os
import random
import sys
import unittest
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import design as dsg  # noqa: E402
from oligoscout_thermo import nn, salt  # noqa: E402

RNG = random.Random(11)
TARGET = "".join(RNG.choices("ACGT", weights=[27, 23, 23, 27], k=400))
CONS = dsg.DesignConstraints(product_min=90, product_max=200, max_pairs=300)
COND = salt.Conditions()


class TestCompositionHelpers(unittest.TestCase):
    def test_homopolymer(self):
        self.assertEqual(dsg.longest_homopolymer("AACCCGT"), 3)
        self.assertEqual(dsg.longest_homopolymer("ACGT"), 1)

    def test_repeat_units(self):
        self.assertGreaterEqual(dsg.max_repeat_units("ATATATAT", 2), 4)
        self.assertEqual(dsg.max_repeat_units("ACGTACGT", 2), 1)

    def test_gc_clamp(self):
        self.assertEqual(dsg.gc_clamp("AAAAAGGCCC"), 5)
        self.assertEqual(dsg.gc_clamp("GGGGGAAAAA"), 0)


class TestCandidateGeneration(unittest.TestCase):
    def setUp(self):
        self.acc, self.rej = dsg.enumerate_candidates(TARGET, "F", CONS, COND)

    def test_every_candidate_satisfies_hard_filters(self):
        for c in self.acc:
            with self.subTest(seq=c.sequence):
                self.assertTrue(CONS.length_min <= len(c.sequence) <= CONS.length_max)
                self.assertTrue(CONS.gc_min <= c.gc <= CONS.gc_max)
                self.assertTrue(CONS.tm_min <= c.tm_c <= CONS.tm_max)
                self.assertTrue(
                    CONS.gc_clamp_min <= dsg.gc_clamp(c.sequence) <= CONS.gc_clamp_max
                )
                self.assertLess(dsg.longest_homopolymer(c.sequence),
                                CONS.max_homopolymer)
                self.assertGreaterEqual(c.self_dimer_dg, CONS.self_dimer_dg)

    def test_every_rejection_has_a_rule_and_value(self):
        self.assertTrue(self.rej)
        for r in self.rej:
            self.assertTrue(r.rule)
            self.assertTrue(r.reason)
            self.assertNotEqual(r.value, "")

    def test_reverse_candidates_are_reverse_complements(self):
        acc_r, _ = dsg.enumerate_candidates(TARGET, "R", CONS, COND)
        for c in acc_r[:20]:
            window = TARGET[c.start:c.end]
            self.assertEqual(c.sequence, nn.reverse_complement(window))

    def test_candidate_and_rejection_counts_are_consistent(self):
        """Her uretilen aday ya kabul edilir ya gerekceyle elenir."""
        total_windows = sum(
            1
            for start in range(len(TARGET))
            for length in range(CONS.length_min, CONS.length_max + 1)
            if start + length <= len(TARGET)
        )
        self.assertEqual(len(self.acc) + len(self.rej), total_windows)


class TestBindingFilter(unittest.TestCase):
    def test_binding_threshold_accepts_equality_and_records_rejection(self):
        primer = "CTGCTACGGTGTGAATGCTAG"
        cons = dsg.DesignConstraints(length_min=len(primer), length_max=len(primer))
        accepted, _ = dsg.enumerate_candidates(primer, "F", cons, COND)
        candidate, = accepted
        at_threshold, _ = dsg.enumerate_candidates(
            primer, "F", replace(cons, theta_min=candidate.theta_ta), COND
        )
        self.assertEqual(at_threshold, accepted)
        below_threshold, rejections = dsg.enumerate_candidates(
            primer, "F", replace(cons, theta_min=candidate.theta_ta + 0.01), COND
        )
        self.assertEqual(below_threshold, [])
        self.assertEqual([r.rule for r in rejections], ["TH-03"])


class TestPairing(unittest.TestCase):
    def setUp(self):
        self.result = dsg.design(TARGET, CONS, COND)

    def test_products_within_window(self):
        for p in self.result["pairs"]:
            self.assertTrue(
                CONS.product_min <= p.product_length <= CONS.product_max
            )
            self.assertLessEqual(p.delta_tm, CONS.max_delta_tm)
            self.assertGreaterEqual(p.hetero_dg, CONS.hetero_dimer_dg)

    def test_primers_do_not_overlap(self):
        for p in self.result["pairs"]:
            self.assertLessEqual(p.left.end, p.right.start)

    def test_pareto_front_is_non_dominated(self):
        pairs = self.result["pairs"]
        front = [p for p in pairs if p.pareto_rank == 1]
        self.assertTrue(front)
        for a in front:
            for b in pairs:
                if b is a:
                    continue
                dominates = (
                    all(x <= y for x, y in zip(b.objectives, a.objectives))
                    and any(x < y for x, y in zip(b.objectives, a.objectives))
                )
                self.assertFalse(dominates, f"{b.left.sequence} dominates front")

    def test_ranking_is_sorted(self):
        ranks = [(p.pareto_rank, p.score) for p in self.result["pairs"]]
        self.assertEqual(ranks, sorted(ranks))

    def test_calibration_status_declared(self):
        status = self.result["calibration_status"]
        self.assertFalse(status["ranking_score_calibrated"])
        self.assertFalse(status["specificity_screened"])

    def test_tight_constraints_yield_nothing_but_log_reasons(self):
        tight = dsg.DesignConstraints(tm_min=90.0, tm_max=95.0)
        out = dsg.design(TARGET, tight, COND)
        self.assertEqual(out["pairs"], [])
        self.assertTrue(any(r.rule == "TH-01" for r in out["rejections"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
