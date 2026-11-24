
"""Regression tests for target coverage and cross-product conflict handling."""
import itertools
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reference_py"))
from multiplex_design import assess_candidates, edge_key, parse_targets, remap_ids, select_pool
from oligoscout_thermo.multiplex import BindingSite


def candidate(pid, target, f="A" * 20, r="T" * 20):
    return {"id": pid, "target_id": target, "F": f, "R": r,
            "design": {"left": {"tm_c": 60}, "right": {"tm_c": 60}, "product_length": 100}}


class MultiplexDesignTests(unittest.TestCase):
    def test_selects_alternative_to_cover_every_target(self):
        a, alternate, b = candidate("A", "T1"), candidate("A2", "T1"), candidate("B", "T2")
        selected, report = select_pool([[a, alternate], [b]], {edge_key("A", "B")})
        self.assertEqual({p["id"] for p in selected}, {"A2", "B"})
        self.assertFalse(report["limit_reached"])

    def test_partial_pool_never_contains_blocked_edges(self):
        groups = [[candidate("A", "T1")], [candidate("B", "T2")], [candidate("C", "T3")]]
        blocked = {edge_key("A", "B"), edge_key("B", "C")}
        selected, _ = select_pool(groups, blocked)
        self.assertEqual({p["id"] for p in selected}, {"A", "C"})

    def test_candidate_failing_its_own_pool_constraints_is_excluded(self):
        a, b = candidate("A", "T1"), candidate("B", "T2")
        selected, _ = select_pool([[a], [b]], {edge_key("A", "A")})
        self.assertEqual([p["id"] for p in selected], ["B"])

    def test_per_target_coordinates_are_used_in_exported_scan_panel(self):
        import csv
        import tempfile
        from local_engine import LocalEngine
        motor = object.__new__(LocalEngine)
        sources = {"T1": {"template": {"chrom": "chr1", "offset": 1000, "conversion": "primary"}},
                   "T2": {"template": {"chrom": "chr2", "offset": 5000, "conversion": "primary"}}}
        source = {**sources["T1"], "target_sources": sources}
        pairs = [{**candidate("A", "T1"), "start": 10, "end": 100},
                 {**candidate("B", "T2"), "start": 20, "end": 120}]
        with tempfile.TemporaryDirectory() as folder:
            path = motor.write_scan_panel(source, pairs, Path(folder))
            with path.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["expected_chrom"], "chr1")
        self.assertEqual(rows[1]["expected_chrom"], "chr2")
        self.assertEqual(rows[1]["left_start"], "5020")
        self.assertEqual(rows[1]["right_start"], "5120")

    def test_empty_target_prevents_full_coverage(self):
        selected, _ = select_pool([[], [candidate("A", "T2")]], set())
        self.assertEqual(len(selected), 1)

    def test_limit_returns_only_a_compatible_partial_pool(self):
        groups = [[candidate("A", "T1")], [candidate("B", "T2")]]
        selected, report = select_pool(groups, set(), limit=1)
        self.assertTrue(report["limit_reached"])
        self.assertLess(len(selected), len(groups))

    def test_duplicate_target_regions_are_rejected(self):
        with self.assertRaises(ValueError):
            parse_targets([{"name": "one", "region": "22:100-500"}, {"name": "two", "region": "chr22:100-500"}])

    def test_duplicate_names_and_unsafe_names_are_rejected(self):
        for name in ("A", "a", "bad\tname"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                parse_targets([{"name": "A", "region": "chr1:100-500"}, {"name": name, "region": "chr2:100-500"}])

    def test_renames_evidence_without_mutating_sequence_or_unrelated_ids(self):
        obj = {"WEB_001": {"primer_id": "WEB_001_F", "sequence": "ACGT", "other": "WEB_0010"}}
        result = remap_ids(obj, {"WEB_001": "T01_WEB_001"})
        self.assertEqual(result["T01_WEB_001"]["primer_id"], "T01_WEB_001_F")
        self.assertEqual(result["T01_WEB_001"]["sequence"], "ACGT")
        self.assertEqual(result["T01_WEB_001"]["other"], "WEB_0010")

    def test_single_primer_genomic_product_excludes_its_candidate(self):
        p = candidate("A", "T1")
        sites = [BindingSite(p["F"], p["F"], "chr1", 100, 120, "+", 0),
                 BindingSite(p["F"], p["F"], "chr1", 160, 180, "-", 0)]
        seqs = sorted({p["F"], p["R"]})
        rows = [{"primer_a": a, "primer_b": b, "dG37_any": 0, "dG37_3prime": 0}
                for a, b in itertools.combinations_with_replacement(seqs, 2)]
        settings = {"min_product": 50, "max_product": 3000, "dimer_any": -6,
                    "dimer_3p": -5, "tm_spread": 3, "product_gap": 0}
        with patch("multiplex_design.mpx.load_bindings", return_value=(sites, {})), \
             patch("multiplex_design.mpx.dimer_matrix", return_value=rows):
            blocked, conflicts, _ = assess_candidates([p], Path("mock"), None, settings, lambda _: None)
        self.assertIn(edge_key("A", "A"), blocked)
        self.assertEqual(conflicts[0]["cross_products"], 1)

    def test_own_tm_spread_is_checked_against_pool_limit(self):
        p = candidate("A", "T1")
        p["design"]["right"]["tm_c"] = 64
        seqs = sorted({p["F"], p["R"]})
        rows = [{"primer_a": a, "primer_b": b, "dG37_any": 0, "dG37_3prime": 0}
                for a, b in itertools.combinations_with_replacement(seqs, 2)]
        settings = {"min_product": 50, "max_product": 3000, "dimer_any": -6,
                    "dimer_3p": -5, "tm_spread": 3, "product_gap": 0}
        with patch("multiplex_design.mpx.dimer_matrix", return_value=rows):
            blocked, conflicts, _ = assess_candidates([p], None, None, settings, lambda _: None)
        self.assertIn(edge_key("A", "A"), blocked)
        self.assertIn("Çiftin kendi Tm aralığı", conflicts[0]["reasons"])

    def test_third_candidate_cannot_hide_cross_product_of_selected_candidates(self):
        a = candidate("A", "T1", "A" * 20, "T" * 20)
        b = candidate("B", "T2", "C" * 20, "G" * 20)
        other = candidate("C", "T3", a["F"], b["R"])
        sites = [
            BindingSite(a["F"], a["F"], "chr1", 100, 120, "+", 0),
            BindingSite(b["R"], b["R"], "chr1", 160, 180, "-", 0),
        ]
        settings = {"min_product": 50, "max_product": 3000, "dimer_any": -6,
                    "dimer_3p": -5, "tm_spread": 3, "product_gap": 0}
        seqs = {p[role] for p in [a, b, other] for role in ["F", "R"]}
        dimers = [{"primer_a": x, "primer_b": y, "dG37_any": 0, "dG37_3prime": 0}
                  for x, y in itertools.combinations_with_replacement(sorted(seqs), 2)]
        with patch("multiplex_design.mpx.load_bindings", return_value=(sites, {})), \
             patch("multiplex_design.mpx.dimer_matrix", return_value=dimers):
            blocked, conflicts, _ = assess_candidates([a, b, other], Path("mock"), None, settings, lambda _: None)
        self.assertIn(edge_key("A", "B"), blocked)
        self.assertEqual(next(c["cross_products"] for c in conflicts if (c["a"], c["b"]) == ("A", "B")), 1)


if __name__ == "__main__":
    unittest.main()
