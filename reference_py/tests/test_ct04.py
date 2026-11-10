"""CT-04 allele orientation, energy and report checks."""
import unittest

from oligoscout_thermo import concentration, engine, nn, params, profiles


class TestCT04(unittest.TestCase):
    primer = "CAGCACTGTGAACAGTATGTTTCT"

    def setUp(self):
        self.engine = engine.ThermoEngine(profile=profiles.get("ctdna").qc)

    def check(self, mismatch, strand="+"):
        if strand == "-":
            return self.engine.evaluate_ct04("a", self.primer, nn.reverse_complement(self.primer),
                                             nn.reverse_complement(mismatch), strand)
        return self.engine.evaluate_ct04("a", self.primer, self.primer, mismatch, strand)

    def test_plus_terminal_record(self):
        row = self.check(self.primer[:-1] + "A")
        self.assertEqual((row.id, row.value["position_from_3prime"]), ("CT-04", 1))
        self.assertEqual(row.value["mismatch_kind"], "terminal")
        self.assertIn("KALIBRE DEGIL", row.threshold["status"])
        self.assertEqual(row.param_set_sha256, params.registry_digest())

    def test_minus_strand_same_metrics(self):
        mismatch = self.primer[:-1] + "G"
        plus, minus = self.check(mismatch), self.check(mismatch, "-")
        self.assertEqual(minus.value["strand"], "-")
        for key in ("theta_matched", "theta_mismatched", "theta_margin", "delta_delta_G37_kcal_mol"):
            self.assertEqual(plus.value[key], minus.value[key])

    def test_penultimate_uses_internal_nn_and_preserves_spc_contract(self):
        mismatch = self.primer[:-2] + "A" + self.primer[-1]
        row = self.check(mismatch)
        self.assertEqual(row.value["position_from_3prime"], 2)
        self.assertEqual(row.model["mismatch_nn"], params.MISMATCH_ALLAWI_PEYRET.provenance.table_id)
        with self.assertRaises(ValueError):
            nn.duplex_with_mismatches(self.primer, mismatch)

    def test_terminal_uses_wc_core_and_source_increment(self):
        matched = nn.watson_crick(self.primer)
        core = nn.watson_crick(self.primer[:-1])
        # Primer ends CT and matched target ends CT in primer orientation.
        # Replacing the last target base with G produces CT/GC.
        mismatch = self.primer[:-1] + "G"
        expected = core.dG_kcal() - 3.9 + 310.15 * 10.6 / 1000.0 - matched.dG_kcal()
        self.assertAlmostEqual(self.check(mismatch).value["delta_delta_G37_kcal_mol"], expected)

    def test_margin_order_follows_weaker_allele_binding(self):
        mismatches = [self.primer[:-1] + b for b in "ACG"]
        ordered = sorted((self.check(a).value for a in mismatches), key=lambda row: row["theta_mismatched"])
        self.assertLess(ordered[0]["theta_mismatched"], ordered[-1]["theta_mismatched"])
        self.assertGreater(ordered[0]["theta_margin"], ordered[-1]["theta_margin"])
        for row in ordered:
            self.assertAlmostEqual(row["theta_margin"], row["theta_matched"] - row["theta_mismatched"])

    def test_unsupported_positions_and_multiple_mismatches_fail(self):
        for mismatch in ("A" + self.primer[1:], self.primer[:-2] + "AG", self.primer):
            with self.subTest(mismatch=mismatch), self.assertRaises(ValueError):
                self.check(mismatch)

    def test_terminal_table_covers_all_48_contexts(self):
        table = params.TERMINAL_MISMATCH_SANTALUCIA_PEYRET.require()
        for previous in "ACGT":
            for terminal in "ACGT":
                p = "CAGCT" + previous + terminal
                for alternate in set("ACGT") - {terminal}:
                    model = concentration.build(p, template_site=p[:-1] + alternate, arms_mismatch=True)
                    self.assertGreater(model.tm_K, 0)
        self.assertEqual(len(table), 48)
