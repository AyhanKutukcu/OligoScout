"""LNA modeli (McTigue, Peterson & Kahn 2004) testleri."""
from __future__ import annotations

import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import concentration, engine, lna, params, salt, seqpair  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestLnaTable(unittest.TestCase):
    def test_32_parameters(self):
        t = params.LNA_MCTIGUE_2004.require()
        self.assertEqual(len(t["lna_5prime"]), 16)
        self.assertEqual(len(t["lna_3prime"]), 16)

    def test_spot_values_from_table4(self):
        t = params.LNA_MCTIGUE_2004.require()
        self.assertEqual(t["lna_5prime"]["GG"], [-2.844, -6.680, -0.700])   # 5'-G^L G
        self.assertEqual(t["lna_3prime"]["AC"], [2.890, 10.576, -0.390])    # 5'-A C^L
        self.assertEqual(t["lna_3prime"]["CT"], [-1.671, -4.070, -0.409])   # 5'-C T^L

    def test_ddg_close_to_ddh_minus_t_dds(self):
        # ddG37 SVD ile ayri kestirildi; makale yalniz yaklasik esitlik bildirir
        t = params.LNA_MCTIGUE_2004.require()
        for part in ("lna_5prime", "lna_3prime"):
            for h, s, g in t[part].values():
                self.assertLess(abs(h - 310.15 * s / 1000 - g), 0.15)


class TestPaperReproduction(unittest.TestCase):
    """Tablo 1: 100 LNA oligosu; makale 32NN combined icin RMSD(dTm)=0.88 C."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, os.path.join(ROOT, "reference_py", "cli"))
        import lna_validate
        cls.summary = lna_validate.summarize(lna_validate.run())

    def test_measured_reference_plus_table_matches_paper(self):
        m3 = self.summary["M3_measured_reference_plus_table"]
        self.assertEqual(self.summary["n"], 98)
        self.assertAlmostEqual(m3["rmsd_dtm_c"], 0.88, delta=0.05)

    def test_full_oligoscout_model(self):
        m4 = self.summary["M4_oligoscout_unified_dna_plus_table"]
        self.assertLess(m4["rmsd_tm_c"], 1.8)
        self.assertLess(m4["rmsd_dtm_c"], 1.0)


class TestDomain(unittest.TestCase):
    S = "ACGTTGCAAGGTCCATTGCA"

    def test_single_internal_in_domain(self):
        inc = lna.increments(self.S, [8])
        self.assertTrue(inc.in_domain)
        self.assertEqual(len(inc.per_position[0]["terms"]), 2)

    def test_terminal_adjacent_and_multiple_flagged(self):
        self.assertFalse(lna.increments(self.S, [1]).in_domain)
        inc = lna.increments(self.S, [8, 9])
        self.assertFalse(inc.in_domain)
        self.assertTrue(any("bitisik" in r for r in inc.out_of_domain))

    def test_mismatch_next_to_lna_skips_term(self):
        site = self.S[:9] + ("A" if self.S[9] != "A" else "C") + self.S[10:]
        inc = lna.increments(self.S, [8], template_site=site)
        applied = [t["applied"] for t in inc.per_position[0]["terms"]]
        self.assertEqual(applied, [True, False])
        self.assertFalse(inc.in_domain)

    def test_lna_raises_tm_for_stabilising_context(self):
        c = salt.Conditions()
        base = concentration.build(self.S, c).tm_c
        with_lna = concentration.build(self.S, c, lna_positions=[10]).tm_c  # 5'-G G^L T
        self.assertGreater(with_lna, base)


class TestIntegration(unittest.TestCase):
    def test_engine_accepts_plus_notation(self):
        e = engine.ThermoEngine()
        a = e.evaluate("p", "ACGTTGCAAGGTCCATTGCA")
        b = e.evaluate("p", "ACGTTGCAAG+GTCCATTGCA")
        self.assertEqual(b["sequence"], a["sequence"])
        self.assertIn("lna", b["thermo"])
        self.assertEqual(b["modifications"]["lna_positions"], [10])
        self.assertNotEqual(a["thermo"]["tm_c"], b["thermo"]["tm_c"])
        self.assertIn("LNA-SS", b["not_evaluated"])

    def test_pair_with_lna(self):
        r = engine.ThermoEngine().evaluate_pair("x", "ACGTTGCAAG+GTCCATTGCA", "TGCAATGGACCTTGCAACGT")
        self.assertTrue(r["checks"])

    def test_ct05_uses_lna(self):
        rec = seqpair.check_ct05("TTGCAAG+GTCCATTG", "TTGCAAGGTCCATTG")
        self.assertTrue(rec["value"]["lna_modelled"])
        self.assertTrue(rec["value"]["lna_in_domain"])


@unittest.skipUnless(shutil.which("node"), "node yok")
class TestBrowserPort(unittest.TestCase):
    def test_engine_js_matches_python(self):
        data = os.path.join(ROOT, "web", "oligoscout_data.json")
        if not os.path.isfile(data):
            self.skipTest("web veri paketi yok")
        with open(data, encoding="utf-8") as fh:
            if "lna" not in json.load(fh).get("engine_params", {}):
                self.skipTest("web veri paketi LNA tablosu olmadan uretilmis")
        rnd = random.Random(7)
        cases = []
        for _ in range(300):
            n = rnd.randint(10, 30)
            s = "".join(rnd.choice("ACGT") for _ in range(n))
            pos = sorted(rnd.sample(range(n), rnd.randint(1, 4)))
            c = salt.Conditions(monovalent_M=0.05, magnesium_M=0.003, dntp_M=0.0006, primer_M=2e-7)
            m = concentration.build(s, c, salt_model="owczarzy_2008", lna_positions=pos)
            cases.append({"s": s, "pos": pos, "tm": m.tm_c, "dom": m.lna["in_domain"]})
        js = ("const fs=require('fs');eval(fs.readFileSync(process.argv[2],'utf8')+';globalThis.OS=OS;');"
              "OS.setParams(JSON.parse(fs.readFileSync(process.argv[3],'utf8')).engine_params);"
              "const C=JSON.parse(fs.readFileSync(process.argv[4],'utf8'));let bad=0;"
              "for(const c of C){const t=OS.thermo(c.s,{na:50,mg:3,dntp:0.6,pnm:200,ta:60},'owczarzy_2008',c.pos);"
              "if(Math.abs(t.tm-c.tm)>1e-9||t.lna.inDomain!==c.dom)bad++;}console.log(bad);")
        with tempfile.TemporaryDirectory() as tmp:
            cp, sp = os.path.join(tmp, "c.json"), os.path.join(tmp, "t.js")
            with open(cp, "w") as fh:
                json.dump(cases, fh)
            with open(sp, "w") as fh:
                fh.write(js)
            out = subprocess.run(["node", sp, os.path.join(ROOT, "web", "src", "engine.js"), data, cp],
                                 capture_output=True, text=True, timeout=60)
        self.assertEqual(out.stdout.strip(), "0", out.stderr)


if __name__ == "__main__":
    unittest.main()
