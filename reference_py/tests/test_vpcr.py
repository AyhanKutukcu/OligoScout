"""Sanal PCR (vpcr): dejenere oligolar, dairesel kalip, kuyruklar, prob, baglanti, Ta."""
from __future__ import annotations

import math
import os
import random
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import salt, vpcr  # noqa: E402

RC = vpcr.revcomp_iupac


def rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


class Plasmid(unittest.TestCase):
    """3000 bp halka; F 400'de (+), R 900'de (-); bir de orijinden gecen cift."""

    @classmethod
    def setUpClass(cls):
        s = list(rnd(3000, 11))
        cls.F = "GACTGCATTGCAGTCAAGCT"
        cls.R = "TTGGCCAAGTCGATCGGAAC"
        s[400:420] = cls.F
        s[880:900] = RC(cls.R)
        # orijinden gecen urun: F2 2900'de, R2 halkanin basinda (100)
        cls.F2 = "CATGCGATAGCTTAGCCGAT"
        cls.R2 = "AGGTCCTTACGGATTCGACA"
        s[2900:2920] = cls.F2
        s[80:100] = RC(cls.R2)
        cls.seq = "".join(s)
        cls.rec = vpcr.Record("pX", cls.seq, circular=True)
        cls.records = {"pX": (3000, True)}

    def scan(self, oligos, rule=vpcr.ScanRule(), circular=True):
        rec = vpcr.Record("pX", self.seq, circular)
        return vpcr.scan_python([rec], oligos, rule)

    def test_linear_product_and_sequence(self):
        ol = [vpcr.parse_oligo("F", self.F, pair="p1"), vpcr.parse_oligo("R", self.R, pair="p1")]
        om = {o.id: o for o in ol}
        prods = vpcr.assemble(self.scan(ol), om, self.records, 50, 2000, "pairs")
        fr = [p for p in prods if p.orientation == "F|R"]
        self.assertEqual(len(fr), 1)
        self.assertEqual((fr[0].start, fr[0].end, fr[0].length), (400, 900, 500))
        seqs = vpcr.amplicon_sequence(fr[0], self.rec, om)
        self.assertEqual(seqs["amplicon_sequence"], self.seq[400:900])

    def test_product_across_origin_only_on_circular(self):
        ol = [vpcr.parse_oligo("F2", self.F2, pair="p2"), vpcr.parse_oligo("R2", self.R2, pair="p2")]
        om = {o.id: o for o in ol}
        circ = vpcr.assemble(self.scan(ol), om, self.records, 50, 2000, "pairs")
        self.assertEqual([(p.start, p.length, p.wraps) for p in circ], [(2900, 200, True)])
        lin = vpcr.assemble(self.scan(ol, circular=False), om, {"pX": (3000, False)}, 50, 2000, "pairs")
        self.assertEqual(lin, [])
        seq = vpcr.amplicon_sequence(circ[0], self.rec, om)["template_sequence"]
        self.assertEqual(seq, self.seq[2900:] + self.seq[:100])

    def test_inverse_pcr(self):
        # disari bakan primerler: R'nin tamamlayicisi ileri, F'nin tamamlayicisi geri
        ol = [vpcr.parse_oligo("iF", RC(self.R), pair="inv"), vpcr.parse_oligo("iR", RC(self.F), pair="inv")]
        om = {o.id: o for o in ol}
        prods = vpcr.assemble(self.scan(ol), om, self.records, 50, 5000, "pairs")
        self.assertEqual([(p.start, p.length, p.wraps) for p in prods], [(880, 3000 - 480 + 20 + 20 - 20 - 20 + 0, True)]
                         if False else [(880, 2540, True)])

    def test_five_prime_tails_enter_amplicon(self):
        ol = [vpcr.parse_oligo("F", "gtaaaacgacggccagt" + self.F, pair="t"),
              vpcr.parse_oligo("R", "caggaaacagctatgac" + self.R, pair="t")]
        om = {o.id: o for o in ol}
        p = [x for x in vpcr.assemble(self.scan(ol), om, self.records, 50, 2000, "pairs") if x.orientation == "F|R"][0]
        self.assertEqual(p.template_length, 500)
        self.assertEqual(p.length, 534)
        amp = vpcr.amplicon_sequence(p, self.rec, om)["amplicon_sequence"]
        self.assertTrue(amp.startswith("GTAAAACGACGGCCAGT" + self.F))
        self.assertTrue(amp.endswith(RC("CAGGAAACAGCTATGAC")))

    def test_three_prime_tail_blocks_extension(self):
        ol = [vpcr.parse_oligo("F", self.F + "ttt", pair="b"), vpcr.parse_oligo("R", self.R, pair="b")]
        om = {o.id: o for o in ol}
        self.assertEqual(ol[0].tail3, "TTT")
        self.assertFalse([p for p in vpcr.assemble(self.scan(ol), om, self.records, 50, 2000, "pairs")
                          if "F" in p.orientation.split("|")])

    def test_probe_inside_and_link(self):
        probe = self.seq[600:625]
        ol = [vpcr.parse_oligo("F", self.F, pair="p1"), vpcr.parse_oligo("R", self.R, pair="p1"),
              vpcr.parse_oligo("P", probe, role="probe")]
        om = {o.id: o for o in ol}
        sites = self.scan(ol)
        prod = [p for p in vpcr.assemble(sites, om, self.records, 50, 2000, "pairs") if p.orientation == "F|R"][0]
        self.assertEqual(vpcr.probes_inside(prod, [s for s in sites if s.oligo == "P"]), ["P+@600"])
        links = vpcr.find_links(sites, "F", "P", 100, 300)
        self.assertEqual([(k.gap, k.insert_start, k.insert_end) for k in links], [(180, 420, 600)])


class Degenerate(unittest.TestCase):
    def test_iupac_primer_matches_any_member(self):
        rec = vpcr.Record("t", rnd(200, 5)[:50] + "CAGATGGACGTCAAGAC" + rnd(200, 6)[:50])
        o = vpcr.parse_oligo("RTp", "CARATGGAYGTNAARAC")
        sites = vpcr.scan_python([rec], [o], vpcr.ScanRule(0, 7, 0))
        self.assertEqual([(s.start, s.strand, s.mismatches) for s in sites], [(50, "+", 0)])
        self.assertEqual(o.degeneracy, 32)

    def test_template_iupac(self):
        rec = vpcr.Record("t", "A" * 20 + "ACGTRCGTACGTAC" + "A" * 20)
        o = vpcr.parse_oligo("o", "ACGTGCGTACGTAC")
        self.assertEqual(len(vpcr.scan_python([rec], [o], vpcr.ScanRule(0, 7, 0))), 1)
        recn = vpcr.Record("t", "A" * 20 + "ACGTNCGTACGTAC" + "A" * 20)
        self.assertEqual(vpcr.scan_python([recn], [o], vpcr.ScanRule(1, 7, 1)), [])
        self.assertEqual(len(vpcr.scan_python([recn], [o], vpcr.ScanRule(1, 7, 1, "mismatch"))), 1)

    def test_g_tolerance(self):
        o = vpcr.parse_oligo("o", "ACGTACGTACGTACGTAC")
        site = list("ACGTACGTACGTACGTAC")
        site[2] = "A"        # oligo G vs template-sense A: oligo G pairs with T -> G.T
        rec = vpcr.Record("t", "C" * 10 + "".join(site) + "C" * 10)
        strict = vpcr.scan_python([rec], [o], vpcr.ScanRule(0, 7, 0))
        loose = [s for s in vpcr.scan_python([rec], [o], vpcr.ScanRule(0, 7, 0, g_tolerant=True, g_max=1))
                 if s.start == 10]
        self.assertEqual(strict, [])
        self.assertEqual([(s.mismatches, s.g_mismatches, s.mask[2]) for s in loose], [(1, 1, "g")])
        site[2] = "T"        # oligo G vs template-sense T: G pairs with A -> G.A (tolerated)
        site[5] = "T"        # oligo C vs template-sense T: C pairs with A -> C.A (counted)
        rec = vpcr.Record("t", "C" * 10 + "".join(site) + "C" * 10)
        loose = [s for s in vpcr.scan_python([rec], [o], vpcr.ScanRule(1, 7, 1, g_tolerant=True, g_max=1))
                 if s.start == 10]
        self.assertEqual([(s.mismatches, s.g_mismatches) for s in loose], [(2, 1)])
        # g_max = 0: the G.A pair counts as an ordinary mismatch -> 2 > 1
        none = [s for s in vpcr.scan_python([rec], [o], vpcr.ScanRule(1, 7, 1, g_tolerant=True, g_max=0))
                if s.start == 10]
        self.assertEqual(none, [])


class FingerprintAndTa(unittest.TestCase):
    def test_single_primer_inverted_repeat(self):
        ltr = "TGGCAACGTTGGCTAGGCAT"
        seq = rnd(100, 1) + ltr + rnd(300, 2) + RC(ltr) + rnd(100, 3)
        o = vpcr.parse_oligo("LTR", ltr)
        sites = vpcr.scan_python([vpcr.Record("g", seq)], [o], vpcr.ScanRule())
        prods = vpcr.assemble(sites, {"LTR": o}, {"g": (len(seq), False)}, 50, 1000, "all")
        self.assertEqual([(p.orientation, p.length) for p in prods], [("LTR|LTR", 340)])

    def test_ta_formula(self):
        self.assertAlmostEqual(vpcr.annealing_temperature(60.0, 62.0, 500), round(60 + math.log(500), 2))
        self.assertIsNone(vpcr.annealing_temperature(None, None, 500))

    def test_site_tm_perfect_vs_mismatch(self):
        o = vpcr.parse_oligo("o", "ACGTTGCAAGGTCCATTGCA")
        perfect = vpcr.Site("o", "r", 0, 20, "+", 0, 0, "." * 20, "ACGTTGCAAGGTCCATTGCA")
        mm = vpcr.Site("o", "r", 0, 20, "+", 1, 0, "." * 9 + "X" + "." * 10, "ACGTTGCAATGTCCATTGCA")
        c = salt.Conditions()
        self.assertGreater(vpcr.site_tm(o, perfect, c)["tm_c"], vpcr.site_tm(o, mm, c)["tm_c"])


@unittest.skipUnless(vpcr.default_scan_binary(), "oligoscout_vpcr_scan derlenmemis")
class CppEqualsPython(unittest.TestCase):
    def test_random_templates(self):
        r = random.Random(9)
        recs = []
        for k in range(3):
            s = list(rnd(4000, 20 + k))
            for _ in range(8):
                pos = r.randrange(3900)
                site = "".join(r.choice(vpcr.IUPAC_SETS[c]) for c in "CARATGGAYGTNAARAC")
                site = list(site if r.random() < .5 else RC(site))
                for _ in range(r.randint(0, 3)):
                    site[r.randrange(len(site))] = r.choice("ACGT")
                s[pos:pos + 17] = site
            if k == 1:
                s[50:60] = "N" * 10
            recs.append(("r%d" % k, "".join(s)))
        ol = [vpcr.parse_oligo("a", "CARATGGAYGTNAARAC"), vpcr.parse_oligo("b", "CATRTCRTCNACRTA")]
        with tempfile.TemporaryDirectory() as tmp:
            fa = os.path.join(tmp, "t.fa")
            with open(fa, "w") as fh:
                fh.write("".join(f">{n}\n{q}\n" for n, q in recs))
            for rule in (vpcr.ScanRule(), vpcr.ScanRule(2, 5, 0, "mismatch", True, 1),
                         vpcr.ScanRule(1, 7, 1, "reject", True, 3)):
                for circ in ((), "all"):
                    py = vpcr.scan_python(vpcr.read_fasta(fa, circ), ol, rule)
                    cp, _ = vpcr.scan_cpp(fa, ol, rule, circ)
                    key = lambda s: (s.oligo, s.record, s.start, s.strand, s.mismatches, s.mask, s.wraps)
                    self.assertEqual(sorted(map(key, py)), sorted(map(key, cp)))


if __name__ == "__main__":
    unittest.main()
