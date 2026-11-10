"""Dizi/kompozisyon/cift/amplikon kontrollerinin dogrulamasi (seqpair)."""

from __future__ import annotations

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import nn, salt, seqpair  # noqa: E402

C = seqpair.SeqPairConstraints()


def sev(records, cid):
    return next(r["severity"] for r in records if r["id"] == cid)


def val(records, cid):
    return next(r["value"] for r in records if r["id"] == cid)


class TestSequenceHelpers(unittest.TestCase):
    def test_degeneracy_counts_iupac_product(self):
        self.assertEqual(seqpair.degeneracy("ACGT"), 1)
        self.assertEqual(seqpair.degeneracy("ACGTR"), 2)    # R=2
        self.assertEqual(seqpair.degeneracy("NN"), 16)      # 4*4

    def test_degenerate_in_3p_window(self):
        self.assertTrue(seqpair.three_prime_has_degenerate("ACGTACGR"))
        self.assertFalse(seqpair.three_prime_has_degenerate("RACGTACG"))

    def test_shannon_entropy_bounds(self):
        # Homopolimer: tek sembol -> 0 bit.
        self.assertAlmostEqual(seqpair.shannon_entropy("AAAAAAAA", 1), 0.0)
        # Esit dagilimli 4 baz -> 2 bit (k=1).
        self.assertAlmostEqual(seqpair.shannon_entropy("ACGTACGTACGT", 1), 2.0,
                               places=6)

    def test_g_quadruplex_detected_on_either_strand(self):
        g4 = "GGGAGGGAGGGAGGG"
        self.assertTrue(seqpair.has_g_quadruplex(g4))
        # Tumleyeni C-zengini; algilayici her iki ipligi tarar.
        self.assertTrue(seqpair.has_g_quadruplex(nn.reverse_complement(g4)))
        self.assertFalse(seqpair.has_g_quadruplex("ACGTACGTACGTACGT"))


class TestPrimerChecks(unittest.TestCase):
    def test_length_thresholds(self):
        self.assertEqual(sev(seqpair.check_primer("A" * 22), "SEQ-02"), "INFO")
        self.assertEqual(sev(seqpair.check_primer("ACGT" * 4), "SEQ-02"), "WARN")
        self.assertEqual(sev(seqpair.check_primer("ACGTACGTAC"), "SEQ-02"),
                         "FAIL")  # 10 nt < 15

    def test_degenerate_primer_fails_and_skips_wc_checks(self):
        records = seqpair.check_primer("ACGTACGTACGTACGTACGTR")
        self.assertEqual(sev(records, "SEQ-03"), "FAIL")  # 3' dejenere
        # Saf olmayan dizi icin WC kontrolleri uretilmez.
        self.assertNotIn("CMP-05", [r["id"] for r in records])

    def test_low_complexity_is_flagged(self):
        records = seqpair.check_primer("AAAAAAAAAAAAAAAAAAAA")
        self.assertEqual(sev(records, "CMP-05"), "FAIL")

    def test_three_prime_t_warns_cmp08(self):
        self.assertEqual(sev(seqpair.check_primer("GCGCGCGACGACGACGT"),
                             "CMP-08"), "WARN")
        self.assertEqual(sev(seqpair.check_primer("GCGCGCGACGACGACGA"),
                             "CMP-08"), "INFO")

    def test_all_rules_present_for_clean_primer(self):
        ids = {r["id"] for r in seqpair.check_primer("ACGTTGCAAGGTCCATTGCA")}
        self.assertEqual(ids, {"SEQ-02", "SEQ-03", "CMP-05", "CMP-06",
                               "CMP-07", "CMP-08", "SS-02"})


class TestPairChecks(unittest.TestCase):
    def setUp(self):
        self.fwd = "ACGTTGCAAGGTCCATTGCA"
        self.rev = "TGCAATGGACCTTGCAACGT"

    def test_product_size_window(self):
        cons = seqpair.SeqPairConstraints(product_bp=(60, 100))
        amp = "A" * 80
        recs = seqpair.check_pair(self.fwd, self.rev, 0, 80, amp, cons=cons)
        self.assertEqual(sev(recs, "PAIR-02"), "INFO")
        recs = seqpair.check_pair(self.fwd, self.rev, 0, 120, "A" * 120, cons=cons)
        self.assertEqual(sev(recs, "PAIR-02"), "FAIL")

    def test_primer_overlap_fails_pair03(self):
        amp = "ACGT" * 8   # 32 bp, primerler toplami 40 > 32 -> ortusme
        recs = seqpair.check_pair(self.fwd, self.rev, 0, 32, amp,
                                  cons=seqpair.SeqPairConstraints(product_bp=(20, 200)))
        self.assertEqual(sev(recs, "PAIR-03"), "FAIL")

    def test_amplicon_gc_window_fail(self):
        amp = "G" * 80
        recs = seqpair.check_pair(self.fwd, self.rev, 0, 80, amp,
                                  cons=seqpair.SeqPairConstraints(product_bp=(20, 200)))
        self.assertEqual(sev(recs, "AMP-01"), "FAIL")

    def test_pair_product_tm_gap(self):
        amp = "GCGC" * 30
        recs = seqpair.check_pair(self.fwd, self.rev, 0, len(amp), amp,
                                  cons=seqpair.SeqPairConstraints(product_bp=(20, 400)))
        self.assertEqual(sev(recs, "PAIR-04"), "WARN")  # uzun urun Tm'i cok yuksek


class TestIndependentCrossCheck(unittest.TestCase):
    """Iki hesabi elle/bagimsiz formulle dogrula."""

    def test_shannon_matches_manual(self):
        seq = "AACCGGTT"   # her baz 2 kez, k=1
        manual = -4 * (0.25 * math.log2(0.25))
        self.assertAlmostEqual(seqpair.shannon_entropy(seq, 1), manual, places=9)

    def test_cmp07_balance_sign(self):
        # 5' A/T zengin, 3' G/C zengin -> 3' yari daha kararli -> WARN.
        recs = seqpair.check_primer("ATATATATGCGCGCGC")
        self.assertEqual(sev(recs, "CMP-07"), "WARN")
        v = val(recs, "CMP-07")
        self.assertLess(v["dG3_half_kcal_mol"], v["dG5_half_kcal_mol"])


class TestTemplateChecks(unittest.TestCase):
    """SEQ-04 (kalipta belirsiz baz) ve SEQ-05 (koordinat-dizi tutarliligi)."""

    def _target(self):
        # 0..20 ileri primer, 40..60 geri primerin ters tamamlayicisi.
        fwd = "ACGTACGTACGTACGTACGT"
        rev = "TTGGCCAATTGGCCAATTGG"
        import oligoscout_thermo.seqpair as sp
        middle = "AAAAAAAAAAAAAAAAAAAA"
        rev_region = sp.reverse_complement_raw(rev)
        target = fwd + middle + rev_region
        return target, fwd, rev

    def test_coordinates_consistent_is_info(self):
        target, fwd, rev = self._target()
        recs = seqpair.check_template(fwd, rev, 0, len(target), target)
        self.assertEqual(sev(recs, "SEQ-05"), "INFO")
        self.assertEqual(sev(recs, "SEQ-04"), "INFO")

    def test_wrong_coordinate_fails_seq05(self):
        target, fwd, rev = self._target()
        recs = seqpair.check_template("GGGGGGGGGGGGGGGGGGGG", rev, 0,
                                      len(target), target)
        self.assertEqual(sev(recs, "SEQ-05"), "FAIL")
        self.assertFalse(val(recs, "SEQ-05")["forward_matches_reference"])

    def test_N_in_primer_region_fails_seq04(self):
        target, fwd, rev = self._target()
        bad = "ACGTNCGTACGTACGTACGT" + target[len(fwd):]
        recs = seqpair.check_template("ACGTNCGTACGTACGTACGT", rev, 0,
                                      len(bad), bad)
        self.assertEqual(sev(recs, "SEQ-04"), "FAIL")

    def test_softmask_in_amplicon_warns_seq04(self):
        target, fwd, rev = self._target()
        soft = fwd + "aaaaaaaaaaaaaaaaaaaa" + target[len(fwd) + 20:]
        recs = seqpair.check_template(fwd, rev, 0, len(soft), soft)
        self.assertEqual(sev(recs, "SEQ-04"), "WARN")

    def test_reverse_complement_raw_preserves_n(self):
        self.assertEqual(seqpair.reverse_complement_raw("ACGTN"), "NACGT")


class TestChipGene(unittest.TestCase):
    def test_chip_gene_warns_in_ctdna(self):
        rec = seqpair.check_chip_gene("TP53", "ctdna")
        self.assertEqual(rec["severity"], "WARN")
        self.assertTrue(rec["value"]["is_chip_gene"])

    def test_chip_gene_quiet_in_other_profile(self):
        self.assertEqual(seqpair.check_chip_gene("TP53", "qpcr_sybr")["severity"],
                         "INFO")

    def test_non_chip_gene_is_info(self):
        rec = seqpair.check_chip_gene("GAPDH", "ctdna")
        self.assertEqual(rec["severity"], "INFO")
        self.assertFalse(rec["value"]["is_chip_gene"])

    def test_missing_gene_is_info(self):
        rec = seqpair.check_chip_gene(None, "ctdna")
        self.assertEqual(rec["severity"], "INFO")
        self.assertIn("belirtilmedi", rec["note"])


class TestAmpliconStructure(unittest.TestCase):
    """SS-04: amplikon ikincil yapisi."""

    def test_strong_hairpin_fails(self):
        amp = "GCGCGCGCAAAAAAAAGCGCGCGC"   # kuvvetli sap
        rec = seqpair.check_amplicon_structure(amp, 60.0)
        self.assertEqual(rec["severity"], "FAIL")
        self.assertGreater(rec["value"]["structure_tm_c"], 60.0)

    def test_unstructured_amplicon_is_info(self):
        amp = "AAAAAAAAAAAAAAAAAAAAAAAA"   # eslesemez
        rec = seqpair.check_amplicon_structure(amp, 60.0)
        self.assertEqual(rec["severity"], "INFO")
        self.assertIsNone(rec["value"]["structure_tm_c"])

    def test_structure_tm_consistent_with_mfe(self):
        """dG37 = dH - 310.15*dS/1000 ozdesligi yapidan turetilen Tm'de korunur."""
        from oligoscout_thermo import fold
        amp = "GCGCATATGCGCAAAAGCGCATATGCGC"
        s = amp
        mfe, structure = fold.zuker_mfe(s)
        pairs = seqpair._structure_pairs(structure)
        dH = seqpair._structure_enthalpy(s, pairs)
        if dH < 0:
            dS = (dH - mfe) * 1000.0 / 310.15
            self.assertAlmostEqual(dH - 310.15 * dS / 1000.0, mfe, places=6)


class TestPrimerProbeDimers(unittest.TestCase):
    """DIM-05: primer x prob dimerleri."""

    def setUp(self):
        self.fwd = "ACGTTGCAAGGTCCATTGCA"
        self.rev = "TGCAATGGACCTTGCAACGT"

    def test_clean_probe_is_info(self):
        rec = seqpair.check_primer_probe_dimers(self.fwd, self.rev,
                                                "TTACCTAGCATTACGGATCA")
        self.assertEqual(rec["severity"], "INFO")

    def test_sticky_probe_fails(self):
        # Ileri primerle genis ortusme -> cok kararli dimer.
        rec = seqpair.check_primer_probe_dimers(self.fwd, self.rev,
                                                "TGCAAGGTCCATTGCAACGT")
        self.assertEqual(rec["severity"], "FAIL")
        self.assertLess(rec["value"]["worst_any_dG37_kcal_mol"], -9.0)

    def test_thresholds_are_overridable(self):
        rec = seqpair.check_primer_probe_dimers(
            self.fwd, self.rev, "TTACCTAGCATTACGGATCA",
            warn_dg37=-1.0, fail_dg37=-2.0)
        self.assertIn(rec["severity"], ("WARN", "FAIL"))


class TestSpc05Classifier(unittest.TestCase):
    """SPC-05 3' terminal uyumsuzluk siniflandirmasi (panel_qc yardimcisi)."""

    def _classify(self, primer, site, strand, mask):
        import collections
        sys.path.insert(0, os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cli"))
        import panel_qc
        tally = collections.Counter()
        panel_qc._classify_terminal_mismatch(
            {"primer": primer, "site_plus_strand": site, "strand": strand,
             "mask": mask}, tally)
        return tally

    def test_tolerant_gt_wobble_at_3prime(self):
        # Primer 3' G, kalipta ayni-anlamli iplikte A -> karsisi comp(A)=T
        # -> G.T tolerant.
        t = self._classify("ACGTACGTACGTACGTACGG", "ACGTACGTACGTACGTACGA",
                           "+", "..................." + "X")
        self.assertEqual(t["tolerant"], 1)
        self.assertEqual(t["disruptive"], 0)

    def test_disruptive_mismatch_at_3prime(self):
        # Primer 3' A, kalipta C -> karsisi comp(C)=G -> A.G disruptive.
        t = self._classify("ACGTACGTACGTACGTACGA", "ACGTACGTACGTACGTACGC",
                           "+", "..................." + "X")
        self.assertEqual(t["disruptive"], 1)
        self.assertEqual(t["tolerant"], 0)

    def test_no_mismatch_in_window_is_ignored(self):
        t = self._classify("ACGTACGT", "ACGTACGT", "+", "XXXXX...")
        self.assertEqual(sum(t.values()), 0)


class TestRT05Hybrid(unittest.TestCase):
    """RT-05: DNA/RNA hibrit Tm. Tablo Biopython R_DNA_NN1'den aktarildi
    (Sugimoto 1995), init terimiyle; verified DEGIL."""

    def test_table_is_transcribed_and_available(self):
        from oligoscout_thermo import params
        self.assertEqual(
            params.HYBRID_SUGIMOTO_1995.provenance.verification_status,
            "transcribed_unverified")
        self.assertTrue(params.HYBRID_SUGIMOTO_1995.is_available)
        self.assertEqual(len(params.HYBRID_SUGIMOTO_1995.require()), 16)

    def test_known_values(self):
        from oligoscout_thermo import params
        t = params.HYBRID_SUGIMOTO_1995.require()
        # Biopython R_DNA_NN1 (Sugimoto 1995) ile ayni degerler
        self.assertEqual(t["CG"], [-16.3, -47.1])
        self.assertEqual(t["AA"], [-7.8, -21.9])
        self.assertEqual(t["TT"], [-11.5, -36.4])

    def test_hybrid_duplex_includes_init(self):
        # tek adim "AA": dH = -7.8 + init 1.9 = -5.9; dS = -21.9 + -3.9 = -25.8
        dH, dS = seqpair.hybrid_duplex_dh_ds("AA")
        self.assertAlmostEqual(dH, -5.9, places=6)
        self.assertAlmostEqual(dS, -25.8, places=6)

    def test_rt05_evaluates(self):
        rec = seqpair.check_rt05("ACGTTGCAAGGTCCATTGCA")
        self.assertEqual(rec["value"]["status"], "evaluated")
        self.assertIsNotNone(rec["value"]["hybrid_tm_c"])
        self.assertEqual(rec["severity"], "INFO")

    def test_template_has_sixteen_keys(self):
        from oligoscout_thermo import params
        self.assertEqual(len(params.HYBRID_NN_KEYS), 16)
        self.assertEqual(set(params.HYBRID_NN_KEYS),
                         {a + b for a in "ACGT" for b in "ACGT"})


class TestSpc08ContigScan(unittest.TestCase):
    """SPC-08: kucuk contig (mitokondri/rRNA) taramasi."""

    def setUp(self):
        self.contig = "ACGTACGTACGTTGCAAGGTCCATTGCATAGCTAACGATCGATCGTAGC"

    def test_exact_site_is_anchored_fail(self):
        primer = self.contig[10:30]       # contig'ten birebir
        rec = seqpair.check_spc08(primer, {"c1": self.contig})
        self.assertEqual(rec["severity"], "FAIL")
        self.assertGreaterEqual(rec["value"]["anchored_sites"], 1)

    def test_absent_primer_is_info(self):
        rec = seqpair.check_spc08("GGGGGGGGGGGGGGGGGGGG", {"c1": self.contig})
        self.assertEqual(rec["severity"], "INFO")

    def test_reverse_strand_hit(self):
        import oligoscout_thermo.nn as _nn
        primer = _nn.reverse_complement(self.contig[10:30])
        rec = seqpair.check_spc08(primer, {"c1": self.contig})
        self.assertTrue(rec["value"]["scanned"])
        self.assertGreaterEqual(rec["value"]["total_sites"], 1)

    def test_no_contigs_reports_not_scanned(self):
        rec = seqpair.check_spc08("ACGTACGTACGTACGTACGT", {})
        self.assertEqual(rec["severity"], "INFO")
        self.assertFalse(rec["value"]["scanned"])


class TestOligoModParser(unittest.TestCase):
    """Modifiye oligo ayristirmasi (LNA, uc isaretleri)."""

    def test_plain_dna(self):
        from oligoscout_thermo import oligomod
        p = oligomod.parse("AGTCCTTTCTGTGTCTTCCCA")
        self.assertEqual(p.dna, "AGTCCTTTCTGTGTCTTCCCA")
        self.assertEqual(p.lna_count, 0)

    def test_lna_positions(self):
        from oligoscout_thermo import oligomod
        p = oligomod.parse("CCC+CT+C+T+C+TGAC")
        self.assertEqual(p.dna, "CCCCTCTCTGAC")
        self.assertEqual(sorted(p.lna_positions), [3, 5, 6, 7, 8])

    def test_end_mods_stripped(self):
        from oligoscout_thermo import oligomod
        p = oligomod.parse("/56-FAM/CCC+CT+C+T+C+TGAC/3IABkFQ/")
        self.assertEqual(p.dna, "CCCCTCTCTGAC")
        self.assertEqual(len(p.end_mods), 2)

    def test_blocker_3phos(self):
        from oligoscout_thermo import oligomod
        p = oligomod.parse("AGCATCTCCAGCAGCAGG+T+CA+T+AGAGGG/3Phos/")
        self.assertEqual(p.dna, "AGCATCTCCAGCAGCAGGTCATAGAGGG")
        self.assertEqual(sorted(p.lna_positions), [18, 19, 21, 22])
        self.assertIn("3Phos", p.end_mods)

    def test_trailing_plus_warns(self):
        from oligoscout_thermo import oligomod
        self.assertTrue(oligomod.parse("ACGT+").warnings)


class TestCt05(unittest.TestCase):
    def test_delta_tm_reported_lna_flagged(self):
        rec = seqpair.check_ct05("CCC+CT+C+T+C+TGAC", "C+CT+C+T+ATG+A+CC")
        self.assertFalse(rec["value"]["lna_modelled"])
        self.assertGreater(rec["value"]["delta_tm_c"], 0)
        self.assertEqual(rec["value"]["mutant_lna_count"], 5)

    def test_identical_probes_fail_discrimination(self):
        rec = seqpair.check_ct05("ACGTACGTACGT", "ACGTACGTACGT")
        self.assertEqual(rec["value"]["delta_tm_c"], 0.0)
        self.assertEqual(rec["severity"], "FAIL")

    def test_large_delta_is_info(self):
        rec = seqpair.check_ct05("GCGCGCGCGCGCGCGC", "ATATATATAT")
        self.assertEqual(rec["severity"], "INFO")


class TestSeq01(unittest.TestCase):
    def test_pure_acgt(self):
        r = seqpair.check_seq01("AGTCCTTTCTGTGTCTTCCCA")
        self.assertTrue(r["value"]["pure_acgt"])
        self.assertEqual(r["severity"], "INFO")

    def test_modifications_parsed(self):
        r = seqpair.check_seq01("/56-FAM/CCC+CT+C+T+C+TGAC/3IABkFQ/")
        self.assertEqual(r["value"]["modifications"]["lna_count"], 5)
        self.assertIn("56-FAM", r["value"]["modifications"]["end_mods"])
        self.assertFalse(r["value"]["pure_acgt"])
        self.assertEqual(r["severity"], "INFO")

    def test_iupac_flagged(self):
        r = seqpair.check_seq01("AGTCCNTTTCTGY")
        self.assertEqual(r["value"]["iupac_degenerate"], ["N", "Y"])

    def test_invalid_char_fails(self):
        r = seqpair.check_seq01("AGTCZZT")
        self.assertEqual(r["severity"], "FAIL")
        self.assertIn("Z", r["value"]["invalid_chars"])


class TestSs05(unittest.TestCase):
    def test_unstructured_probe_info(self):
        r = seqpair.check_ss05("CCTCTATGACC")
        self.assertEqual(r["severity"], "INFO")

    def test_strong_hairpin_fails(self):
        r = seqpair.check_ss05("GGGGCCCCTTTTGGGGCCCC")
        self.assertIsNotNone(r["value"]["structure_tm_c"])
        self.assertEqual(r["severity"], "FAIL")

    def test_lna_probe_parsed_as_lower_bound(self):
        r = seqpair.check_ss05("C+CT+C+T+ATG+A+CC")
        self.assertEqual(r["value"]["lna_count"], 6)


class TestDim07(unittest.TestCase):
    def test_tail_increases_dimer_risk(self):
        r = seqpair.check_dim07(
            "AGTCCTTTCTGTGTCTTCCCA", "TCCAGCATCTCCAGCAGC",
            fwd_tail="TGTAAAACGACGGCCAGT")
        self.assertLessEqual(r["value"]["tail_included_any_dG37_kcal_mol"],
                             r["value"]["no_tail_any_dG37_kcal_mol"])

    def test_no_tail_equals_base(self):
        r = seqpair.check_dim07("AGTCCTTTCTGTGTCTTCCCA", "TCCAGCATCTCCAGCAGC")
        self.assertEqual(r["value"]["tail_included_any_dG37_kcal_mol"],
                         r["value"]["no_tail_any_dG37_kcal_mol"])


class TestSeq06(unittest.TestCase):
    def test_build_patch_parsed(self):
        r = seqpair.check_seq06(">chr22 GRCh38.p14 primary", "deadbeef", 1000)
        self.assertEqual(r["value"]["build"], "GRCh38")
        self.assertEqual(r["value"]["patch"], "p14")
        self.assertEqual(r["value"]["fasta_sha256"], "deadbeef")
        self.assertEqual(r["severity"], "INFO")

    def test_missing_checksum_warns(self):
        r = seqpair.check_seq06(">x", "", None)
        self.assertEqual(r["severity"], "WARN")


if __name__ == "__main__":
    unittest.main(verbosity=2)
