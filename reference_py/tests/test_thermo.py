"""Termodinamik motorun dogrulama paketi (yalniz standart kutuphane).

Kapsam (QC spesifikasyonu, bolum 5 "Yazilimin kendi dogrulamasi"):
  * analitik ozdeslikler   theta(Tm) = 0.5, 1e-9 tolerans
  * tureyen isaretler      dTm/dlogC > 0 , dTm/d[Na+] > 0
  * degismezlikler         ters-tamamlayici, simetri duzeltmesi
  * toplanabilirlik        NN yigin toplami elle hesapla birebir
  * girdi sozlesmesi       alfabe disi karakter reddedilir (SEQ-01)
  * parametre kayit defteri checksum determinizmi, eksik tablo davranisi
  * bagimsiz oracle        Biopython varsa WC tablosuna karsi diferansiyel
                           test (yalniz gelistirme; urun bagimliligi degil)

Calistirma:  python3 -m unittest discover -s tests -v
"""

from __future__ import annotations

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import concentration, engine, nn, params, salt  # noqa: E402

SEQS = [
    "GTCACTGGCATTCAAGAGCTACAT",
    "ACCTGTGACATTGCCTGGAGATTA",
    "CAAGTGCCTAAGGTCATCGGATTC",
    "ATGCGC",
    "GGGGCCCC",          # kendine tamamlayici
    "AAAAAAAAAAAAAAAAAAAA",
]


class TestAnalyticIdentities(unittest.TestCase):
    def test_theta_at_tm_is_half(self):
        """Her konsantrasyon modelinde theta(Tm) = 0.5 (1e-9)."""
        cond = salt.Conditions()
        for seq in SEQS:
            for model in concentration.MODELS:
                if model == "self_complementary" and not nn.is_self_complementary(seq):
                    continue
                if model == "equal":
                    m = concentration.build(seq, cond=cond, conc_model=model,
                                            total_strand_M=4 * cond.primer_M)
                else:
                    m = concentration.build(seq, cond=cond, conc_model=model)
                with self.subTest(seq=seq, model=model):
                    self.assertAlmostEqual(m.theta(m.tm_c), 0.5, delta=1e-9)

    def test_theta_monotone_decreasing(self):
        m = concentration.build(SEQS[0])
        temps = [40, 50, 55, 60, 65, 70, 80]
        vals = [m.theta(t) for t in temps]
        self.assertTrue(all(a > b for a, b in zip(vals, vals[1:])), vals)
        self.assertGreater(vals[0], 0.99)
        self.assertLess(vals[-1], 0.01)

    def test_tm_increases_with_concentration(self):
        """dTm/dlogC > 0 (primer fazlasi modeli)."""
        prev = None
        for p in (1e-8, 1e-7, 1e-6, 1e-5):
            m = concentration.build(SEQS[0], cond=salt.Conditions(primer_M=p))
            if prev is not None:
                self.assertGreater(m.tm_c, prev)
            prev = m.tm_c

    def test_tm_increases_with_sodium(self):
        """dTm/d[Na+] > 0."""
        prev = None
        for na in (0.01, 0.05, 0.1, 0.5):
            m = concentration.build(
                SEQS[0],
                cond=salt.Conditions(monovalent_M=na, magnesium_M=0.0, dntp_M=0.0),
                salt_model="owczarzy_2004",
            )
            if prev is not None:
                self.assertGreater(m.tm_c, prev)
            prev = m.tm_c

    def test_salt_correction_lowers_tm_below_1M(self):
        m = concentration.build(SEQS[0], cond=salt.Conditions())
        self.assertLess(m.tm_c, m.tm_1m_K - params.KELVIN)


class TestNearestNeighbour(unittest.TestCase):
    def test_reverse_complement_invariance(self):
        """Duplex ayni oldugu icin dH/dS ters tamamlayicida degismemeli."""
        for seq in SEQS:
            a = nn.watson_crick(seq)
            b = nn.watson_crick(nn.reverse_complement(seq))
            with self.subTest(seq=seq):
                self.assertAlmostEqual(a.dH_kcal, b.dH_kcal, places=10)
                self.assertAlmostEqual(a.dS_cal, b.dS_cal, places=10)

    def test_symmetry_correction_applied(self):
        """Kendine tamamlayici dizide dS'e -1.4 cal/(mol K) eklenir."""
        selfc = "GGGGCCCC"
        self.assertTrue(nn.is_self_complementary(selfc))
        got = nn.watson_crick(selfc)
        table = params.WC_SANTALUCIA_1998.data
        manual_dS = sum(table["nn"][s if s in table["nn"]
                                    else nn.reverse_complement(s)][1]
                        for s in (selfc[i:i + 2] for i in range(len(selfc) - 1)))
        manual_dS += table["init"]["GC"][1] * 2
        self.assertAlmostEqual(got.dS_cal, manual_dS + table["symmetry"][1], places=10)

    def test_manual_additivity(self):
        """ATGC: AT + TG + GC yiginlari + iki uc baslatmasi."""
        seq = "ATGC"
        t = params.WC_SANTALUCIA_1998.data
        steps = ["AT", "TG", "GC"]
        dH = dS = 0.0
        for s in steps:
            key = s if s in t["nn"] else nn.reverse_complement(s)
            dH += t["nn"][key][0]
            dS += t["nn"][key][1]
        dH += t["init"]["AT"][0] + t["init"]["GC"][0]
        dS += t["init"]["AT"][1] + t["init"]["GC"][1]
        got = nn.watson_crick(seq)
        self.assertAlmostEqual(got.dH_kcal, dH, places=10)
        self.assertAlmostEqual(got.dS_cal, dS, places=10)

    def test_dG_sign_and_magnitude(self):
        """Uzun GC-zengin duplex belirgin negatif dG37 vermeli."""
        d = nn.watson_crick("GCGCGCGCGCGCGCGCGCGC")
        self.assertLess(d.dG_kcal(37.0), -20.0)

    def test_alphabet_contract(self):
        for bad in ("ACGU", "ACGTN", "acgt-", "ACG T"):
            with self.subTest(seq=bad):
                with self.assertRaises(ValueError):
                    nn.validate(bad)
        self.assertEqual(nn.validate("acgt"), "ACGT")

    def test_too_short(self):
        with self.assertRaises(ValueError):
            nn.watson_crick("A")


class TestSalt(unittest.TestCase):
    def test_spread_is_positive_and_finite(self):
        m = concentration.build(SEQS[0])
        sp = salt.spread(m.tm_1m_K, m.f_gc, m.duplex.n_bp, m.cond)
        self.assertGreater(sp["spread_c"], 0.0)
        self.assertTrue(all(math.isfinite(v) for v in sp.values()))

    def test_free_magnesium_floor(self):
        c = salt.Conditions(magnesium_M=0.0005, dntp_M=0.0008)
        self.assertEqual(c.free_magnesium_M, 0.0)

    def test_regime_selection(self):
        """R = sqrt([Mg])/[Mon] kurali: dusuk R'de tek degerlikli rejime duser."""
        _, regime = salt.owczarzy_2008_divalent(
            340.0, mg_M=1e-8, mon_M=0.5, f_gc=0.5, n_bp=20)
        self.assertEqual(regime, "monovalent_dominant")
        _, regime = salt.owczarzy_2008_divalent(
            340.0, mg_M=0.003, mon_M=0.05, f_gc=0.5, n_bp=20)
        self.assertEqual(regime, "mixed_mono_divalent")
        _, regime = salt.owczarzy_2008_divalent(
            340.0, mg_M=0.003, mon_M=0.001, f_gc=0.5, n_bp=20)
        self.assertEqual(regime, "divalent_dominant")

    def test_von_ahsen_equivalent(self):
        """Katsayi mM birimindedir; molar girdi/cikti icin cevrilmeli."""
        na = salt.von_ahsen_2001_na_eq(0.05, 0.003, 0.0006)
        self.assertAlmostEqual(na, 0.05 + 120 * math.sqrt(2.4) / 1e3, places=12)
        # fiziksel akil saglamasi: 50 mM Na + 2.4 mM serbest Mg -> 0.1-0.5 M arasi
        self.assertTrue(0.1 < na < 0.5, na)

    def test_salt_models_agree_within_few_degrees(self):
        """TH-02: uc tuz modeli tipik PCR kosullarinda birkac derece icinde olmali."""
        m = concentration.build("GTCACTGGCATTCAAGAGCTACAT")
        sp = salt.spread(m.tm_1m_K, m.f_gc, m.duplex.n_bp, m.cond)
        self.assertLess(sp["spread_c"], 5.0, sp)


class TestParameterRegistry(unittest.TestCase):
    def test_wc_table_shape(self):
        d = params.WC_SANTALUCIA_1998.data
        self.assertEqual(len(d["nn"]), 10)
        self.assertEqual(set(d["init"]), {"AT", "GC"})
        self.assertEqual(d["symmetry"], [0.0, -1.4])
        for v in d["nn"].values():
            self.assertLess(v[0], 0.0)   # tum WC yiginlari ekzotermik
            self.assertLess(v[1], 0.0)

    def test_checksum_deterministic(self):
        a = params.WC_SANTALUCIA_1998.sha256
        b = params.WC_SANTALUCIA_1998.sha256
        self.assertEqual(a, b)
        self.assertEqual(len(a), 64)
        self.assertEqual(len(params.registry_digest()), 64)

    def test_empty_table_raises(self):
        """Aktarilmamis tablo sessizce varsayilan uretmez, hata firlatir."""
        table = params.ParameterTable(params.MISMATCH_ALLAWI_PEYRET.provenance)
        self.assertFalse(table.is_available)
        with self.assertRaises(params.ParameterTableUnavailable):
            table.require()

    def test_dangling_table_is_available_with_unverified_provenance(self):
        table = params.DANGLING_BOMMARITO_2000
        self.assertEqual(len(table.require()), 32)
        row = next(row for row in params.registry_summary()
                   if row["table_id"] == table.provenance.table_id)
        self.assertEqual(row["entries"], "32")
        self.assertEqual(row["verification_status"], "transcribed_unverified")

    def test_mismatch_duplex_refuses_terminal_and_multiple_mismatches(self):
        with self.assertRaises(ValueError):
            nn.duplex_with_mismatches("ACGTACGTAC", "ACGTACGTAG")
        with self.assertRaises(ValueError):
            nn.duplex_with_mismatches("ACGTACGTAC", "AAGTACGTAC")
        with self.assertRaises(ValueError):
            nn.duplex_with_mismatches("ACGTACGTAC", "ACATTCGTAC")

    def test_single_internal_mismatch_matches_published_trimer(self):
        # AGC/TTG katkisi Allawi ve SantaLucia (1997), Tablo 4.
        primer, site = "CAGCG", "CAACG"
        duplex = nn.duplex_with_mismatches(primer, site)
        self.assertAlmostEqual(duplex.dH_kcal, -22.3)
        self.assertAlmostEqual(duplex.dS_cal, -66.9)
        rotated = nn.duplex_with_mismatches(nn.reverse_complement(site),
                                          nn.reverse_complement(primer))
        self.assertAlmostEqual(duplex.dG_kcal(), rotated.dG_kcal())

    def test_perfect_match_passes_through(self):
        a = nn.duplex_with_mismatches("ACGTACGTAC", "ACGTACGTAC")
        b = nn.watson_crick("ACGTACGTAC")
        self.assertEqual(a.dH_kcal, b.dH_kcal)

    def test_mismatch_positions(self):
        self.assertEqual(nn.mismatch_positions("ACGTAC", "ACCTAG"), [2, 5])


class TestEngineRecords(unittest.TestCase):
    def setUp(self):
        self.eng = engine.ThermoEngine()

    def test_record_schema(self):
        r = self.eng.evaluate("P1", SEQS[0])
        self.assertEqual({c["id"] for c in r["checks"]},
                         {"TH-01", "TH-02", "TH-03", "TH-04", "TH-05",
                          "DIM-01", "DIM-02", "SS-01", "SS-03"})
        for c in r["checks"]:
            for field in ("id", "target", "value", "threshold", "severity",
                          "model", "conditions", "param_set_sha256",
                          "engine_version", "profile"):
                self.assertIn(field, c)
            self.assertIn(c["severity"], ("INFO", "WARN", "FAIL"))
            self.assertEqual(len(c["param_set_sha256"]), 64)

    def test_internal_consistency_check_passes(self):
        for seq in SEQS:
            r = self.eng.evaluate("x", seq)
            th05 = [c for c in r["checks"] if c["id"] == "TH-05"][0]
            with self.subTest(seq=seq):
                self.assertEqual(th05["severity"], "INFO")
                self.assertLess(th05["value"]["abs_error"], 1e-9)

    def test_not_evaluated_is_explicit(self):
        r = self.eng.evaluate("P1", SEQS[0])
        self.assertIn("SPC-04", r["not_evaluated"])
        self.assertFalse(r["calibration_status"]["parameter_tables_verified"])

    def test_fold_checks_report_their_own_reference_conditions(self):
        result = self.eng.evaluate("hairpin", "GGGAAACCC")
        checks = {row["id"]: row for row in result["checks"]}
        self.assertEqual(checks["SS-01"]["value"]["dot_bracket"], "(((...)))")
        for identifier in ("SS-01", "SS-03"):
            self.assertEqual(checks[identifier]["severity"], "INFO")
            self.assertEqual(checks[identifier]["conditions"],
                             {"temperature_c": 37.0, "monovalent_mM": 1000.0})
        unpaired = checks["SS-03"]["value"]["p_unpaired_by_base"]
        self.assertEqual(unpaired[4], 1.0)
        self.assertLess(unpaired[0], unpaired[4])

    def test_determinism(self):
        """REP-01: ayni girdi -> ayni cikti."""
        import json
        a = json.dumps(self.eng.evaluate("P1", SEQS[0]), sort_keys=True)
        b = json.dumps(engine.ThermoEngine().evaluate("P1", SEQS[0]), sort_keys=True)
        self.assertEqual(a, b)

    def test_extreme_sequence_flagged(self):
        """Poli-A 20-mer qPCR penceresinin cok altinda kalmali -> FAIL."""
        r = self.eng.evaluate("polyA", "AAAAAAAAAAAAAAAAAAAA")
        th01 = [c for c in r["checks"] if c["id"] == "TH-01"][0]
        self.assertEqual(th01["severity"], "FAIL")


class TestOracleDifferential(unittest.TestCase):
    """Gelistirme oracle'i: bagimsiz bir uygulamaya karsi WC tablosu kontrolu.

    Amac transkripsiyon hatasi yakalamaktir. Oracle urun bagimliligi degildir;
    kurulu degilse test atlanir.
    """

    def test_against_biopython(self):
        try:
            from Bio.SeqUtils import MeltingTemp as mt  # type: ignore
        except Exception:
            self.skipTest("Biopython yok; oracle testi atlandi")

        conc_nM = 250.0
        cond = salt.Conditions(primer_M=conc_nM * 1e-9)
        for seq in SEQS:
            if nn.is_self_complementary(seq):
                continue
            m = concentration.build(
                seq, cond=cond, conc_model="equal",
                total_strand_M=4 * cond.primer_M,
            )
            ours = m.tm_1m_K - params.KELVIN
            theirs = mt.Tm_NN(
                seq, nn_table=mt.DNA_NN3, dnac1=2 * conc_nM, dnac2=2 * conc_nM,
                selfcomp=False, saltcorr=0, Na=50, Mg=0, dNTPs=0,
            )
            with self.subTest(seq=seq):
                self.assertAlmostEqual(ours, theirs, delta=0.05)


class TestModelConstants(unittest.TestCase):
    """Beyan edilmis model sabitleri kayitli ve denetlenebilir olmali.

    Bu sabitler literatur tablosu DEGILDIR; bu yuzden ayri bir kayitta
    dururlar. Ama uretilen sayilari degistirdikleri icin `registry_digest`e
    girmek ve her rapora yazilmak zorundadirlar (REP-02).
    """

    def test_at_bulge_closure_is_registered(self):
        const = params.MODEL_CONSTANTS["at_bulge_closure_dg37"]
        self.assertEqual(const, params.AT_BULGE_CLOSURE)
        self.assertEqual(const.status, "project_modelling_constant")
        # Literatur tablosu olmadigini acikca soylemeli.
        self.assertIn("degil", const.rationale.lower())

    def test_align_uses_the_registered_value_not_a_literal(self):
        from oligoscout_thermo import align
        self.assertEqual(align._AT_BULGE_CLOSURE_DG37,
                         params.AT_BULGE_CLOSURE.value)

    def test_constant_enters_the_registry_digest(self):
        """Degeri degistirmek `param_set_sha256`yi degistirmeli."""
        before = params.registry_digest()
        original = params.AT_BULGE_CLOSURE
        patched = params.ModelConstant(
            const_id=original.const_id, value=original.value + 0.1,
            units=original.units, applies_to=original.applies_to,
            rationale=original.rationale, status=original.status,
        )
        params.MODEL_CONSTANTS[original.const_id] = patched
        try:
            self.assertNotEqual(params.registry_digest(), before)
        finally:
            params.MODEL_CONSTANTS[original.const_id] = original
        self.assertEqual(params.registry_digest(), before)

    def test_summary_rows_are_complete(self):
        rows = params.model_constants_summary()
        self.assertTrue(rows)
        for row in rows:
            for field in ("const_id", "value", "units", "applies_to",
                          "rationale", "status"):
                self.assertTrue(row[field], f"{field} bos")

    def test_dimer_model_metadata_declares_the_constant(self):
        from oligoscout_thermo import align
        meta = align.model_metadata(allow_bulge=True)
        ids = [c["const_id"] for c in meta["model_constants"]]
        self.assertIn("at_bulge_closure_dg37", ids)
        self.assertNotIn("model_constants", align.model_metadata(False))


if __name__ == "__main__":
    unittest.main(verbosity=2)
