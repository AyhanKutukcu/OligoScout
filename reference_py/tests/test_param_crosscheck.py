"""Parametre tablolarinin bagimsiz kaynaklara karsi diferansiyel denetimi.

Kaynak (primer3) bulunamazsa testler ATLANIR: referans katmani cevrimdisi
ve bagimsiz calisabilir kalmalidir. Kaynak varsa dort tablonun tamami
tolerans icinde olmali ve kunyedeki `cross_check` alani canli sayimla
birebir tutmalidir (bayat kunye yakalanir).
"""

from __future__ import annotations

import os
import re
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import crosscheck as xc, params  # noqa: E402

SOURCES = xc.oracle_sources()
ENTRIES = xc.run_crosscheck() if SOURCES else []
SUMMARY = xc.summarise(ENTRIES) if ENTRIES else []


@unittest.skipUnless(SOURCES, "primer3 kaynagi yok; diferansiyel denetim atlandi")
class TestOracleParsing(unittest.TestCase):
    """Kaynak dosyalarin sozlesmesi beklendigi gibi mi."""

    def test_thal_tables_have_expected_sizes(self):
        if "config_dir" not in SOURCES:
            self.skipTest("primer3_config yok")
        thal = xc.load_thal_oracle(SOURCES["config_dir"])
        # 16 din. kombinasyon, 96 fiziksel tekli uyumsuzluk yigini,
        # 2 uc x 4 kapanis cifti x 4 sarkan baz = 32 sarkan uc.
        self.assertEqual(len(thal["wc"]), 16)
        self.assertEqual(len(thal["mismatch"]), 96)
        self.assertEqual(len(thal["dangling"]), 32)
        self.assertGreaterEqual(len(thal["loops"]), 47)

    def test_known_values_are_recovered(self):
        """Elle dogrulanmis birkac capa deger; ayristirma yonu dogru mu."""
        if "config_dir" not in SOURCES:
            self.skipTest("primer3_config yok")
        thal = xc.load_thal_oracle(SOURCES["config_dir"])
        self.assertAlmostEqual(thal["wc"]["AA"][0], -7.9, places=4)
        self.assertAlmostEqual(thal["wc"]["CG"][1], -27.2, places=4)
        self.assertAlmostEqual(thal["mismatch"]["AA/AT"][0], 4.7, places=4)
        self.assertAlmostEqual(thal["dangling"]["D3:A:C-G"][0], -5.9, places=4)
        self.assertAlmostEqual(thal["dangling"]["D5:T:A-T"][0], -6.9, places=4)
        self.assertAlmostEqual(thal["loops"]["BULGE:1"], 4.0, places=2)

    def test_oligotm_initiation_terms(self):
        if "oligotm_c" not in SOURCES:
            self.skipTest("oligotm.c yok")
        tm = xc.load_oligotm_oracle(SOURCES["oligotm_c"])
        self.assertAlmostEqual(tm["init"]["AT"][0], 2.3, places=4)
        self.assertAlmostEqual(tm["init"]["GC"][1], -2.8, places=4)
        self.assertAlmostEqual(tm["symmetry"][1], -1.4, places=4)

    def test_mismatch_rotation_is_an_involution(self):
        key = "AA/CT"
        self.assertEqual(xc._rotate_mismatch_key(xc._rotate_mismatch_key(key)),
                         key)


@unittest.skipUnless(SOURCES, "primer3 kaynagi yok; diferansiyel denetim atlandi")
class TestCrossCheck(unittest.TestCase):

    def test_no_entry_exceeds_tolerance(self):
        bad = [e for e in ENTRIES if e.status == "DIFF"]
        self.assertEqual(bad, [], f"{len(bad)} girdi toleransi asiyor")

    def test_no_entry_is_missing_on_either_side(self):
        missing = [e for e in ENTRIES
                   if e.status in ("MISSING_OURS", "MISSING_ORACLE")]
        self.assertEqual(missing, [],
                         f"{len(missing)} girdi iki kaynaktan birinde yok")

    def test_every_registry_table_is_covered(self):
        covered = {row["table_id"] for row in SUMMARY}
        # 'pending' tablolarin bagimsiz oracle'i yoktur; capraz denetime
        # girmezler. Ayrica DNA/RNA hibrit NN tablosunun (Sugimoto 1995)
        # primer3'te bagimsiz oracle'i YOKTUR (primer3 yalniz DNA); aktarilmis
        # olsa bile capraz denetime giremez.
        no_oracle = {params.HYBRID_SUGIMOTO_1995.provenance.table_id}
        expected = {tid for tid, t in params.REGISTRY.items()
                    if t.provenance.verification_status != "pending"
                    and tid not in no_oracle}
        if "terminal_dir" not in SOURCES:
            expected.remove(params.TERMINAL_MISMATCH_SANTALUCIA_PEYRET.provenance.table_id)
        self.assertEqual(covered, expected)

    def test_terminal_source_notes_are_explicit(self):
        if "terminal_dir" not in SOURCES:
            self.skipTest("tstack2 yok")
        notes = {(e.key, e.quantity) for e in ENTRIES if e.status == "SOURCE_NOTE"}
        self.assertEqual(notes, {("AC/TC", "dS"), ("GG/CG", "dH")})

    def test_provenance_cross_check_matches_live_counts(self):
        """Kunyedeki sayimlar canli denetimle tutmali; bayat kunye hatadir."""
        live: dict[str, list[str]] = {}
        for row in SUMMARY:
            live.setdefault(row["table_id"], []).append(
                f"{row['oracle']} {row['match']}/{row['compared']} "
                f"max|d|={float(row['max_abs_delta'] or 0.0):.4f}"
            )
        for table_id, parts in live.items():
            recorded = params.REGISTRY[table_id].provenance.cross_check
            self.assertTrue(recorded,
                            f"{table_id}: cross_check alani bos")
            self.assertEqual(recorded, "; ".join(parts),
                             f"{table_id}: kunye canli sayimla tutmuyor")

    def test_cross_check_does_not_promote_verification_status(self):
        """Diferansiyel denetim 'verified' vermez; tanim degismedi."""
        for table in params.REGISTRY.values():
            self.assertNotEqual(table.provenance.verification_status,
                                "verified")

    def test_summary_rows_are_machine_readable(self):
        for row in SUMMARY:
            self.assertTrue(re.fullmatch(r"\d+", row["compared"]))
            self.assertGreater(int(row["compared"]), 0)


if __name__ == "__main__":
    unittest.main()
