"""panel_qc: tek primer urunleri (F-F / R-R, oligoscout_panel >= 0.4).

Gercek ctDNA panelinden (DESIGN_007) alinmis kucuk bir fikstur: hedef F-R
urunu + chr9'da yalniz R primerinin urettigi 515 bp R-R urunu.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CLI = os.path.join(ROOT, "cli", "panel_qc.py")

F = "CCTTATCAAAGCGGCTCTTG"
R = "AGATTCCGTGTCTGGTGAG"

BINDINGS = [
    ["pair_id", "role", "primer", "chrom", "start", "end_exclusive", "strand",
     "mismatches", "mask_bits", "mask", "site_plus_strand"],
    ["DESIGN_007", "R", R, "chr9", "26191557", "26191576", "+", "3", "74",
     ".X.X..X............", "AAACTCAGTGTCTGGTGAG"],
    ["DESIGN_007", "R", R, "chr9", "26192053", "26192072", "-", "1", "64",
     "......X............", "CTCACCAGACACTGAATCT"],
    ["DESIGN_007", "F", F, "chr22", "20000753", "20000773", "+", "0", "0",
     "....................", F],
    ["DESIGN_007", "R", R, "chr22", "20000830", "20000849", "-", "0", "0",
     "...................", "CTCACCAGACACGGAATCT"],
]
HEAD = ["pair_id", "chrom", "amplicon_start", "amplicon_end_exclusive",
        "amplicon_length", "left_position", "right_position", "left_mismatches",
        "right_mismatches", "total_mismatches", "left_mask", "right_mask"]
RR = ["DESIGN_007", "chr9", "26191557", "26192072", "515", "26191557", "26192053",
      "3", "1", "4", ".X.X..X............", "......X............"]
FR = ["DESIGN_007", "chr22", "20000753", "20000849", "96", "20000753", "20000830",
      "0", "0", "0", "....................", "..................."]


def write(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)


def run(tmp, products, panel=True):
    prefix = os.path.join(tmp, "panel")
    write(prefix + ".bindings.tsv", BINDINGS)
    write(prefix + ".products.tsv", products)
    write(prefix + ".contract.tsv", [["field", "value"], ["anchor_length", "12"]])
    args = [sys.executable, CLI, "--prefix", prefix, "--out", tmp]
    if panel:
        pp = os.path.join(tmp, "design_panel.tsv")
        write(pp, [["pair_id", "primer1", "primer2", "left_start", "right_start",
                    "expected_amplicon_length"],
                   ["DESIGN_007", F, R, "20000753", "20000849", "96"]])
        args += ["--panel", pp]
    subprocess.run(args, check=True, capture_output=True, cwd=ROOT)
    with open(os.path.join(tmp, "panel_qc.tsv"), encoding="utf-8") as fh:
        qc = next(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(tmp, "offtarget_products.tsv"), encoding="utf-8") as fh:
        off = list(csv.DictReader(fh, delimiter="\t"))
    return qc, off


class SinglePrimerProducts(unittest.TestCase):
    def test_rr_product_is_offtarget(self):
        rows = [HEAD + ["left_role", "right_role", "orientation"],
                RR + ["R", "R", "R-R"], FR + ["F", "R", "F-R"]]
        with tempfile.TemporaryDirectory() as tmp:
            qc, off = run(tmp, rows)
        self.assertEqual(qc["products_total"], "2")
        self.assertEqual(qc["products_intended"], "1")
        self.assertEqual(qc["offtarget_total"], "1")
        self.assertEqual(qc["offtarget_single_primer"], "1")
        self.assertEqual(qc["spc03"], "WARN")
        self.assertEqual(off[0]["orientation"], "R-R")
        self.assertEqual((off[0]["left_role"], off[0]["right_role"]), ("R", "R"))

    def test_exact_single_primer_never_intended(self):
        # panel koordinati yokken "0 uyumsuzluk = hedef" kurali tek primer
        # urununu hedef saymamali
        rr_exact = RR[:7] + ["0", "0", "0", "...................", "..................."]
        rows = [HEAD + ["left_role", "right_role", "orientation"],
                rr_exact + ["R", "R", "R-R"], FR + ["F", "R", "F-R"]]
        with tempfile.TemporaryDirectory() as tmp:
            qc, _ = run(tmp, rows, panel=False)
        self.assertEqual(qc["products_intended"], "1")
        self.assertEqual(qc["offtarget_exact"], "1")
        self.assertEqual(qc["spc03"], "FAIL")

    def test_legacy_products_without_orientation(self):
        with tempfile.TemporaryDirectory() as tmp:
            qc, off = run(tmp, [HEAD, FR])
        self.assertEqual(qc["products_total"], "1")
        self.assertEqual(qc["offtarget_single_primer"], "0")
        self.assertEqual(off, [])


if __name__ == "__main__":
    unittest.main()
