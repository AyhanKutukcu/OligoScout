"""QC katalogu YAML'i gecerli ve uretilen TSV'lerle tutarli olmali."""

from __future__ import annotations

import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

try:
    import yaml  # noqa: F401
    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False


@unittest.skipUnless(HAVE_YAML, "PyYAML yok (yalniz katalog araci icin gerekir)")
class CatalogYaml(unittest.TestCase):
    def test_valid_and_in_sync(self):
        if not os.path.exists(os.path.join(ROOT, "out", "durum", "qc_katalog_kapsami.tsv")):
            self.skipTest("uretilen katalog TSV yok")
        done = subprocess.run([sys.executable, "cli/catalog_yaml.py", "check"], cwd=ROOT,
                              capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("91 kural", done.stdout)

    def test_every_rule_declares_calibration(self):
        doc = yaml.safe_load(open(os.path.join(ROOT, "catalog", "qc_catalog.yaml"), encoding="utf-8"))
        for rule in doc["rules"]:
            self.assertIn(rule["calibration"]["status"], ("uncalibrated", "calibrated"), rule["id"])


if __name__ == "__main__":
    unittest.main()
