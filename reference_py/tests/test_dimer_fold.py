"""Dimer hizalama ve katlanma modullerinin dogrulamasi.

Kabul olcutu kaba kuvvettir: kucuk dizilerde tum olasiliklar tek tek sayilip
DP sonucuyla karsilastirilir (QC spesifikasyonu, bolum 5).
"""

from __future__ import annotations

import math
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import align, fold, nn, params  # noqa: E402

E = fold.DeclaredLoopEnergies()

SHORT = ["GGGAAACCC", "GCGAAAGCTT", "ACGTAAAACGT", "GGATTTTATCC",
         "GCGCAAAAGCGC", "ATATAAAATATA", "GGGGAAAACCCC", "ACGTACGTACGT"]

DIMER_PAIRS = [("ACGTACGT", "ACGTACGT"), ("GGCCTTAA", "TTAAGGCC"),
               ("ATGCATGC", "GCATGCAT"), ("AAAAAAAA", "TTTTTTTT"),
               ("TTCCAGCGGTT", "AACCGTTGGAA"), ("CCAGCGG", "CCGTTGG"),
               ("GGCAAGTGCC", "GGCATTCGCC")]

_PAIRS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}


def brute_contiguous(a, b, mode):
    """Kesintisiz helis / tek ic uyumsuzluk sinifinin bagimsiz sayimi."""
    rb = b[::-1]
    energies = []
    sym_h, sym_s = params.WC_SANTALUCIA_1998.require()["symmetry"]
    for offset in range(1 - len(rb), len(a)):
        lo, hi = max(0, offset), min(len(a), offset + len(rb))
        for start in range(lo, hi - 1):
            for end in range(start + 2, hi + 1):
                if mode == "3prime" and start != offset and end != len(a):
                    continue
                top = a[start:end]
                site = nn.reverse_complement(rb[start - offset:end - offset][::-1])
                try:
                    duplex = nn.duplex_with_mismatches(top, site)
                except ValueError:
                    continue
                energy = duplex.dG_kcal()
                if duplex.self_complementary:
                    energy -= sym_h - 310.15 * sym_s / 1000.0
                energies.append(energy)
    return min(energies, default=None)


def brute_single_bulge(a, b, mode):
    """Tek bazlik ic bulge sinifinin bagimsiz sayimi.

    `align._scan_bulge`ten farkli bir parametrelendirme kullanir: burada
    yapilar dogrudan (sol helis baslangici, uzunluk, bulge ipligi, sag
    helis uzunlugu) olarak sayilir, kayma aritmetigi yoktur.
    """
    rb = b[::-1]
    table = params.WC_SANTALUCIA_1998.require()
    stacks, init = table["nn"], table["init"]
    bulge_dg = params.DNA_LOOP_SANTALUCIA_HICKS_2004.require()["BULGE:1"][0]

    def step(pair):
        return stacks.get(pair) or stacks[nn.reverse_complement(pair)]

    def helix(seq):
        dh = ds = 0.0
        for i in range(len(seq) - 1):
            h, s = step(seq[i:i + 2])
            dh += h
            ds += s
        return dh, ds

    def matches(i, j, length):
        return all((a[i + t], rb[j + t]) in _PAIRS for t in range(length))

    na, nb = len(a), len(rb)
    best = None
    for i in range(na):
        for j in range(nb):
            for l1 in range(2, min(na - i, nb - j) + 1):
                if not matches(i, j, l1):
                    break
                for strand in ("a", "b"):
                    i2 = i + l1 + 1 if strand == "a" else i + l1
                    j2 = j + l1 if strand == "a" else j + l1 + 1
                    if i2 >= na or j2 >= nb:
                        continue
                    for l2 in range(2, min(na - i2, nb - j2) + 1):
                        if not matches(i2, j2, l2):
                            break
                        terminal = (i2 + l2 - 1 == na - 1) or (j == 0)
                        if mode == "3prime" and not terminal:
                            continue
                        dh, ds = 0.0, 0.0
                        for seq in (a[i:i + l1], a[i2:i2 + l2]):
                            h, s = helix(seq)
                            dh += h
                            ds += s
                        for base in (a[i], a[i2 + l2 - 1]):
                            h, s = init["AT" if base in "AT" else "GC"]
                            dh += h
                            ds += s
                        closing = (params.AT_BULGE_CLOSURE.value
                                   * sum(base in "AT"
                                         for base in (a[i + l1 - 1], a[i2])))
                        ds -= 1000.0 * (bulge_dg + closing) / 310.15
                        if strand == "b":
                            h, s = step(a[i + l1 - 1] + a[i + l1])
                        else:
                            h, s = step(rb[j + l1] + rb[j + l1 - 1])
                        dh += h
                        ds += s
                        energy = dh - 310.15 * ds / 1000.0
                        if best is None or energy < best:
                            best = energy
    return best


class TestDimerAlignment(unittest.TestCase):
    def test_symmetry(self):
        """DIM-08: dimer(a,b) ile dimer(b,a) ayni enerjiyi vermeli."""
        pairs = [("GTCACTGGCATTCAAGAGCTACAT", "ACCTGTGACATTGCCTGGAGATTA"),
                 ("GCGCGCGC", "GCGCGCGC"), ("AAAATTTT", "AAAATTTT"),
                 ("CCAGCGG", "CCGTTGG")]
        for a, b in pairs:
            for mode in align.MODES:
                with self.subTest(a=a, b=b, mode=mode):
                    self.assertAlmostEqual(align.best_dimer(a, b, mode).dG37_kcal,
                                           align.best_dimer(b, a, mode).dG37_kcal, places=9)

    def test_brute_force_offsets(self):
        """Tum gecerli alt helisler, iki 3' uc kurali icin bagimsiz sayilir.

        Bulge kapali: bu test yalniz kesintisiz helis / tek ic uyumsuzluk
        sinifini kapsar. Bulge dahil tarama TestSingleBaseBulge'dedir.
        """
        for a, b in DIMER_PAIRS:
            for mode in align.MODES:
                with self.subTest(a=a, b=b, mode=mode):
                    actual = align.best_dimer(a, b, mode, allow_bulge=False)
                    expected = brute_contiguous(a, b, mode)
                    self.assertEqual(actual.found, expected is not None)
                    if actual.found:
                        self.assertAlmostEqual(actual.dG37_kcal, expected, places=9)

    def test_single_mismatch_is_marked_in_eligible_3prime_dimer(self):
        result = align.best_dimer("CCAGCGG", "CCGTTGG", "3prime")
        self.assertEqual(result.match, "|||X|||")
        self.assertEqual(result.mismatch_count, 1)
        self.assertEqual(result.helix_len, 6)
        self.assertTrue(result.three_prime_paired)
        self.assertAlmostEqual(result.dH_kcal, -38.3)
        self.assertAlmostEqual(result.dS_cal, -106.7)

    def test_perfect_complement_is_most_stable(self):
        a = "GCGCATATGCGC"
        res = align.best_dimer(a, align.reverse_complement(a))
        self.assertEqual(res.helix_len, len(a))
        self.assertLess(res.dG37_kcal, -10.0)

    def test_three_prime_mode_requires_terminal_pair(self):
        a, b = "AAAAAAGGGGCC", "GGCCCCAAAAAA"
        any_mode = align.best_dimer(a, b, "any")
        tp_mode = align.best_dimer(a, b, "3prime")
        self.assertTrue(any_mode.found)
        if tp_mode.found:
            self.assertTrue(tp_mode.three_prime_paired)
            self.assertGreaterEqual(tp_mode.dG37_kcal, any_mode.dG37_kcal - 1e-9)

    def test_no_dimer_reports_not_found(self):
        self.assertFalse(align.best_dimer("AAAAAAAA", "AAAAAAAA").found)

    def test_tm_increases_with_stability(self):
        weak = align.best_dimer("ACGTAAAA", "TTTTACGT")
        strong = align.best_dimer("GCGCGCGC", "GCGCGCGC")
        self.assertGreater(strong.tm_c(2e-7), weak.tm_c(2e-7))


class TestSingleBaseBulge(unittest.TestCase):
    a, b = "GCGCTGCGC", "GCGCGCGC"

    def test_production_marks_gap_in_both_orders(self):
        for a, b, strand in ((self.a, self.b, "a"), (self.b, self.a, "b")):
            result = align.best_dimer(a, b)
            self.assertEqual((result.bulge_count, result.bulge_strand), (1, strand))
            self.assertEqual(result.match, "||||B||||")
            self.assertEqual(result.as_dict()["bulge_count"], 1)
            self.assertEqual(result.helix_len, 8)

    def test_bulge_switch_preserves_baseline(self):
        baseline = align.best_dimer(self.a, self.b, allow_bulge=False)
        self.assertEqual(baseline.bulge_count, 0)
        self.assertLess(align.best_dimer(self.a, self.b).dG37_kcal, baseline.dG37_kcal)
        self.assertEqual(align.model_metadata(False)["model"], align.MODEL_NAME_NO_BULGE)

    def test_bulge_symmetry_in_each_mode(self):
        for mode in align.MODES:
            self.assertAlmostEqual(align.best_dimer(self.a, self.b, mode).dG37_kcal,
                                   align.best_dimer(self.b, self.a, mode).dG37_kcal, places=9)

    def test_single_wc_pair_on_one_side_is_insufficient(self):
        self.assertFalse(align._scan_bulge("CAT", "GA", "any").found)

    def test_three_prime_anchor_excludes_internal_bulge(self):
        a, rb = "A" + self.a + "A", ("A" + self.b + "A")[::-1]
        self.assertTrue(align._scan_bulge(a, rb, "any").found)
        self.assertFalse(align._scan_bulge(a, rb, "3prime").found)

    def test_against_exhaustive_short_dimer_enumeration(self):
        rng = random.Random(53)
        pairs = [(self.a, self.b), ("ATCGACCG", "CGGTGAT")]
        pairs.extend(("".join(rng.choices("ACGT", k=9)), "".join(rng.choices("ACGT", k=8)))
                     for _ in range(70))
        for a, b in pairs:
            for mode in align.MODES:
                energies = [v for v in (brute_contiguous(a, b, mode), brute_single_bulge(a, b, mode))
                            if v is not None]
                result = align.best_dimer(a, b, mode)
                self.assertEqual(result.found, bool(energies), (a, b, mode))
                if energies:
                    self.assertAlmostEqual(result.dG37_kcal, min(energies), places=9)


class TestFoldAgainstBruteForce(unittest.TestCase):
    def test_mfe_matches_enumeration(self):
        for seq in SHORT:
            with self.subTest(seq=seq):
                dp, structure = fold.zuker_mfe(seq, E)
                brute = min(fold.structure_energy(seq, st, E)
                            for st in fold.enumerate_structures(seq))
                self.assertAlmostEqual(dp, min(brute, 0.0), places=9)
                self.assertEqual(len(structure), len(seq))
                self.assertEqual(structure.count("("), structure.count(")"))

    def test_traceback_structure_has_mfe_energy(self):
        for seq in SHORT:
            dp, structure = fold.zuker_mfe(seq, E)
            pairs, stack = [], []
            for idx, ch in enumerate(structure):
                if ch == "(":
                    stack.append(idx)
                elif ch == ")":
                    pairs.append((stack.pop(), idx))
            with self.subTest(seq=seq):
                self.assertAlmostEqual(fold.structure_energy(seq, pairs, E),
                                       dp, places=9)

    def test_partition_function_matches_enumeration(self):
        rt = params.R_CAL * (37.0 + params.KELVIN) / 1000.0
        for model in (E, fold.DnaLoopEnergies()):
            for seq in SHORT + ["GGAAACGAAACC"]:
                with self.subTest(seq=seq, model=model.name):
                    z_dp, _ = fold.mccaskill(seq, model)
                    z_brute = sum(math.exp(-fold.structure_energy(seq, st, model) / rt)
                                  for st in fold.enumerate_structures(seq))
                    self.assertAlmostEqual(z_dp / z_brute, 1.0, places=9)

    def test_base_pair_probabilities_match_enumeration(self):
        rt = params.R_CAL * (37.0 + params.KELVIN) / 1000.0
        for seq in SHORT[:5]:
            z_dp, probs = fold.mccaskill(seq, E)
            n = len(seq)
            acc = [[0.0] * n for _ in range(n)]
            for st in fold.enumerate_structures(seq):
                w = math.exp(-fold.structure_energy(seq, st, E) / rt)
                for (i, j) in st:
                    acc[i][j] += w
            with self.subTest(seq=seq):
                for i in range(n):
                    for j in range(i + 1, n):
                        self.assertAlmostEqual(probs[i][j], acc[i][j] / z_dp,
                                               places=9)

    def test_probabilities_bounded(self):
        _, probs = fold.mccaskill("GGGAAACCCGGGAAACCC", E)
        for row in probs:
            self.assertLessEqual(sum(row), 1.0 + 1e-9)
            self.assertTrue(all(0.0 <= p <= 1.0 for p in row))

    def test_unpaired_probability_complement(self):
        seq = "GCGAAAGCTT"
        _, probs = fold.mccaskill(seq, E)
        unp = fold.unpaired_probability(seq, E)
        for i, u in enumerate(unp):
            self.assertAlmostEqual(u, 1.0 - sum(probs[i]), places=9)

    def test_ensemble_free_energy_below_mfe(self):
        """Toplulugun serbest enerjisi MFE'den dusuk ya da esit olmali."""
        for seq in SHORT:
            mfe, _ = fold.zuker_mfe(seq, E)
            with self.subTest(seq=seq):
                self.assertLessEqual(fold.ensemble_free_energy(seq, E), mfe + 1e-9)


class TestFoldModels(unittest.TestCase):
    def test_fold_without_loop_parameters_raises(self):
        for fn in (fold.zuker_mfe, fold.mccaskill, fold.unpaired_probability):
            with self.subTest(fn=fn.__name__):
                with self.assertRaises(params.ParameterTableUnavailable):
                    fn("GGGAAACCCGGG", fold.PendingLoopEnergies())

    def test_default_dna_model_uses_published_loop_lengths(self):
        model = fold.DnaLoopEnergies()
        self.assertEqual(model.hairpin(3), 3.5)
        self.assertEqual(model.bulge(2), 2.9)
        self.assertAlmostEqual(model.internal(1, 2), 3.5)
        self.assertTrue(math.isinf(model.internal(1, 1)))
        self.assertTrue(math.isinf(model.ml_a()))
        self.assertLess(model.hairpin(10), model.hairpin(11))
        self.assertLess(model.hairpin(11), model.hairpin(12))
        self.assertGreater(model.hairpin(31), model.hairpin(30))
        seq = "GGGAAACCC"
        energy, structure = fold.zuker_mfe(seq)
        # Iki GG/CC yigininin dG37 degeri ve 3 bazlik sac tokasi cezasi.
        self.assertAlmostEqual(energy, 3.5 + 2 * (-8.0 + 310.15 * 19.9 / 1000))
        self.assertEqual(structure, "(((...)))")
        self.assertEqual(fold.structure_energy("AAAAT", [(0, 4)], model), 4.0)
        self.assertAlmostEqual(fold.structure_energy("GAGAAACC", [(0, 7), (2, 6)], model),
                               7.5 - 8.0 + 310.15 * 19.9 / 1000)
        with self.assertRaises(ValueError):
            fold.mccaskill(seq, temp_c=60.0)

    def test_declared_model_is_labelled(self):
        self.assertIn("not_literature", E.name)

    def test_dimer_limitation_is_declared(self):
        self.assertIn("kesintisiz", align.__doc__)


if __name__ == "__main__":
    unittest.main(verbosity=2)
