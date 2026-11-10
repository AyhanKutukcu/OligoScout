"""Ikincil yapi: Zuker MFE ve McCaskill bolusum fonksiyonu.

Varsayilan DNA modeli 37 C ve 1 M NaCl'de Watson-Crick yiginlarini,
SantaLucia & Hicks Tablo 4'teki uzunluk terimlerini ve ic dongu asimetrisini
kullanir. Cok dalli ve 1x1 ic donguler kapsam disindadir. Uc bazli sac
tokasinda A-T kapanisi ve tek bazlik bulge icin WC yigini uygulanir.
Diger uc terimleri ve kisa sac tokasi bonuslari kapsam disindadir.
Algoritma acikca verilen alternatif enerji modellerini de kabul eder.

Enerji ayristirmasi (DP, kaba kuvvet ve bolusum fonksiyonu ayni tanimi
kullanir, belirsizliksiz ayristirma):
  hairpin(i,j)   dongu boyutuna bagli ceza
  stack(i,j)     WC NN dG37
  bulge/ic dongu tek cocuklu durumlar
  multiloop      a + b*(1+cocuk) + c*eslesmemis
  dis dongu      eslesmemis baz ve dal basina ucretsiz

Dogrulama (tests/test_dimer_fold.py): n <= 14 icin tum pseudoknot'suz yapilar
sayilir; MFE, Z ve baz-cifti olasiliklari birebir karsilastirilir.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Iterator

from .nn import reverse_complement, validate
from .params import (
    DNA_LOOP_SANTALUCIA_HICKS_2004, KELVIN, R_CAL, ParameterTableUnavailable,
    WC_SANTALUCIA_1998,
)

_PAIR = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
MIN_HAIRPIN = 3
INF = float("inf")


def can_pair(x: str, y: str) -> bool:
    return (x, y) in _PAIR


@lru_cache(maxsize=256)
def _dna_loop_length_dg37(kind: str, size: int) -> float:
    table = DNA_LOOP_SANTALUCIA_HICKS_2004.require()
    if f"{kind}:{size}" in table:
        return table[f"{kind}:{size}"][0]
    sizes = sorted(int(key.split(":")[1]) for key in table if key.startswith(f"{kind}:"))
    if size < sizes[0]:
        return INF
    lower = max(n for n in sizes if n < size)
    g_lower = table[f"{kind}:{lower}"][0]
    if size > sizes[-1]:
        return g_lower + 2.44 * R_CAL * (37.0 + KELVIN) / 1000.0 * math.log(size / lower)
    upper = min(n for n in sizes if n > size)
    g_upper = table[f"{kind}:{upper}"][0]
    fraction = math.log(size / lower) / math.log(upper / lower)
    return g_lower + fraction * (g_upper - g_lower)


@dataclass(frozen=True)
class DnaLoopEnergies:
    """DNA icin 37 C'de uzunluga dayali, cok dalsiz model."""

    name: str = "dna_loop_length_37c"

    def hairpin(self, size: int) -> float:
        return _dna_loop_length_dg37("HAIRPIN", size)

    def bulge(self, size: int) -> float:
        return _dna_loop_length_dg37("BULGE", size)

    def internal(self, s1: int, s2: int) -> float:
        return _dna_loop_length_dg37("INTERNAL", s1 + s2) + 0.3 * abs(s1 - s2)

    def ml_a(self) -> float: return INF  # cok dalli yapinin Boltzmann agirligi sifir
    def ml_b(self) -> float: return 0.0
    def ml_c(self) -> float: return 0.0


@dataclass(frozen=True)
class PendingLoopEnergies:
    """Aktarilmamis dongu parametreleri. Her cagri hata firlatir."""

    name: str = "pending_loop_parameters"
    citation: str = ("SantaLucia & Hicks 2004 (DNA dongu parametreleri); "
                     "Turner/Mathews 2004 (RNA)")

    def _fail(self) -> float:
        raise ParameterTableUnavailable(
            "Dongu (hairpin/bulge/ic dongu/multiloop) parametreleri "
            f"aktarilmadi. Kaynak: {self.citation}. Katlanma enerjisi "
            "hesaplanamaz; once parametreleri girin ya da acik bir model verin."
        )

    def hairpin(self, size: int) -> float: return self._fail()
    def bulge(self, size: int) -> float: return self._fail()
    def internal(self, s1: int, s2: int) -> float: return self._fail()
    def ml_a(self) -> float: return self._fail()
    def ml_b(self) -> float: return self._fail()
    def ml_c(self) -> float: return self._fail()


@dataclass(frozen=True)
class DeclaredLoopEnergies:
    """Acikca beyan edilmis, LITERATUR OLMAYAN dongu modeli.

    Yalniz algoritma dogrulamasi icindir; bilimsel sonuc uretmez. Adi
    (`name`) her rapora yazilir ki bir literatur parametresiyle karistirilmasin.
    """

    name: str = "declared_test_model_not_literature"
    a: float = 3.4
    b: float = 0.4
    c: float = 0.1

    def hairpin(self, size: int) -> float:
        return 4.0 + 1.2 * math.log(max(size, 3) / 3.0)

    def bulge(self, size: int) -> float:
        return 3.5 + 0.5 * math.log(max(size, 1))

    def internal(self, s1: int, s2: int) -> float:
        return 2.5 + 0.4 * (s1 + s2) + 0.3 * abs(s1 - s2)

    def ml_a(self) -> float: return self.a
    def ml_b(self) -> float: return self.b
    def ml_c(self) -> float: return self.c


def _wc_step_dg37(step: str) -> float:
    table = WC_SANTALUCIA_1998.require()["nn"]
    h, s = table.get(step) or table[reverse_complement(step)]
    return h - 310.15 * s / 1000.0


def stack_dG37(seq: str, i: int, j: int) -> float:
    return _wc_step_dg37(seq[i : i + 2])


def _hairpin_pen(e, seq, i, j):
    size = j - i - 1
    penalty = e.hairpin(size)
    if isinstance(e, DnaLoopEnergies) and size == 3 and seq[i] in "AT":
        penalty += 0.5  # denklem 8: A-T kapanis cezasi
    return penalty


def _inner_pen(e, seq, i, j, p, q) -> float:
    s1, s2 = p - i - 1, j - q - 1
    if s1 == 0 and s2 == 0:
        return stack_dG37(seq, i, j)
    if s1 == 0 or s2 == 0:
        penalty = e.bulge(s1 + s2)
        if isinstance(e, DnaLoopEnergies) and s1 + s2 == 1:
            # Denklem 13: tek bazlik bulge uzerinden komsu WC yigini.
            penalty += _wc_step_dg37(seq[i] + seq[p])
        return penalty
    return e.internal(s1, s2)


def _inner_pairs(n, i, j, max_loop):
    for p in range(i + 1, min(i + max_loop + 2, j)):
        for q in range(max(p + MIN_HAIRPIN + 1, j - max_loop - 1), j):
            if (p - i - 1) + (j - q - 1) <= max_loop:
                yield p, q


# ---------------------------------------------------------------------------
# MFE (Zuker)
# ---------------------------------------------------------------------------
def _fill_mfe(s, e, max_loop):
    n = len(s)
    V = [[INF] * n for _ in range(n)]
    M1 = [[INF] * n for _ in range(n)]
    M = [[INF] * n for _ in range(n)]
    W = [[0.0] * (n + 1) for _ in range(n + 1)]   # W[i][j], j+1 ofsetli

    def w(i, j):
        return 0.0 if i > j else W[i][j]

    for span in range(0, n):
        for i in range(0, n - span):
            j = i + span
            if span > MIN_HAIRPIN and can_pair(s[i], s[j]):
                best = _hairpin_pen(e, s, i, j)
                for p, q in _inner_pairs(n, i, j, max_loop):
                    if V[p][q] < INF:
                        best = min(best, _inner_pen(e, s, i, j, p, q) + V[p][q])
                for k in range(i + 2, j):
                    if M[i + 1][k - 1] < INF and M1[k][j - 1] < INF:
                        best = min(best, e.ml_a() + e.ml_b()
                                   + M[i + 1][k - 1] + M1[k][j - 1])
                V[i][j] = best

            m1 = INF
            for k in range(i + MIN_HAIRPIN + 1, j + 1):
                if V[i][k] < INF:
                    m1 = min(m1, V[i][k] + e.ml_b() + e.ml_c() * (j - k))
            M1[i][j] = m1

            m = INF
            for k in range(i, j + 1):
                if M1[k][j] == INF:
                    continue
                m = min(m, e.ml_c() * (k - i) + M1[k][j])
                if k - 1 >= i and M[i][k - 1] < INF:
                    m = min(m, M[i][k - 1] + M1[k][j])
            M[i][j] = m

            best_w = w(i, j - 1)
            for k in range(i, j + 1):
                if V[k][j] < INF:
                    best_w = min(best_w, w(i, k - 1) + V[k][j])
            W[i][j] = min(best_w, 0.0)
    return V, M1, M, W


def zuker_mfe(seq: str, energies=None, max_loop: int = 30) -> tuple[float, str]:
    """(MFE kcal/mol, nokta-parantez yapisi)."""
    s = validate(seq)
    e = DnaLoopEnergies() if energies is None else energies
    n = len(s)
    if n < MIN_HAIRPIN + 2:
        return 0.0, "." * n
    V, M1, M, W = _fill_mfe(s, e, max_loop)
    out = ["."] * n
    _tb_w(s, V, M1, M, W, 0, n - 1, e, out, max_loop)
    return W[0][n - 1], "".join(out)


def _w(W, i, j):
    return 0.0 if i > j else W[i][j]


def _tb_w(s, V, M1, M, W, i, j, e, out, ml):
    if i >= j:
        return
    t = W[i][j]
    if abs(t - _w(W, i, j - 1)) < 1e-9:
        _tb_w(s, V, M1, M, W, i, j - 1, e, out, ml)
        return
    for k in range(i, j + 1):
        if V[k][j] < INF and abs(t - (_w(W, i, k - 1) + V[k][j])) < 1e-9:
            out[k], out[j] = "(", ")"
            _tb_v(s, V, M1, M, W, k, j, e, out, ml)
            _tb_w(s, V, M1, M, W, i, k - 1, e, out, ml)
            return


def _tb_v(s, V, M1, M, W, i, j, e, out, ml):
    t = V[i][j]
    if abs(t - _hairpin_pen(e, s, i, j)) < 1e-9:
        return
    for p, q in _inner_pairs(len(s), i, j, ml):
        if V[p][q] < INF and abs(t - (_inner_pen(e, s, i, j, p, q) + V[p][q])) < 1e-9:
            out[p], out[q] = "(", ")"
            _tb_v(s, V, M1, M, W, p, q, e, out, ml)
            return
    for k in range(i + 2, j):
        if M[i + 1][k - 1] < INF and M1[k][j - 1] < INF and abs(
            t - (e.ml_a() + e.ml_b() + M[i + 1][k - 1] + M1[k][j - 1])
        ) < 1e-9:
            _tb_m(s, V, M1, M, W, i + 1, k - 1, e, out, ml)
            _tb_m1(s, V, M1, M, W, k, j - 1, e, out, ml)
            return


def _tb_m1(s, V, M1, M, W, i, j, e, out, ml):
    t = M1[i][j]
    for k in range(i + MIN_HAIRPIN + 1, j + 1):
        if V[i][k] < INF and abs(t - (V[i][k] + e.ml_b() + e.ml_c() * (j - k))) < 1e-9:
            out[i], out[k] = "(", ")"
            _tb_v(s, V, M1, M, W, i, k, e, out, ml)
            return


def _tb_m(s, V, M1, M, W, i, j, e, out, ml):
    t = M[i][j]
    for k in range(i, j + 1):
        if M1[k][j] == INF:
            continue
        if abs(t - (e.ml_c() * (k - i) + M1[k][j])) < 1e-9:
            _tb_m1(s, V, M1, M, W, k, j, e, out, ml)
            return
        if k - 1 >= i and M[i][k - 1] < INF and \
           abs(t - (M[i][k - 1] + M1[k][j])) < 1e-9:
            _tb_m(s, V, M1, M, W, i, k - 1, e, out, ml)
            _tb_m1(s, V, M1, M, W, k, j, e, out, ml)
            return


# ---------------------------------------------------------------------------
# McCaskill: ic (inside) + dis (outside)
# ---------------------------------------------------------------------------
def _inside(s, e, rt, max_loop):
    n = len(s)
    bz = lambda g: math.exp(-g / rt)  # noqa: E731
    Zb = [[0.0] * n for _ in range(n)]
    Z1 = [[0.0] * n for _ in range(n)]
    Zm = [[0.0] * n for _ in range(n)]
    Z = [[1.0] * (n + 1) for _ in range(n + 1)]

    def z(i, j):
        return 1.0 if i > j else Z[i][j]

    for span in range(0, n):
        for i in range(0, n - span):
            j = i + span
            if span > MIN_HAIRPIN and can_pair(s[i], s[j]):
                tot = bz(_hairpin_pen(e, s, i, j))
                for p, q in _inner_pairs(n, i, j, max_loop):
                    if Zb[p][q] > 0.0:
                        tot += bz(_inner_pen(e, s, i, j, p, q)) * Zb[p][q]
                for k in range(i + 2, j):
                    tot += bz(e.ml_a() + e.ml_b()) * Zm[i + 1][k - 1] * Z1[k][j - 1]
                Zb[i][j] = tot

            acc1 = 0.0
            for k in range(i + MIN_HAIRPIN + 1, j + 1):
                if Zb[i][k] > 0.0:
                    acc1 += Zb[i][k] * bz(e.ml_b() + e.ml_c() * (j - k))
            Z1[i][j] = acc1

            accm = 0.0
            for k in range(i, j + 1):
                if Z1[k][j] == 0.0:
                    continue
                accm += bz(e.ml_c() * (k - i)) * Z1[k][j]
                if k - 1 >= i:
                    accm += Zm[i][k - 1] * Z1[k][j]
            Zm[i][j] = accm

            acc = z(i, j - 1)
            for k in range(i, j + 1):
                if Zb[k][j] > 0.0:
                    acc += z(i, k - 1) * Zb[k][j]
            Z[i][j] = acc
    return Zb, Z1, Zm, Z, z


def mccaskill(seq: str, energies=None, temp_c: float = 37.0,
              max_loop: int = 30) -> tuple[float, list[list[float]]]:
    """(Z, baz-cifti olasilik matrisi). Dis gecis tam yapilir."""
    s = validate(seq)
    e = DnaLoopEnergies() if energies is None else energies
    if isinstance(e, DnaLoopEnergies) and temp_c != 37.0:
        raise ValueError("DNA uzunluk modeli yalniz 37 C icin tanimlidir")
    n = len(s)
    rt = R_CAL * (temp_c + KELVIN) / 1000.0
    bz = lambda g: math.exp(-g / rt)  # noqa: E731
    if n < MIN_HAIRPIN + 2:
        return 1.0, [[0.0] * n for _ in range(n)]

    Zb, Z1, Zm, Z, z = _inside(s, e, rt, max_loop)
    total = Z[0][n - 1]

    Zh = [[0.0] * n for _ in range(n)]   # outside
    for span in range(n - 1, MIN_HAIRPIN, -1):
        for i in range(0, n - span):
            j = i + span
            if Zb[i][j] == 0.0:
                continue
            out = z(0, i - 1) * z(j + 1, n - 1)
            for p in range(0, i):
                for q in range(j + 1, n):
                    if Zh[p][q] == 0.0:
                        continue
                    s1, s2 = i - p - 1, q - j - 1
                    if s1 + s2 <= max_loop:
                        out += bz(_inner_pen(e, s, p, q, i, j)) * Zh[p][q]
                    left = Zm[p + 1][i - 1] if i - 1 >= p + 1 else 0.0
                    right = Zm[j + 1][q - 1] if q - 1 >= j + 1 else 0.0
                    cl = bz(e.ml_c() * s1)
                    cr = bz(e.ml_c() * s2)
                    ways = left * cr + cl * right + left * right
                    if ways > 0.0:
                        out += Zh[p][q] * bz(e.ml_a() + 2 * e.ml_b()) * ways
            Zh[i][j] = out

    probs = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if Zb[i][j] > 0.0 and total > 0.0:
                p = Zb[i][j] * Zh[i][j] / total
                probs[i][j] = probs[j][i] = min(1.0, max(0.0, p))
    return total, probs


def unpaired_probability(seq: str, energies=None, temp_c: float = 37.0) -> list[float]:
    """Her konumun eslesmemis olma olasiligi (erisilebilirlik, SS-03)."""
    _, probs = mccaskill(seq, energies, temp_c)
    return [max(0.0, 1.0 - sum(row)) for row in probs]


def ensemble_free_energy(seq: str, energies=None, temp_c: float = 37.0) -> float:
    """-RT ln Z (kcal/mol)."""
    total, _ = mccaskill(seq, energies, temp_c)
    rt = R_CAL * (temp_c + KELVIN) / 1000.0
    return -rt * math.log(total)


# ---------------------------------------------------------------------------
# Kaba kuvvet dogrulama
# ---------------------------------------------------------------------------
def enumerate_structures(seq: str) -> Iterator[list[tuple[int, int]]]:
    """Tum pseudoknot'suz ikincil yapilar (kucuk n icin)."""
    s = validate(seq)
    n = len(s)

    @lru_cache(maxsize=None)
    def rec(i: int, j: int) -> tuple[tuple[tuple[int, int], ...], ...]:
        if j - i < MIN_HAIRPIN + 1:
            return ((),)
        out = list(rec(i + 1, j))
        for k in range(i + MIN_HAIRPIN + 1, j + 1):
            if not can_pair(s[i], s[k]):
                continue
            for inner in rec(i + 1, k - 1):
                for rest in rec(k + 1, j):
                    out.append(tuple(sorted(((i, k),) + inner + rest)))
        return tuple(dict.fromkeys(out))

    for st in rec(0, n - 1):
        yield list(st)


def structure_energy(seq: str, pairs, energies) -> float:
    """Bir yapinin enerjisi; DP ile ayni ayristirma."""
    s = validate(seq)
    if not pairs:
        return 0.0
    children: dict[tuple[int, int], list[tuple[int, int]]] = {}
    stack: list[tuple[int, int]] = []
    for (i, j) in sorted(pairs):
        while stack and not (stack[-1][0] < i and j < stack[-1][1]):
            stack.pop()
        if stack:
            children.setdefault(stack[-1], []).append((i, j))
        children.setdefault((i, j), [])
        stack.append((i, j))

    total = 0.0
    for (i, j), kids in children.items():
        if not kids:
            total += _hairpin_pen(energies, s, i, j)
        elif len(kids) == 1:
            k, m = kids[0]
            total += _inner_pen(energies, s, i, j, k, m)
        else:
            unpaired = (j - i - 1) - sum(m - k + 1 for k, m in kids)
            total += (energies.ml_a()
                      + energies.ml_b() * (1 + len(kids))
                      + energies.ml_c() * unpaired)
    return total
