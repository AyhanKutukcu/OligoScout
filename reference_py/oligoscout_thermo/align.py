"""Termodinamik hizalama: iki oligo arasindaki en kararli dupleks (dimer).

QC karsiligi: DIM-01..DIM-04, DIM-06 ve DIM-08.

Model ve sinirlari (acikca beyan edilir):
  * Taranan uzay : iki oligonun tum goreli kaymalari (offset).
  * Puanlanan yapi: kesintisiz Watson-Crick helisi, ya da EN FAZLA BIR kusur:
      - tek bir ic uyumsuzluk (her iki yaninda >= 2 WC cifti), veya
      - tek bazlik bir ic bulge (her iki yaninda >= 2 WC cifti).
    Coklu uyumsuzluk, birden buyuk bulge, ayni yapida uyumsuzluk + bulge ve
    diziye bagli uc terimleri kapsam disindadir.
  * Enerji : iki uc baslatmasi ve NN yiginlari, 37 C ve 1 M NaCl'de.
  * Mod    : "any"      -> en dusuk dG'li helis
             "3prime"   -> en az bir oligonun 3' terminal bazi eslesmis olmali
                           (polimerazla uzatilabilir dimer)

Tek bazlik bulge enerjisi (SantaLucia & Hicks 2004, denklem 13):
BULGE:1 + kesintisiz iplikteki komsu WC yigini + her A-T kapanisina
0.5 kcal/mol. Uzunluk tablosu yalniz dG37 tasir ve
kaynak varsayimi dH = 0'dir; bu yuzden bulge terimi dH = 0,
dS = -1000 * dG37 / 310.15 olarak dagitilir. Yigin terimi dH/dS olarak
dogrudan WC tablosundan gelir.

`allow_bulge=False` ile bulge kapatilir; o zaman model, bulge oncesi
surumle birebir ayni sonucu verir (karsilastirma: cli/bulge_compare.py).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .nn import _mismatch_stack_adjustment, reverse_complement, validate
from .params import (
    AT_BULGE_CLOSURE,
    DNA_LOOP_SANTALUCIA_HICKS_2004,
    KELVIN,
    MISMATCH_ALLAWI_PEYRET,
    R_CAL,
    WC_SANTALUCIA_1998,
)

_PAIR = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
MODES = ("any", "3prime")
MODEL_NAME = "nn_single_internal_mismatch_or_bulge"
MODEL_NAME_NO_BULGE = "nn_single_internal_mismatch"
_T37 = 310.15
# Beyan edilmis model sabiti; deger ve gerekce params.MODEL_CONSTANTS'ta
# tutulur ve registry_digest'e girer (REP-02).
_AT_BULGE_CLOSURE_DG37 = AT_BULGE_CLOSURE.value


def model_metadata(allow_bulge: bool = True) -> dict:
    provenance = MISMATCH_ALLAWI_PEYRET.provenance
    loop = DNA_LOOP_SANTALUCIA_HICKS_2004.provenance
    meta = {
        "model": MODEL_NAME if allow_bulge else MODEL_NAME_NO_BULGE,
        "mismatch_table": provenance.table_id,
        "mismatch_table_status": provenance.verification_status,
        "reference_conditions": "37 C, 1 M NaCl",
        "limitation": (
            "en fazla bir kusur: tekli ic uyumsuzluk VEYA tek bazlik ic "
            "bulge (her iki yanda >=2 WC); coklu kusur ve uc terimleri kapsam disi"
            if allow_bulge else
            "tekli ic uyumsuzluk; bulge ve uc terimleri kapsam disi"
        ),
    }
    if allow_bulge:
        meta["loop_table"] = loop.table_id
        meta["loop_table_status"] = loop.verification_status
        meta["model_constants"] = [AT_BULGE_CLOSURE.as_dict()]
    return meta


@dataclass(frozen=True)
class DimerResult:
    """Bir oligo ciftinin en kararli (puanlanabilir) dimeri."""

    dH_kcal: float
    dS_cal: float
    dG37_kcal: float
    offset: int            # b'nin a'ya gore kaymasi
    helix_len: int         # eslesen baz cifti sayisi
    mode: str
    top: str               # hizalamanin ust satiri (a, 5'->3')
    match: str             # '|' eslesme, 'X' ic uyumsuzluk, 'B' bulge
    bottom: str            # alt satir (b, 3'->5')
    three_prime_paired: bool
    found: bool = True
    mismatch_count: int = 0
    bulge_count: int = 0
    bulge_strand: str = ""   # "" | "a" | "b" -- bulge hangi ipliktedir

    def tm_c(self, total_strand_M: float) -> float:
        """Iki molekullu, esit konsantrasyon modeliyle dimer Tm'i (C)."""
        if not self.found or self.helix_len < 2:
            return float("-inf")
        denom = self.dS_cal + R_CAL * math.log(total_strand_M / 4.0)
        if denom >= 0:
            return float("-inf")
        return self.dH_kcal * 1000.0 / denom - KELVIN

    def as_dict(self, total_strand_M: float = 2.0e-7) -> dict:
        return {
            "dG37_kcal_mol": round(self.dG37_kcal, 4),
            "dH_kcal_mol": round(self.dH_kcal, 4),
            "dS_cal_mol_K": round(self.dS_cal, 4),
            "tm_c": round(self.tm_c(total_strand_M), 3),
            "helix_len": self.helix_len,
            "offset": self.offset,
            "mode": self.mode,
            "three_prime_paired": self.three_prime_paired,
            "alignment": [self.top, self.match, self.bottom],
            "mismatch_count": self.mismatch_count,
            "bulge_count": self.bulge_count,
            "bulge_strand": self.bulge_strand,
        }


_EMPTY = DimerResult(0.0, 0.0, 0.0, 0, 0, "any", "", "", "", False, found=False)


def _nn_step(step: str) -> tuple[float, float]:
    """Bir din. yiginin (dH, dS)'i; anahtar yoksa ters tamamlayicisi alinir."""
    nn = WC_SANTALUCIA_1998.require()["nn"]
    h, s = nn.get(step) or nn[reverse_complement(step)]
    return h, s


def _stack_energy(top_segment: str) -> tuple[float, float]:
    """Yalniz NN yiginlari; uc baslatmasi EKLENMEZ."""
    dH = dS = 0.0
    for i in range(len(top_segment) - 1):
        h, s = _nn_step(top_segment[i : i + 2])
        dH += h
        dS += s
    return dH, dS


def _init_energy(base: str) -> tuple[float, float]:
    init = WC_SANTALUCIA_1998.require()["init"]
    return init["AT" if base in "AT" else "GC"]


def _helix_energy(top_segment: str) -> tuple[float, float]:
    """Kesintisiz, tam eslesen helisin (dH, dS)'i. Ust iplik 5'->3'."""
    dH, dS = _stack_energy(top_segment)
    for base in (top_segment[0], top_segment[-1]):
        h, s = _init_energy(base)
        dH += h
        dS += s
    return dH, dS


def _bulge1_energy() -> tuple[float, float]:
    """Tek bazlik bulge uzunluk teriminin (dH, dS)'i.

    Tablo yalniz dG37 tasir ve kaynak varsayimi dH = 0'dir (bkz. kunye);
    terim bu yuzden tamamen entropiye yuklenir.
    """
    dg37 = DNA_LOOP_SANTALUCIA_HICKS_2004.require()["BULGE:1"][0]
    return 0.0, -1000.0 * dg37 / _T37


def best_dimer(a: str, b: str, mode: str = "any",
               allow_bulge: bool = True) -> DimerResult:
    """a ve b arasinda en dusuk dG37'li, en fazla bir kusurlu dupleks.

    `a` ve `b` 5'->3' verilir. Hizalamada b ters cevrilir (3'->5'),
    boylece ayni sutunda karsilikli bazlar bulunur.

    `allow_bulge=False` tek bazlik bulge taramasini kapatir.
    """
    if mode not in MODES:
        raise ValueError(f"bilinmeyen mod: {mode}")
    sa, sb = validate(a), validate(b)
    rb = sb[::-1]                      # b, 3'->5'
    best = _scan_contiguous(sa, rb, mode)
    if allow_bulge:
        bulged = _scan_bulge(sa, rb, mode)
        if bulged.found and (not best.found or bulged.dG37_kcal < best.dG37_kcal):
            best = bulged
    return best


def _scan_contiguous(sa: str, rb: str, mode: str) -> DimerResult:
    """Kesintisiz helis ya da tek ic uyumsuzluk; tum kaymalar taranir."""
    na, nb = len(sa), len(rb)
    best = _EMPTY

    for offset in range(-(nb - 1), na):
        # sutun k: sa[k] ile rb[k - offset]
        lo = max(0, offset)
        hi = min(na, offset + nb)
        if hi - lo < 2:
            continue
        paired = [(sa[k], rb[k - offset]) in _PAIR for k in range(lo, hi)]

        start = previous_run = None
        for idx in range(len(paired) + 1):
            if idx < len(paired) and paired[idx]:
                if start is None:
                    start = idx
                continue
            if start is None:
                continue
            run_lo, run_hi = lo + start, lo + idx      # [run_lo, run_hi)
            start = None
            candidates = []
            if run_hi - run_lo >= 2:
                candidates.append((run_lo, None))
                if previous_run:
                    prev_lo, prev_hi = previous_run
                    if run_lo == prev_hi + 1 and prev_hi - prev_lo >= 2:
                        candidates.append((prev_lo, prev_hi))
            previous_run = (run_lo, run_hi)

            # Ek WC yigini, uc baslatma farkindan daha cok stabilizasyon saglar.
            for left, mismatch in candidates:
                terminal = run_hi == na or left == offset
                if mode == "3prime" and not terminal:
                    continue
                top_segment = sa[left:run_hi]
                dH, dS = _helix_energy(top_segment)
                if mismatch is not None:
                    bottom_segment = rb[left - offset:run_hi - offset]
                    delta_h, delta_s = _mismatch_stack_adjustment(top_segment, bottom_segment,
                                                                mismatch - left)
                    dH += delta_h
                    dS += delta_s
                dG = dH - 310.15 * dS / 1000.0
                if best.found and dG >= best.dG37_kcal:
                    continue
                pad = max(0, -offset)
                top_line = " " * pad + sa
                bot_line = " " * max(0, offset) + rb
                width = max(len(top_line), len(bot_line))
                marks = [" "] * width
                marks[pad + left:pad + run_hi] = "|" * (run_hi - left)
                count = int(mismatch is not None)
                if mismatch is not None:
                    marks[pad + mismatch] = "X"
                best = DimerResult(
                    dH_kcal=dH, dS_cal=dS, dG37_kcal=dG, offset=offset,
                    helix_len=run_hi - left - count, mode=mode,
                    top=top_line.ljust(width), match="".join(marks), bottom=bot_line.ljust(width),
                    three_prime_paired=terminal, mismatch_count=count,
                )
    return best


def _scan_bulge(sa: str, rb: str, mode: str) -> DimerResult:
    """Tek bazlik ic bulge tasiyan en kararli dupleks.

    Geometri: sol segment `offset` kaymasinda, sag segment `offset -+ 1`
    kaymasindadir. Fark 1 oldugu icin tam olarak bir baz eslesmeden kalir:
    `strand == "b"` ise bu baz b ipliginde, `strand == "a"` ise a
    ipligindedir. Her iki segment de en az iki Watson-Crick cifti
    icermelidir.

    Her iki segment de BULGE'DAN UZAGA dogru buyutulur (sol segment geriye,
    sag segment ileriye). Bu sozlesme aynali simetriyi garanti eder:
    `best_dimer(a, b)` ile `best_dimer(b, a)` ayni enerjiyi verir, cunku
    iplikleri degistirmek hizalamayi yatay olarak aynalar ve sol/sag
    segmentlerin rollerini takas eder.
    """
    na, nb = len(sa), len(rb)
    best = _EMPTY
    bulge_h, bulge_s = _bulge1_energy()

    for offset in range(-(nb - 1), na):
        lo, hi = max(0, offset), min(na, offset + nb)
        for strand in ("b", "a"):
            o2 = offset - 1 if strand == "b" else offset + 1
            lo2, hi2 = max(0, o2), min(na, o2 + nb)
            for h1 in range(lo + 2, hi + 1):
                if (sa[h1 - 1], rb[h1 - 1 - offset]) not in _PAIR:
                    continue
                left = h1 - 1
                while left - 1 >= lo and (sa[left - 1],
                                          rb[left - 1 - offset]) in _PAIR:
                    left -= 1
                if h1 - left < 2:
                    continue
                l2 = h1 if strand == "b" else h1 + 1
                if l2 < lo2 or l2 >= hi2:
                    continue
                if (sa[l2], rb[l2 - o2]) not in _PAIR:
                    continue
                h2 = l2
                while h2 < hi2 and (sa[h2], rb[h2 - o2]) in _PAIR:
                    h2 += 1
                if h2 - l2 < 2:
                    continue

                dH, dS = 0.0, 0.0
                for segment in (sa[left:h1], sa[l2:h2]):
                    sh, ss = _stack_energy(segment)
                    dH += sh
                    dS += ss
                for base in (sa[left], sa[h2 - 1]):
                    ih, isx = _init_energy(base)
                    dH += ih
                    dS += isx
                dH += bulge_h
                dS += bulge_s
                at_closures = (sa[h1 - 1] in "AT") + (sa[l2] in "AT")
                dS -= 1000.0 * _AT_BULGE_CLOSURE_DG37 * at_closures / _T37
                # Denklem 13: bulge'un uzerinden gecen komsu WC yigini,
                # kesintisiz iplikte okunur.
                if strand == "b":
                    sh, ss = _nn_step(sa[h1 - 1] + sa[h1])
                else:
                    j = h1 - 1 - offset
                    sh, ss = _nn_step(rb[j + 1] + rb[j])
                dH += sh
                dS += ss
                dG = dH - _T37 * dS / 1000.0

                terminal = h2 == na or left == offset
                if mode == "3prime" and not terminal:
                    continue
                if best.found and dG >= best.dG37_kcal:
                    continue

                pad = max(0, -offset)
                if strand == "b":
                    top_line = " " * pad + sa[:h1] + "-" + sa[h1:]
                    bot_line = " " * max(0, offset) + rb
                    right_shift = 1
                else:
                    m = h1 - offset
                    top_line = " " * pad + sa
                    bot_line = " " * max(0, offset) + rb[:m] + "-" + rb[m:]
                    right_shift = 0
                width = max(len(top_line), len(bot_line))
                marks = [" "] * width
                for k in range(left, h1):
                    marks[pad + k] = "|"
                for k in range(l2, h2):
                    marks[pad + k + right_shift] = "|"
                marks[pad + h1] = "B"

                best = DimerResult(
                    dH_kcal=dH, dS_cal=dS, dG37_kcal=dG, offset=offset,
                    helix_len=(h1 - left) + (h2 - l2), mode=mode,
                    top=top_line.ljust(width), match="".join(marks),
                    bottom=bot_line.ljust(width), three_prime_paired=terminal,
                    bulge_count=1, bulge_strand=strand,
                )
    return best


def self_dimer(a: str, mode: str = "any", allow_bulge: bool = True) -> DimerResult:
    return best_dimer(a, a, mode, allow_bulge=allow_bulge)


def hairpin_stem(a: str, min_loop: int = 3, mode: str = "any") -> DimerResult:
    """Tek molekullu sap (stem) taramasi — yalnizca yigin enerjisi.

    Bu tarama yalniz sap enerjisini verir. SS-01 icin `fold.zuker_mfe`
    dongu terimlerini de uygular.
    """
    s = validate(a)
    best = _EMPTY
    n = len(s)
    for i in range(n):
        for j in range(i + min_loop + 1, n):
            length = 0
            while (
                i - length >= 0
                and j + length < n
                and (s[i - length], s[j + length]) in _PAIR
            ):
                length += 1
            if length < 2:
                continue
            seg = s[i - length + 1 : i + 1]
            dH, dS = _helix_energy(seg)
            dG = dH - 310.15 * dS / 1000.0
            if not best.found or dG < best.dG37_kcal:
                best = DimerResult(
                    dH_kcal=dH, dS_cal=dS, dG37_kcal=dG, offset=i,
                    helix_len=length, mode=mode, top=seg, match="|" * length,
                    bottom=seg[::-1], three_prime_paired=(j + length - 1 == n - 1),
                )
    return best
