"""Nearest-neighbor duplex termodinamigi (Watson-Crick).

Model:  dH_toplam = sum(NN yiginlari) + baslatma(5' uc) + baslatma(3' uc)
        dS_toplam = ayni + simetri duzeltmesi (kendine tamamlayici ise)
        dG(T)     = dH - T * dS

Birimler: dH kcal/mol, dS cal/(mol K), dG kcal/mol, T Kelvin.
"""

from __future__ import annotations

from dataclasses import dataclass

from .params import (
    KELVIN,
    MISMATCH_ALLAWI_PEYRET,
    TERMINAL_MISMATCH_SANTALUCIA_PEYRET,
    WC_SANTALUCIA_1998,
)

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
VALID = set("ACGT")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def validate(seq: str) -> str:
    """Diziyi buyuk harfe cevirir ve alfabe disinda karakter varsa hata verir (SEQ-01)."""
    s = seq.strip().upper()
    if not s:
        raise ValueError("bos dizi")
    bad = sorted(set(s) - VALID)
    if bad:
        raise ValueError(f"tanimsiz karakter: {', '.join(bad)} (yalniz A/C/G/T)")
    return s


def gc_fraction(seq: str) -> float:
    s = validate(seq)
    return (s.count("G") + s.count("C")) / len(s)


def is_self_complementary(seq: str) -> bool:
    s = validate(seq)
    return s == reverse_complement(s)


@dataclass(frozen=True)
class DuplexThermo:
    """Bir dupleksin standart durum termodinamigi (1 M NaCl)."""

    dH_kcal: float
    dS_cal: float
    self_complementary: bool
    n_bp: int

    def dG_kcal(self, temp_c: float = 37.0) -> float:
        t = temp_c + KELVIN
        return self.dH_kcal - t * self.dS_cal / 1000.0

    def as_dict(self, temp_c: float = 37.0) -> dict[str, float | bool | int]:
        return {
            "dH_kcal_mol": round(self.dH_kcal, 4),
            "dS_cal_mol_K": round(self.dS_cal, 4),
            f"dG{int(temp_c)}_kcal_mol": round(self.dG_kcal(temp_c), 4),
            "self_complementary": self.self_complementary,
            "n_bp": self.n_bp,
        }


def _init_key(base: str) -> str:
    return "AT" if base in "AT" else "GC"


def _wc_stack(stacks, step: str) -> list[float]:
    return stacks[step] if step in stacks else stacks[reverse_complement(step)]


def _mismatch_stack_adjustment(top: str, bottom: str, position: int) -> tuple[float, float]:
    """Tek uyumsuzlugun iki komsu yigini icin WC'ye gore dH/dS farki."""
    stacks = MISMATCH_ALLAWI_PEYRET.require()
    wc = WC_SANTALUCIA_1998.require()["nn"]
    delta_h = delta_s = 0.0
    for i in (position - 1, position):
        top_step, bottom_step = top[i:i + 2], bottom[i:i + 2]
        key = min(f"{top_step}/{bottom_step}", f"{bottom_step[::-1]}/{top_step[::-1]}")
        wc_h, wc_s = _wc_stack(wc, top_step)
        h, s = stacks[key]
        delta_h += h - wc_h
        delta_s += s - wc_s
    return delta_h, delta_s


def watson_crick(seq: str) -> DuplexThermo:
    """Tam eslesen duplexin dH/dS'i. `seq` ust iplik, 5'->3'."""
    s = validate(seq)
    if len(s) < 2:
        raise ValueError("NN modeli icin en az 2 nt gerekir")

    table = WC_SANTALUCIA_1998.require()
    nn = table["nn"]
    init = table["init"]
    sym_dH, sym_dS = table["symmetry"]

    dH = dS = 0.0
    for i in range(len(s) - 1):
        h, sv = _wc_stack(nn, s[i : i + 2])
        dH += h
        dS += sv

    for end_base in (s[0], s[-1]):
        h, sv = init[_init_key(end_base)]
        dH += h
        dS += sv

    selfc = is_self_complementary(s)
    if selfc:
        dH += sym_dH
        dS += sym_dS

    return DuplexThermo(dH_kcal=dH, dS_cal=dS, self_complementary=selfc, n_bp=len(s))


def duplex_with_mismatches(primer: str, template_site: str) -> DuplexThermo:
    """Uyumsuzluk iceren duplex (SPC-04, VAR-04 icin gerekli).

    `primer`        : 5'->3' primer dizisi
    `template_site` : primerin bagladigi kalip bolgesinin, primerle ayni
                      yonde yazilmis tamamlayicisi (yani tam eslesmede
                      template_site == primer olur)

    Parametreler tek bir ic uyumsuzluk icindir. Uc ve uca bitisik
    konumlar kaynak tablolardaki uygulama alaninin disindadir.
    """
    p = validate(primer)
    t = validate(template_site)
    if len(p) != len(t):
        raise ValueError("primer ve kalip bolgesi ayni uzunlukta olmali")
    if p == t:
        return watson_crick(p)

    positions = [i for i, (a, b) in enumerate(zip(p, t)) if a != b]
    if len(positions) != 1:
        raise ValueError("yalniz tek uyumsuzluk desteklenir")
    position = positions[0]
    if not 2 <= position < len(p) - 2:
        raise ValueError("uc ve uca bitisik uyumsuzluk desteklenmez")

    wc = WC_SANTALUCIA_1998.require()
    baseline = watson_crick(p)
    delta_h, delta_s = _mismatch_stack_adjustment(p, t.translate(_COMPLEMENT), position)

    self_complementary = p == reverse_complement(t)
    symmetry_change = int(self_complementary) - int(baseline.self_complementary)
    sym_h, sym_s = wc["symmetry"]
    return DuplexThermo(baseline.dH_kcal + delta_h + symmetry_change * sym_h,
                        baseline.dS_cal + delta_s + symmetry_change * sym_s,
                        self_complementary, len(p))


def duplex_with_arms_mismatch(primer: str, template_site: str) -> DuplexThermo:
    """CT-04: exactly one mismatch at the primer's last or penultimate base."""
    p, t = validate(primer), validate(template_site)
    positions = mismatch_positions(p, t)
    if len(p) < 3 or len(positions) != 1 or positions[0] < len(p) - 2:
        raise ValueError("CT-04: en az 3 nt ve son iki bazda tek uyumsuzluk gerekir")
    position = positions[0]
    bottom = t.translate(_COMPLEMENT)
    if position == len(p) - 1:
        baseline = watson_crick(p[:-1])
        h, s = TERMINAL_MISMATCH_SANTALUCIA_PEYRET.require()[f"{p[-2:]}/{bottom[-2:]}"]
    else:
        baseline = watson_crick(p)
        h, s = _mismatch_stack_adjustment(p, bottom, position)
    if baseline.self_complementary:
        sym_h, sym_s = WC_SANTALUCIA_1998.require()["symmetry"]
        h -= sym_h
        s -= sym_s
    return DuplexThermo(baseline.dH_kcal + h, baseline.dS_cal + s, False, len(p))


def mismatch_positions(primer: str, template_site: str) -> list[int]:
    """Uyumsuz konumlarin 0-tabanli indeksleri (primerin 5'->3' yonunde)."""
    p = validate(primer)
    t = validate(template_site)
    if len(p) != len(t):
        raise ValueError("primer ve kalip bolgesi ayni uzunlukta olmali")
    return [i for i, (a, b) in enumerate(zip(p, t)) if a != b]
