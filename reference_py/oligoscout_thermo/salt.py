"""Tuz duzeltmeleri.

Uc model, hepsi 1/Tm uzerinde calisir; hangisinin kullanildigi her rapor
satirinda beyan edilir (TH-01, TH-02, COND-01, COND-02).

  owczarzy_2004_monovalent  Tek degerlikli katyon (Na+/K+), GC bagimli.
  owczarzy_2008_divalent    Mg2+ varliginda; tek/iki degerlikli oran kurali.
  von_ahsen_2001_na_eq      Mg2+'yi Na-esdegerine cevirip tek degerlikli
                            modeli uygular; serbest Mg = [Mg] - [dNTP].

Not: Bu formullerin katsayilari da parametre tablolari gibi birincil
makalelerden aktarilmistir ve ayni dogrulama kurallarina tabidir
(bkz. FORMULA_PROVENANCE). Hicbiri henuz birincil PDF'e karsi ikinci kez
aktarilip diff alinmadi.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

FORMULA_PROVENANCE: dict[str, dict[str, str]] = {
    "owczarzy_2004_monovalent": {
        "citation": "Owczarzy et al., Biochemistry 2004;43:3537-3554",
        "form": "1/Tm(Na) = 1/Tm(1M) + (4.29*fGC - 3.95)e-5*ln[Na] + 9.40e-6*(ln[Na])^2",
        "verification_status": "transcribed_unverified",
    },
    "owczarzy_2008_divalent": {
        "citation": "Owczarzy et al., Biochemistry 2008;47:5336-5353",
        "form": (
            "1/Tm = 1/Tm(1M Na) + a + b*ln[Mg] + fGC*(c + d*ln[Mg]) "
            "+ 1/(2(N-1))*(e + f*ln[Mg] + g*(ln[Mg])^2); "
            "R = sqrt([Mg])/[Mon] kuralina gore a, d, g duzeltilir"
        ),
        "verification_status": "transcribed_unverified",
    },
    "von_ahsen_2001_na_eq": {
        "citation": "von Ahsen, Wittwer & Schutz, Clin Chem 2001;47:1956-1961",
        "form": "[Na]esdeger(mM) = [Mon](mM) + 120*sqrt([Mg]-[dNTP] in mM)",
        "verification_status": "transcribed_unverified",
    },
}

# Owczarzy 2008 katsayilari (denklem 16)
_A, _B, _C, _D, _E, _F, _G = (
    3.92e-5,
    -9.11e-6,
    6.26e-5,
    1.42e-5,
    -4.82e-4,
    5.25e-4,
    8.31e-5,
)


@dataclass(frozen=True)
class Conditions:
    """Reaksiyon kosullari (COND-01). Konsantrasyonlar molar."""

    monovalent_M: float = 0.050      # Na+ + K+
    magnesium_M: float = 0.0015
    dntp_M: float = 0.0006
    primer_M: float = 2.0e-7
    template_M: float = 0.0
    anneal_c: float = 60.0

    @property
    def free_magnesium_M(self) -> float:
        """Serbest Mg2+ (COND-02): dNTP'ler Mg'yi sekestre eder."""
        return max(0.0, self.magnesium_M - self.dntp_M)

    def as_dict(self) -> dict[str, float]:
        return {
            "monovalent_mM": round(self.monovalent_M * 1e3, 4),
            "magnesium_mM": round(self.magnesium_M * 1e3, 4),
            "free_magnesium_mM": round(self.free_magnesium_M * 1e3, 4),
            "dntp_mM": round(self.dntp_M * 1e3, 4),
            "primer_nM": round(self.primer_M * 1e9, 4),
            "anneal_c": self.anneal_c,
        }


def owczarzy_2004_monovalent(tm_1m_K: float, na_M: float, f_gc: float) -> float:
    """Tek degerlikli katyon duzeltmesi; Tm Kelvin."""
    if na_M <= 0:
        raise ValueError("[Na+] pozitif olmali")
    ln_na = math.log(na_M)
    inv = 1.0 / tm_1m_K + (4.29 * f_gc - 3.95) * 1e-5 * ln_na + 9.40e-6 * ln_na**2
    return 1.0 / inv


def von_ahsen_2001_na_eq(mon_M: float, mg_M: float, dntp_M: float) -> float:
    """Mg2+ -> Na esdegeri. Girdi ve cikti molar.

    Formulun katsayisi (120) milimolar birim icindir:
        [Na]esdeger(mM) = [Mon](mM) + 120 * sqrt([Mg] - [dNTP] in mM)
    Bu yuzden icerde mM'ye cevrilir. Serbest Mg negatifse 0 alinir (COND-02).
    """
    free_mg_mM = max(0.0, (mg_M - dntp_M) * 1e3)
    return mon_M + 120.0 * math.sqrt(free_mg_mM) / 1e3


def owczarzy_2008_divalent(
    tm_1m_K: float, mg_M: float, mon_M: float, f_gc: float, n_bp: int
) -> tuple[float, str]:
    """Mg2+ duzeltmesi. (Tm Kelvin, uygulanan rejim) dondurur."""
    if mg_M <= 0:
        if mon_M <= 0:
            raise ValueError("ne Mg2+ ne tek degerlikli katyon var")
        return owczarzy_2004_monovalent(tm_1m_K, mon_M, f_gc), "monovalent_only"

    ratio = math.sqrt(mg_M) / mon_M if mon_M > 0 else float("inf")
    if ratio < 0.22:
        return owczarzy_2004_monovalent(tm_1m_K, mon_M, f_gc), "monovalent_dominant"

    a, d, g = _A, _D, _G
    if ratio < 6.0:
        ln_mon = math.log(mon_M)
        a = _A * (0.843 - 0.352 * math.sqrt(mon_M) * ln_mon)
        d = _D * (1.279 - 4.03e-3 * ln_mon - 8.03e-3 * ln_mon**2)
        g = _G * (0.486 - 0.258 * ln_mon + 5.25e-3 * ln_mon**3)
        regime = "mixed_mono_divalent"
    else:
        regime = "divalent_dominant"

    ln_mg = math.log(mg_M)
    inv = (
        1.0 / tm_1m_K
        + a
        + _B * ln_mg
        + f_gc * (_C + d * ln_mg)
        + (1.0 / (2.0 * (n_bp - 1))) * (_E + _F * ln_mg + g * ln_mg**2)
    )
    return 1.0 / inv, regime


def correct(
    tm_1m_K: float,
    f_gc: float,
    n_bp: int,
    cond: Conditions,
    model: str = "owczarzy_2008",
) -> tuple[float, dict[str, str | float]]:
    """Secilen modele gore duzeltilmis Tm (Kelvin) ve beyan kaydi."""
    if model == "owczarzy_2004":
        tm = owczarzy_2004_monovalent(tm_1m_K, cond.monovalent_M, f_gc)
        detail: dict[str, str | float] = {"regime": "monovalent"}
    elif model == "owczarzy_2008":
        tm, regime = owczarzy_2008_divalent(
            tm_1m_K, cond.free_magnesium_M, cond.monovalent_M, f_gc, n_bp
        )
        detail = {"regime": regime}
    elif model == "von_ahsen_2001":
        na_eq = von_ahsen_2001_na_eq(cond.monovalent_M, cond.magnesium_M, cond.dntp_M)
        tm = owczarzy_2004_monovalent(tm_1m_K, na_eq, f_gc)
        detail = {"regime": "na_equivalent", "na_equivalent_M": round(na_eq, 6)}
    else:
        raise ValueError(f"bilinmeyen tuz modeli: {model}")

    detail["model"] = model
    detail["citation"] = FORMULA_PROVENANCE[
        {"owczarzy_2004": "owczarzy_2004_monovalent",
         "owczarzy_2008": "owczarzy_2008_divalent",
         "von_ahsen_2001": "von_ahsen_2001_na_eq"}[model]
    ]["citation"]
    return tm, detail


MODELS = ("owczarzy_2004", "owczarzy_2008", "von_ahsen_2001")


def applicable_models(cond: Conditions) -> tuple[str, ...]:
    """Beyan edilen kosullara uygulanabilir tuz modelleri.

    Mg2+ varken `owczarzy_2004` modeli Mg'yi hic gormez; onu yayilima katmak
    model belirsizligini degil *kosul* farkini olcer. Bu yuzden serbest Mg > 0
    oldugunda yalniz Mg'yi hesaba katan iki model karsilastirilir. Mg yokken
    zaten owczarzy_2008 ve von_ahsen_2001 tek degerlikli rejime duser.
    """
    if cond.free_magnesium_M > 0:
        return ("owczarzy_2008", "von_ahsen_2001")
    return MODELS


def spread(tm_1m_K: float, f_gc: float, n_bp: int, cond: Conditions) -> dict[str, float]:
    """TH-02: uygulanabilir tuz modelleri altinda Tm yayilimi (Celsius)."""
    out: dict[str, float] = {}
    for m in MODELS:
        tm, _ = correct(tm_1m_K, f_gc, n_bp, cond, m)
        out[m] = tm - 273.15
    used = applicable_models(cond)
    vals = [out[m] for m in used]
    out["spread_c"] = max(vals) - min(vals)
    out["n_models_compared"] = float(len(used))
    return out
