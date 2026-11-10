"""Konsantrasyon modelleri, Tm ve baglanma kesri theta(T).

Uc model acikca secilir ve her raporda beyan edilir (QC: TH-01, TH-03, TH-05):

  equal            Kendine tamamlayici olmayan duplex, iki iplik esit
                   toplam konsantrasyonda (C_etkin = C_T/4).
  self_complementary  Kendine tamamlayici oligo (C_etkin = C_T).
  primer_excess    Primer kalibin cok uzerinde; [P] sabit kabul edilir
                   (C_etkin = [P]). PCR icin fiziksel olarak dogru olan budur.

Hepsinde Tm su ortak formdan gelir:
      Tm = dH / (dS + R ln C_etkin)
ve theta(Tm) = 0.5 ozdesligi her uc model icin analitik olarak saglanir;
tests/test_thermo.py bunu 1e-9 toleransla dogrular.

Tuz duzeltmesi modelleme karari (aciktan beyan edilir):
  Tuz duzeltmesi Tm uzerinde tanimlidir. theta(T) egrisini tutarli kilmak
  icin dH sabit tutulur ve dS, duzeltilmis Tm'i tam olarak verecek sekilde
  yeniden olculur:  dS_etkin = dH / Tm_duzeltilmis - R ln C_etkin.
  Bu sayede theta(Tm_duzeltilmis) = 0.5 ozdesligi tuz duzeltmesi altinda da
  korunur. Rapor kaydinda `ds_rescaled_to_salt_corrected_tm: true` olarak
  gorunur.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .nn import DuplexThermo, duplex_with_arms_mismatch, duplex_with_mismatches, gc_fraction, watson_crick
from .params import KELVIN, R_CAL
from .salt import Conditions, correct

MODELS = ("equal", "self_complementary", "primer_excess")


def effective_concentration(model: str, cond: Conditions, total_strand_M: float) -> float:
    if model == "equal":
        return total_strand_M / 4.0
    if model == "self_complementary":
        return total_strand_M
    if model == "primer_excess":
        return cond.primer_M
    raise ValueError(f"bilinmeyen konsantrasyon modeli: {model}")


def default_model(duplex: DuplexThermo) -> str:
    return "self_complementary" if duplex.self_complementary else "primer_excess"


@dataclass(frozen=True)
class MeltingModel:
    """Bir oligo icin tam tanimli erime modeli."""

    duplex: DuplexThermo
    cond: Conditions
    f_gc: float
    conc_model: str
    salt_model: str
    c_eff_M: float
    tm_1m_K: float
    tm_K: float
    dS_eff_cal: float
    salt_detail: dict

    # -- tureyen degerler ---------------------------------------------------
    @property
    def tm_c(self) -> float:
        return self.tm_K - KELVIN

    def dG_cal(self, temp_c: float) -> float:
        """Tuz-duzeltilmis dG(T), cal/mol."""
        t = temp_c + KELVIN
        return self.duplex.dH_kcal * 1000.0 - t * self.dS_eff_cal

    def dG_kcal(self, temp_c: float = 37.0) -> float:
        return self.dG_cal(temp_c) / 1000.0

    def equilibrium_constant(self, temp_c: float) -> float:
        t = temp_c + KELVIN
        return math.exp(-self.dG_cal(temp_c) / (R_CAL * t))

    def theta(self, temp_c: float) -> float:
        """Baglanmis iplik kesri (0..1)."""
        k = self.equilibrium_constant(temp_c)
        if self.conc_model == "primer_excess":
            kp = k * self.cond.primer_M
            return kp / (1.0 + kp)
        if self.conc_model == "equal":
            # A + B <-> AB, her iplik C_T/2; x = [AB]
            a = self.c_eff_M * 2.0  # C_T/2  (c_eff = C_T/4)
            disc = math.sqrt(max(0.0, (2 * k * a + 1) ** 2 - 4 * (k * a) ** 2))
            x = 2 * k * a * a / ((2 * k * a + 1) + disc)
            return x / a
        if self.conc_model == "self_complementary":
            # 2A <-> A2, toplam iplik C_T; x = [A2]
            s = self.c_eff_M  # C_T
            disc = math.sqrt(max(0.0, (4 * k * s + 1) ** 2 - 16 * (k * s) ** 2))
            x = 2 * k * s * s / ((4 * k * s + 1) + disc)
            return 2 * x / s
        raise ValueError(self.conc_model)

    def melting_curve(
        self, t_min: float = 20.0, t_max: float = 95.0, step: float = 0.25
    ) -> list[tuple[float, float]]:
        """[(T_C, theta)] listesi (TH-05)."""
        out = []
        n = int(round((t_max - t_min) / step))
        for i in range(n + 1):
            t = t_min + i * step
            out.append((round(t, 4), self.theta(t)))
        return out

    def as_dict(self) -> dict:
        return {
            "conc_model": self.conc_model,
            "salt_model": self.salt_model,
            "c_eff_M": self.c_eff_M,
            "tm_1M_c": round(self.tm_1m_K - KELVIN, 4),
            "tm_c": round(self.tm_c, 4),
            "dH_kcal_mol": round(self.duplex.dH_kcal, 4),
            "dS_1M_cal_mol_K": round(self.duplex.dS_cal, 4),
            "dS_eff_cal_mol_K": round(self.dS_eff_cal, 4),
            "ds_rescaled_to_salt_corrected_tm": True,
            "dG37_kcal_mol": round(self.dG_kcal(37.0), 4),
            "salt_regime": self.salt_detail.get("regime"),
        }


def build(
    seq: str,
    cond: Conditions | None = None,
    conc_model: str | None = None,
    salt_model: str = "owczarzy_2008",
    total_strand_M: float | None = None,
    *,
    template_site: str | None = None,
    arms_mismatch: bool = False,
) -> MeltingModel:
    """Bir oligo icin erime modelini kurar."""
    cond = cond or Conditions()
    if arms_mismatch:
        if template_site is None:
            raise ValueError("CT-04 uyumsuz alel dizisi gerekli")
        duplex = duplex_with_arms_mismatch(seq, template_site)
    else:
        duplex = (watson_crick(seq) if template_site is None
                  else duplex_with_mismatches(seq, template_site))
    f_gc = gc_fraction(seq)
    conc_model = conc_model or default_model(duplex)
    if total_strand_M is None:
        total_strand_M = cond.primer_M

    c_eff = effective_concentration(conc_model, cond, total_strand_M)
    if c_eff <= 0:
        raise ValueError("etkin konsantrasyon pozitif olmali")

    denom = duplex.dS_cal + R_CAL * math.log(c_eff)
    if denom >= 0:
        raise ValueError("gecersiz termodinamik: dS + R ln C >= 0")
    tm_1m_K = duplex.dH_kcal * 1000.0 / denom

    tm_K, detail = correct(tm_1m_K, f_gc, duplex.n_bp, cond, salt_model)
    dS_eff = duplex.dH_kcal * 1000.0 / tm_K - R_CAL * math.log(c_eff)

    return MeltingModel(
        duplex=duplex,
        cond=cond,
        f_gc=f_gc,
        conc_model=conc_model,
        salt_model=salt_model,
        c_eff_M=c_eff,
        tm_1m_K=tm_1m_K,
        tm_K=tm_K,
        dS_eff_cal=dS_eff,
        salt_detail=detail,
    )
