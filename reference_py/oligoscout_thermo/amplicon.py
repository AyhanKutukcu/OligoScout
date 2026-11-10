"""Amplikon erime profili (QPCR-02): Poland tipi kismi erime modeli.

Bir PCR urunu tek bir "iki durumlu" molekul gibi erimez: AT yogun bolgeler
once acilir, GC yogun bolgeler sonra. SYBR gibi cift iplik boyasiyla okunan
bir qPCR'de bu, erime egrisinde BIRDEN COK PIK demektir ve QPCR-02 tam olarak
bunu sorar: urun tek pik veriyor mu?

Model (Zimm & Bragg 1959 tipi isbirlikli helis-kangal; Blake & Delcourt 1998):
  * Her baz cifti ya helis icinde ya da kangalda.
  * Helis yigin agirligi s_i = exp(-dG_i(T)/RT), dG_i WC NN tablosundan.
  * Yeni helis baslatmak sigma ile cezalandirilir (isbirligi).
  * Uclar serbest (fraying); transfer matrisi uzerinden her baz ciftinin
    eslesme olasiligi ve theta(T) hesaplanir.

ACIKCA BEYAN EDILEN SINIRLAR
  1. sigma ve alpha LITERATUR MODEL SABITLERIDIR, bu veri uzerinde
     KALIBRE EDILMEMISTIR (sigma=1e-5, alpha=1.75).
  2. Poland dongu ENTROPISI (xi(n) = n^-alpha) bu surumde YOKTUR; dongu
     parametre tablolari `pending` oldugu icin isbirligi tek bir sigma ile
     temsil edilir. Bunun etkisi gecislerin gercekte olduklarindan daha
     yumusak cikmasidir; pik SAYISI ve siralamasi korunur, pik genisligi
     oldugundan buyuk tahmin edilir. Tablolar girildiginde Poland dongu
     terimi eklenecektir.
  3. Tuz duzeltmesi baz cifti basina entropi kaymasidir
     (SantaLucia 1998: dS += 0.368 * ln[Na+] cal/(mol K) / baz cifti);
     Owczarzy duzeltmeleri oligo Tm'i icindir, buraya uygulanmaz.
  4. Boya baglanmasi (SYBR) modellenmez; olculen egri ile tahmin arasinda
     sistematik kayma beklenir. Karar pik SAYISI ve pikler arasi FARK
     uzerinden verilir, mutlak Tm uzerinden degil.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from . import nn, salt
from .params import KELVIN, R_CAL, WC_SANTALUCIA_1998

SIGMA = 1.0e-5        # dongu acilisi isbirligi sabiti (model sabiti)
ALPHA = 1.75          # Poland dongu usu — tablolar gelince kullanilacak


@dataclass(frozen=True)
class MeltProfile:
    temperatures: list[float]
    helicity: list[float]          # theta(T): eslesmis baz cifti orani
    derivative: list[float]        # -dtheta/dT
    peaks: list[tuple[float, float]]   # (T, -dtheta/dT) yerel tepe noktalari
    tm_c: float                    # theta = 0.5 gectigi sicaklik
    n_bp: int

    def as_dict(self) -> dict:
        return {
            "tm_c": round(self.tm_c, 3),
            "n_bp": self.n_bp,
            "n_peaks": len(self.peaks),
            "peaks_c": [round(t, 2) for t, _ in self.peaks],
            "peak_heights": [round(h, 5) for _, h in self.peaks],
            "model": "zimm_bragg_cooperative",
            "model_constants": {"sigma": SIGMA, "calibrated": False},
            "loop_entropy": "yok (dongu tablolari pending)",
        }


def _stack_energies(seq: str, cond: salt.Conditions) -> list[tuple[float, float]]:
    """Her komsu baz cifti basamagi icin (dH kcal/mol, dS cal/mol/K).

    Tuz duzeltmesi basamak basina entropiye eklenir (SantaLucia 1998).
    """
    table = WC_SANTALUCIA_1998.require()["nn"]
    na = salt.equivalent_sodium(cond) if hasattr(salt, "equivalent_sodium") \
        else cond.monovalent_M
    na = max(na, 1e-6)
    ds_salt = 0.368 * math.log(na)          # cal/(mol K) / baz cifti
    out = []
    for i in range(len(seq) - 1):
        step = seq[i:i + 2]
        dH, dS = table.get(step) or table[nn.reverse_complement(step)]
        out.append((dH, dS + ds_salt))
    return out


def _theta_at(stacks: list[tuple[float, float]], t_c: float) -> float:
    """Eslesmis baz cifti orani theta(T), transfer matrisiyle (Zimm-Bragg).

    Durumlar: her baz cifti helis (h) ya da kangal (c). Basamak agirligi
    s_i = exp(-dG_i(T)/RT); yeni bir helis baslatmak sigma ile cezalandirilir.
    Ileri/geri vektorlerin carpimi her baz ciftinin eslesme olasiligini verir;
    her adimda normalize edildigi icin uzun amplikonlarda tasma olmaz.
    Olcek carpanlari oranda sadelesir.
    """
    T = t_c + KELVIN
    RT = R_CAL * T / 1000.0                       # kcal/mol
    s = [math.exp(-(dH - T * dS / 1000.0) / RT) for dH, dS in stacks]
    n = len(s) + 1

    fwd = [[SIGMA, 1.0]]
    for si in s:
        h, c = fwd[-1]
        nh = h * si + c * SIGMA * si
        nc = h + c
        m = max(nh, nc, 1e-300)
        fwd.append([nh / m, nc / m])

    bwd = [[1.0, 1.0] for _ in range(n)]
    for i in range(n - 2, -1, -1):
        bh, bc = bwd[i + 1]
        si = s[i]
        nh = si * bh + bc
        nc = SIGMA * si * bh + bc
        m = max(nh, nc, 1e-300)
        bwd[i] = [nh / m, nc / m]

    total = 0.0
    for i in range(n):
        ph = fwd[i][0] * bwd[i][0]
        pc = fwd[i][1] * bwd[i][1]
        total += ph / (ph + pc) if (ph + pc) > 0 else 0.0
    return total / n


def melt_profile(sequence: str, cond: salt.Conditions | None = None,
                 t_min: float = 50.0, t_max: float = 98.0,
                 step: float = 0.25) -> MeltProfile:
    """Amplikonun erime profili, turevi ve pikleri."""
    cond = cond or salt.Conditions()
    seq = nn.validate(sequence)
    if len(seq) < 2:
        raise ValueError("amplikon en az 2 nt olmali")
    stacks = _stack_energies(seq, cond)

    temps, theta = [], []
    t = t_min
    while t <= t_max + 1e-9:
        temps.append(round(t, 4))
        theta.append(_theta_at(stacks, t))
        t += step

    deriv = [0.0] * len(temps)
    for i in range(1, len(temps) - 1):
        deriv[i] = -(theta[i + 1] - theta[i - 1]) / (temps[i + 1] - temps[i - 1])

    peaks = []
    for i in range(2, len(temps) - 2):
        if deriv[i] >= deriv[i - 1] and deriv[i] > deriv[i + 1] and deriv[i] > 0:
            peaks.append((temps[i], deriv[i]))
    if peaks:
        tallest = max(h for _, h in peaks)
        peaks = [(t_, h) for t_, h in peaks if h >= 0.15 * tallest]

    tm = float("nan")
    for i in range(1, len(temps)):
        if theta[i - 1] >= 0.5 > theta[i]:
            f = (theta[i - 1] - 0.5) / (theta[i - 1] - theta[i])
            tm = temps[i - 1] + f * (temps[i] - temps[i - 1])
            break

    return MeltProfile(temperatures=temps, helicity=theta, derivative=deriv,
                       peaks=peaks, tm_c=tm, n_bp=len(seq))


def qpcr02(profile: MeltProfile, min_separation_c: float = 2.0) -> dict:
    """QPCR-02: tek erime piki beklenir.

    Birden cok pik, urunun kismi erime alan adanlari oldugunu gosterir;
    SYBR okumasinda bu, tek urun olsa bile coklu pik demektir.
    """
    peaks = profile.peaks
    sev = "INFO"
    note = "tek pik"
    if len(peaks) > 1:
        gaps = [b[0] - a[0] for a, b in zip(peaks, peaks[1:])]
        sep = max(gaps) if gaps else 0.0
        sev = "WARN" if sep >= min_separation_c else "INFO"
        note = (f"{len(peaks)} pik, en buyuk ayrim {sep:.2f} C"
                if sev == "WARN" else
                f"{len(peaks)} pik, hepsi {min_separation_c} C icinde")
    return {"check": "QPCR-02", "severity": sev, "n_peaks": len(peaks),
            "peaks_c": [round(t, 2) for t, _ in peaks],
            "tm_c": round(profile.tm_c, 3), "note": note,
            "threshold": f"warn: >1 pik ve ayrim >= {min_separation_c} C"}
