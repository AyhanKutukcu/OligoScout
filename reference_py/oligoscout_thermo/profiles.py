"""Uygulama profilleri: her deney tipinin kendi kisitlari ve esikleri.

QC katalogundaki uygulama bolumlerinin (qPCR, prob, mRNA/RT-qPCR, ctDNA,
bisulfit, Sanger, NGS amplikon, multipleks) makineye okunur hali. Bir profil
uc sey tasir:
  * tasarim kisitlari   (design.DesignConstraints)
  * QC esikleri         (engine.Profile)
  * profile ozgu kontroller ve bunlarin durumu

Her esik literatur kilavuz degeridir ve kalibre DEGILDIR; profil adi ve esik
degerleri her rapor satirinda beyan edilir.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .design import DesignConstraints
from .engine import Profile


@dataclass(frozen=True)
class ApplicationProfile:
    name: str
    description: str
    constraints: DesignConstraints
    qc: Profile
    notes: tuple[str, ...] = ()
    pending_checks: tuple[str, ...] = ()

    def as_dict(self) -> dict:
        return {
            "profile": self.name,
            "description": self.description,
            "product_bp": [self.constraints.product_min,
                           self.constraints.product_max],
            "tm_window_c": list(self.qc.tm_window_c),
            "length_nt": [self.constraints.length_min,
                          self.constraints.length_max],
            "gc_fraction": [self.constraints.gc_min, self.constraints.gc_max],
            "notes": list(self.notes),
            "pending_checks": list(self.pending_checks),
        }


# ---------------------------------------------------------------------------
# ctDNA tespit edilebilirligi (CT-02, CT-03)
# ---------------------------------------------------------------------------
def intact_fragment_probability(amplicon_bp: int,
                                fragment_bp: int = 166) -> float:
    """CT-02: L uzunlugundaki bir fragmanin A uzunlugundaki amplikonu
    butun olarak icerme olasiligi = max(0, L - A + 1) / L.

    cfDNA fragman uzunlugu ~166 bp'de tepe yapar (Snyder ve ark. 2016).
    Dagilim verilirse `detectability` uzerinden agirlikli ortalama alinir.
    """
    if amplicon_bp <= 0 or fragment_bp <= 0:
        raise ValueError("uzunluklar pozitif olmali")
    return max(0.0, fragment_bp - amplicon_bp + 1) / fragment_bp


def genome_equivalents(cfdna_ng: float) -> float:
    """1 haploid genom ~3.3 pg; GE = ng / 3.3e-3."""
    return cfdna_ng / 3.3e-3


def detectability(amplicon_bp: int, vaf: float, cfdna_ng: float = 10.0,
                  fragment_bp: int = 166,
                  fragment_distribution: dict[int, float] | None = None
                  ) -> dict[str, float]:
    """CT-03: en az bir mutant kopya yakalama olasiligi (Poisson).

    lambda = GE * VAF * P(butun fragman);  P(>=1) = 1 - exp(-lambda)
    """
    if fragment_distribution:
        total = sum(fragment_distribution.values())
        p_intact = sum(
            w * intact_fragment_probability(amplicon_bp, L)
            for L, w in fragment_distribution.items()
        ) / total
    else:
        p_intact = intact_fragment_probability(amplicon_bp, fragment_bp)

    ge = genome_equivalents(cfdna_ng)
    lam = ge * vaf * p_intact
    return {
        "amplicon_bp": amplicon_bp,
        "p_intact_fragment": round(p_intact, 6),
        "genome_equivalents": round(ge, 1),
        "vaf": vaf,
        "lambda": round(lam, 4),
        "p_detect": round(1.0 - math.exp(-lam), 6),
    }


# ---------------------------------------------------------------------------
# Profiller
# ---------------------------------------------------------------------------
def _profile(name, desc, *, product, tm, length=(18, 27), gc=(0.40, 0.60),
             notes=(), pending=()) -> ApplicationProfile:
    return ApplicationProfile(
        name=name, description=desc,
        constraints=DesignConstraints(
            length_min=length[0], length_max=length[1],
            tm_min=tm[0], tm_max=tm[1], tm_opt=(tm[0] + tm[1]) / 2,
            gc_min=gc[0], gc_max=gc[1],
            product_min=product[0], product_max=product[1],
        ),
        qc=Profile(name=name, tm_window_c=tm),
        notes=notes, pending_checks=pending,
    )


PROFILES: dict[str, ApplicationProfile] = {
    "qpcr_sybr": _profile(
        "qpcr_sybr", "SYBR tabanli qPCR",
        product=(70, 150), tm=(58.0, 62.0),
        notes=("QPCR-01: amplikon 70-150 bp (MIQE)",
               "QPCR-02: tek erime piki beklenir; erime modeli "
               "oligoscout_thermo.amplicon, CLI: cli/amplicon_melt.py",
               "QPCR-02 model sabiti (sigma) kalibre degil, dongu "
               "entropisi tablolari bekliyor"),
    ),
    "rt_qpcr": _profile(
        "rt_qpcr", "mRNA / RT-qPCR",
        product=(70, 150), tm=(58.0, 62.0),
        notes=("RT-01: ekson-ekson kavsagi tercih edilir",
               "RT-02: gDNA urunu ayrica taranmali"),
        pending=("RT-01/RT-03: transkriptom indeksi henuz kurulmadi",
                 "RT-04: RNA dongu parametreleri aktarilmadi"),
    ),
    "ctdna": _profile(
        "ctdna", "ctDNA / cfDNA, kisa amplikon",
        product=(60, 100), tm=(58.0, 62.0), length=(18, 25),
        notes=("CT-01: cfDNA ~166 bp tepe; amplikon <=100 bp",
               "CT-02/CT-03: butun fragman ve Poisson tespit olasiligi",
               "CT-04: ARMS alel baglanma ayrimi, KALIBRE DEGIL"),
    ),
    "ddpcr": _profile(
        "ddpcr", "ddPCR",
        product=(60, 120), tm=(58.0, 62.0),
        notes=("PRB-01: prob Tm'i primerlerden 6-10 C yuksek",
               "PRB-02: 5' uc G yasak; PRB-03: C sayisi > G sayisi",
               "PRB-04: prob primer baglanma yerleriyle cakisamaz",
               "PRB-05: prob self/x-primer dimerleri (tekli ic uyumsuzluk dahil)",
               "prob modulu: oligoscout_thermo.probe, CLI: cli/probe_qc.py"),
    ),
    "sanger": _profile(
        "sanger", "Sanger dizileme",
        product=(300, 800), tm=(55.0, 65.0), length=(18, 30),
        notes=("SAN-01: ilgi bolgesi primere >=60 nt uzak olmali",
               "SAN-02: amplikon 300-800 bp"),
    ),
    "ngs_amplicon": _profile(
        "ngs_amplicon", "NGS amplikon paneli",
        product=(120, 280), tm=(58.0, 62.0),
        notes=("NGS-01: adaptor kuyruklari dimer analizine dahil edilmeli",),
        pending=("NGS-01: kuyruk dahil analiz (DIM-07) CLI'ye baglanmadi",),
    ),
    "standard_pcr": _profile(
        "standard_pcr", "Standart PCR",
        product=(100, 1000), tm=(58.0, 66.0), length=(18, 30),
    ),
    "bisulfite": _profile(
        "bisulfite", "Bisulfit donusturulmus DNA",
        product=(100, 300), tm=(50.0, 60.0), length=(20, 30), gc=(0.20, 0.50),
        notes=("BIS-01: primerde CpG bulunmamali (MSP disi)",),
        pending=("BIS-01..BIS-04: donusturulmus genom indeksleri kurulmadi",),
    ),
}


def get(name: str) -> ApplicationProfile:
    if name not in PROFILES:
        raise KeyError(f"bilinmeyen profil: {name}. "
                       f"Secenekler: {', '.join(sorted(PROFILES))}")
    return PROFILES[name]
