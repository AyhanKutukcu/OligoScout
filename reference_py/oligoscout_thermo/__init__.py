"""OligoScout termodinamik cekirdegi — Python referans uygulamasi.

Yalnizca standart kutuphane kullanir (bagimlilik yok). C++20 surumune
tasindiginda ikisi arasinda diferansiyel test calistirilacaktir.

Moduller:
  params         parametre tablolari, kunye ve SHA-256
  nn             Watson-Crick nearest-neighbor duplex termodinamigi
  salt           Owczarzy 2004 / 2008, von Ahsen 2001 tuz duzeltmeleri
  concentration  uc konsantrasyon modeli, Tm, theta(T), erime egrisi
  engine         QC kontrolleri ve REP-02 kanit kayitlari
  align          dimer ve sac tokasi (tek surekli WC sarmali)
  fold           Zuker MFE ve McCaskill bolusum fonksiyonu
  design         primer tasarimi, Pareto siralama
  multiplex      MPX-01..MPX-03 capraz reaktivite
  amplicon       QPCR-02 amplikon erime profili
  probe          PRB-01..PRB-05 hidroliz probu (TaqMan/ddPCR)
  variants       VAR-01..VAR-08 varyant/alel dusmesi katmani
  profiles       uygulama profilleri (qPCR, ctDNA, Sanger, ...)
"""

from .amplicon import MeltProfile, melt_profile, qpcr02
from .align import DimerResult, best_dimer, self_dimer
from .concentration import MeltingModel, build
from .fold import (
    DeclaredLoopEnergies,
    PendingLoopEnergies,
    ensemble_free_energy,
    mccaskill,
    unpaired_probability,
    zuker_mfe,
)
from .engine import ENGINE_VERSION, Profile, ThermoEngine
from .multiplex import (
    BindingSite,
    CrossProduct,
    assemble_cross_products,
    dimer_matrix,
    iter_products,
    load_bindings,
    size_separation,
    tm_homogeneity,
)
from .variants import (
    Site,
    Variant,
    VariantConstraints,
    evaluate_amplicon,
    evaluate_site,
)
from .probe import (
    ProbeCandidate,
    ProbeConstraints,
    enumerate_probes,
    evaluate_probe,
    pick_probe,
)
from .profiles import (
    PROFILES,
    ApplicationProfile,
    detectability,
    genome_equivalents,
    intact_fragment_probability,
)
from .nn import DuplexThermo, gc_fraction, reverse_complement, watson_crick
from .params import ParameterTableUnavailable, registry_digest, registry_summary
from .salt import Conditions

__all__ = [
    "ApplicationProfile",
    "BindingSite",
    "Conditions",
    "CrossProduct",
    "PROFILES",
    "DuplexThermo",
    "ENGINE_VERSION",
    "MeltProfile",
    "MeltingModel",
    "ParameterTableUnavailable",
    "Profile",
    "ProbeCandidate",
    "ProbeConstraints",
    "Site",
    "Variant",
    "VariantConstraints",
    "ThermoEngine",
    "assemble_cross_products",
    "build",
    "detectability",
    "dimer_matrix",
    "gc_fraction",
    "genome_equivalents",
    "intact_fragment_probability",
    "iter_products",
    "enumerate_probes",
    "evaluate_amplicon",
    "evaluate_probe",
    "evaluate_site",
    "load_bindings",
    "melt_profile",
    "qpcr02",
    "pick_probe",
    "size_separation",
    "tm_homogeneity",
    "registry_digest",
    "registry_summary",
    "reverse_complement",
    "watson_crick",
]
