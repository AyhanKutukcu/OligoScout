"""Bisulfit-donusturulmus DNA primer kontrolleri (BIS-01..BIS-04).

Konvansiyon: primerler bisulfit-DONUSTURULMUS hedef iplik uzerinde tasarlanir
(MethPrimer vb. ile ayni). Bisulfit muamelesi metillenmemis sitozini urasile
(dizilemede T) cevirir; CpG disi tum C'ler T olur, CpG'deki C metilasyon
durumuna gore degisken kalir.

  BIS-01  primer / CpG icerigi        (primer dizisinden; MSP disi CpG istenmez)
  BIS-02  primer / Donusum seciciligi (ORIJINAL hedef footprint'inde CpG disi C)
  BIS-03  primer / Donusmus karmasiklik (CMP-05 entropisi donusmus primerde)
  BIS-04  cift   / Iplik ozgullugu    (CT/GA donusmus genom taramasindan)

BIS-01/02/03 saf dizi kontrolleridir; BIS-04 donusmus genom FM-index'ini (CT/GA)
ve bir OligoScout taramasini gerektirir. Esikler literatur baslangic degerleridir
ve KALIBRE DEGILDIR (REP-03).
"""

from __future__ import annotations

from . import nn, seqpair


def _rec(cid, value, threshold, severity, note=""):
    return {"id": cid, "value": value, "threshold": threshold,
            "severity": severity, "note": note}


def count_cpg(s: str) -> int:
    """Dizideki CpG (5'-CG-3') dinukleotid sayisi."""
    return sum(1 for i in range(len(s) - 1) if s[i] == "C" and s[i + 1] == "G")


def nearest_cpg_distance_from_3p(s: str):
    """3' uca en yakin CpG'nin C'sinin 3' uctan uzakligi (nt). CpG yoksa None."""
    best = None
    for i in range(len(s) - 1):
        if s[i] == "C" and s[i + 1] == "G":
            d = len(s) - 1 - i          # C konumunun 3' uca uzakligi
            best = d if best is None else min(best, d)
    return best


def count_non_cpg_c(s: str) -> int:
    """CpG baglaminda OLMAYAN C sayisi (bisulfit'te T'ye donusen C'ler)."""
    n = 0
    for i, ch in enumerate(s):
        if ch == "C" and (i + 1 >= len(s) or s[i + 1] != "G"):
            n += 1
    return n


def check_bis01(primer: str, msp: bool = False, warn_3p_nt: int = 5) -> dict:
    """BIS-01: primerde CpG. MSP (metilasyona-ozgu PCR) disinda CpG istenmez;
    bir CpG 3' uca yakinsa risk daha yuksektir (metilasyon durumu ayrimi)."""
    s = nn.validate(primer)
    n = count_cpg(s)
    d3 = nearest_cpg_distance_from_3p(s)
    if msp:
        sev = "INFO"
    elif n == 0:
        sev = "INFO"
    elif d3 is not None and d3 < warn_3p_nt:
        sev = "FAIL"
    else:
        sev = "WARN"
    return _rec("BIS-01",
                {"cpg_count": n, "nearest_cpg_dist_from_3p": d3, "msp": msp},
                {"warn": "CpG var (MSP disi)",
                 "fail": f"CpG 3' son {warn_3p_nt} nt icinde"},
                sev,
                "MSP disi primerde CpG metilasyon durumuna gore ayrimli baglanir")


def check_bis02(original_footprint: str, warn_min: int = 2) -> dict:
    """BIS-02: donusum seciciligi. Primerin ORIJINAL (donusum oncesi) hedef
    footprint'indeki CpG disi C sayisi. Bu C'ler donusmus sablonda T olur;
    primer bu ayrimli konumlari kapsadigi olcude DONUSMUS DNA'ya ozgu olur
    (donusmemis DNA'da hala C vardir). Az sayida ayrimli konum -> dusuk
    donusum seciciligi (donusmemis DNA'yi da cogaltabilir)."""
    s = original_footprint.upper()
    n = count_non_cpg_c(s)
    sev = "WARN" if n < warn_min else "INFO"
    return _rec("BIS-02",
                {"non_cpg_c_count": n, "footprint_len": len(s)},
                {"warn": f"< {warn_min} CpG-disi C (dusuk donusum seciciligi)"},
                sev,
                "donusmus DNA'ya ozguluk icin primer donusen (C->T) konum icermeli")


def check_bis03(converted_primer: str, warn_bits: float = 1.2,
                fail_bits: float = 0.9) -> dict:
    """BIS-03: donusmus dizide karmasiklik (CMP-05 entropisi). Bisulfit sonrasi
    alfabe pratik olarak 3 harfe iner (C seyrekler), bu yuzden karmasiklik
    dogal olarak dusuktur; bu kontrol asiri dusuk karmasikligi (yanlis onceleme
    / off-target riski) yakalar. Donusmus primer uzerinde hesaplanir."""
    s = nn.validate(converted_primer)
    h1 = seqpair.shannon_entropy(s, 1)
    h2 = seqpair.shannon_entropy(s, 2)
    if h1 < fail_bits:
        sev = "FAIL"
    elif h1 < warn_bits:
        sev = "WARN"
    else:
        sev = "INFO"
    return _rec("BIS-03",
                {"shannon_k1_bit": round(h1, 4), "shannon_k2_bit": round(h2, 4)},
                {"warn": f"k1 < {warn_bits} bit", "fail": f"k1 < {fail_bits} bit"},
                sev,
                "donusmus dizi 3-harfli; esik bisulfit icin gevsetilmistir")


def check_bis04(ct_intended: int, ct_total: int,
                ga_intended: int, ga_total: int,
                intended_conversion: str = "CT") -> dict:
    """BIS-04: iplik ozgullugu. Cift, DONUSTURULMUS genomun dort ipliginde
    (CT ve GA indeksleri, her biri iki iplik) arandiktan sonra: amaclanan
    donusum ipliginde 1 hedef urun beklenir; diger donusum ipliginde urun
    iplik-capraz baglanmadir (risk).

    ct_* / ga_*: CT ve GA indekslerinde bulunan (amaclanan, toplam) urun
    sayilari. intended_conversion: "CT" ust iplik, "GA" alt iplik tasarimi.
    """
    if intended_conversion.upper() == "CT":
        intended_total, off_total = ct_total, ga_total
    else:
        intended_total, off_total = ga_total, ct_total
    off_extra = off_total                      # diger donusumde herhangi urun
    intended_extra = max(0, intended_total - 1)  # amaclanan disi fazla urun
    sev = "INFO"
    if off_extra > 0:
        sev = "FAIL"            # yanlis donusum ipliginde urun
    elif intended_extra > 0:
        sev = "WARN"            # amaclanan iplikte fazladan urun
    return _rec("BIS-04",
                {"intended_conversion": intended_conversion.upper(),
                 "ct_products_total": ct_total, "ga_products_total": ga_total,
                 "intended_products_total": intended_total,
                 "offstrand_products_total": off_total},
                {"warn": "amaclanan donusum ipliginde fazladan urun",
                 "fail": "yanlis donusum ipliginde urun (iplik-capraz)"},
                sev,
                "donusmus genom CT/GA taramasindan; donusmemis genom ayri")
