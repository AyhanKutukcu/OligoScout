"""Dizi, kompozisyon, cift uyumu ve amplikon kontrolleri (kendine yeten set).

QC kataloğundan, DIS VERI ya da yeni parametre tablosu GEREKTIRMEYEN
kurallar. Hepsi mevcut Watson-Crick tablosu ve dizi analizinden hesaplanir:

  SEQ-02  primer uzunlugu
  SEQ-03  dejenerasyon katsayisi (IUPAC)
  CMP-05  dusuk karmasiklik (Shannon entropisi k=1,2)
  CMP-06  G-quadruplex motifi
  CMP-07  5'/3' kararlilik dengesi
  CMP-08  3' terminal baz
  SS-02   3' uclu hairpin
  PAIR-02 urun boyu (profil penceresi)
  PAIR-03 primer cakismasi / arasi mesafe
  PAIR-04 urun Tm'i - primer Tm'i
  AMP-01  amplikon GC (global + 50 nt pencere)
  AMP-03  amplikonda tekrar / dusuk karmasiklik (CMP-05/06 amplikona)

Kapsam disi (baglam ya da veri gerektirir, burada DEGIL):
  SEQ-04/05 (kalip N / koordinat tutarliligi; kullanici koordinati gerekir),
  SPC-* (ozgulluk; genom taramasi + anotasyon), CT-05/06 (blocker / CHIP gen
  listesi), DIM-05/07 (prob / kuyruk baglami).

ESIKLER QC kataloğunun literatur baslangic degerleridir ve KALIBRE
DEGILDIR (REP-03). Esikler `SeqPairConstraints` alanlarindan degistirilir;
kural kimlikleri ve kayit semasi sabit kalir.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

from . import align, concentration, nn, oligomod, salt

# IUPAC dejenere kodlarin genislemesi (SEQ-03).
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


@dataclass(frozen=True)
class SeqPairConstraints:
    """Profil tarafindan ezilebilen esikler (REP-03, kalibre degil)."""

    length_warn: tuple[int, int] = (18, 30)
    length_fail: tuple[int, int] = (15, 35)
    degeneracy_warn: int = 32
    degeneracy_fail: int = 256
    entropy_warn: float = 1.5       # Shannon (k=1), bit; altinda WARN
    entropy_fail: float = 1.2
    product_bp: tuple[int, int] = (60, 200)
    product_edge_frac: float = 0.10
    pair_gap_min_nt: int = 10
    product_primer_tm_gap_warn: float = 30.0
    amplicon_gc_warn: tuple[float, float] = (0.35, 0.65)
    amplicon_gc_fail_window: float = 0.85
    gc_window_nt: int = 50
    ss_3p_warn_dg37: float = -2.0


def _rec(cid, value, threshold, severity, note=""):
    return {"id": cid, "value": value, "threshold": threshold,
            "severity": severity, "note": note}


# ---------------------------------------------------------------------------
# Dizi yardimcilari
# ---------------------------------------------------------------------------
def degeneracy(seq: str) -> int:
    """SEQ-03: her IUPAC konumundaki baz sayisinin carpimi."""
    product = 1
    for base in seq.strip().upper():
        if base not in IUPAC:
            raise ValueError(f"tanimsiz baz: {base}")
        product *= len(IUPAC[base])
    return product


def three_prime_has_degenerate(seq: str, window: int = 3) -> bool:
    return any(len(IUPAC.get(b, "")) > 1 for b in seq.strip().upper()[-window:])


def shannon_entropy(seq: str, k: int = 1) -> float:
    """k-mer Shannon entropisi (bit). Dusuk = tekrarli / dusuk karmasiklik."""
    s = seq.strip().upper()
    if len(s) < k:
        return 0.0
    counts: dict[str, int] = {}
    total = 0
    for i in range(len(s) - k + 1):
        counts[s[i:i + k]] = counts.get(s[i:i + k], 0) + 1
        total += 1
    return -sum((n / total) * math.log2(n / total) for n in counts.values())


_G4 = None


def has_g_quadruplex(seq: str) -> bool:
    """CMP-06: G>=3 N1-7 dort kez, her iki iplikte."""
    global _G4
    if _G4 is None:
        import re
        _G4 = re.compile(r"G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}")
    s = seq.strip().upper()
    return bool(_G4.search(s) or _G4.search(nn.reverse_complement(s)))


# ---------------------------------------------------------------------------
# Tek primer kontrolleri
# ---------------------------------------------------------------------------
def check_primer(seq: str, cons: SeqPairConstraints | None = None) -> list[dict]:
    cons = cons or SeqPairConstraints()
    raw = seq.strip().upper()
    out: list[dict] = []

    # SEQ-02 uzunluk
    n = len(raw)
    lw, hw = cons.length_warn
    lf, hf = cons.length_fail
    sev = "INFO"
    if n < lf or n > hf:
        sev = "FAIL"
    elif n < lw or n > hw:
        sev = "WARN"
    out.append(_rec("SEQ-02", {"length_nt": n},
                    {"warn": f"{lw}-{hw} disi", "fail": f"{lf}-{hf} disi"}, sev))

    # SEQ-03 dejenerasyon
    deg = degeneracy(raw)
    deg3 = three_prime_has_degenerate(raw)
    sev = "INFO"
    if deg > cons.degeneracy_fail or deg3:
        sev = "FAIL"
    elif deg > cons.degeneracy_warn:
        sev = "WARN"
    out.append(_rec("SEQ-03",
                    {"degeneracy": deg, "degenerate_in_3p_last3": deg3},
                    {"warn": f">{cons.degeneracy_warn}",
                     "fail": f">{cons.degeneracy_fail} veya 3' son 3 nt dejenere"},
                    sev,
                    "" if deg > 1 else "dejenere baz yok"))

    # Kalan kontroller yalniz saf A/C/G/T diziler icin tanimli.
    if deg > 1:
        return out
    s = nn.validate(raw)

    # CMP-05 dusuk karmasiklik
    h1 = shannon_entropy(s, 1)
    h2 = shannon_entropy(s, 2)
    sev = "INFO"
    if h1 < cons.entropy_fail:
        sev = "FAIL"
    elif h1 < cons.entropy_warn:
        sev = "WARN"
    out.append(_rec("CMP-05", {"shannon_k1_bit": round(h1, 4),
                               "shannon_k2_bit": round(h2, 4)},
                    {"warn": f"k1 < {cons.entropy_warn} bit",
                     "fail": f"k1 < {cons.entropy_fail} bit"}, sev))

    # CMP-06 G-quadruplex
    g4 = has_g_quadruplex(s)
    out.append(_rec("CMP-06", {"g_quadruplex": g4},
                    {"fail": "G>=3 N1-7 dort tekrar (primer)"},
                    "FAIL" if g4 else "INFO"))

    # CMP-07 5'/3' kararlilik dengesi
    half = len(s) // 2
    dg5 = nn.watson_crick(s[:half]).dG_kcal() if half >= 2 else 0.0
    dg3 = nn.watson_crick(s[-half:]).dG_kcal() if half >= 2 else 0.0
    # 3' yari belirgin daha kararli (daha negatif) ise WARN.
    sev = "WARN" if dg3 < dg5 - 2.0 else "INFO"
    out.append(_rec("CMP-07",
                    {"dG5_half_kcal_mol": round(dg5, 3),
                     "dG3_half_kcal_mol": round(dg3, 3),
                     "delta_kcal_mol": round(dg3 - dg5, 3)},
                    {"warn": "3' yari > 2 kcal/mol daha kararli"}, sev,
                    "3' ucun asiri kararli olmasi yanlis onceleme riski"))

    # CMP-08 3' terminal baz
    last = s[-1]
    out.append(_rec("CMP-08", {"three_prime_base": last},
                    {"warn": "3'-T (uyumsuzluga toleransli, Kwok 1990)"},
                    "WARN" if last == "T" else "INFO"))

    # SS-02 3' uclu hairpin
    hp = align.hairpin_stem(s)
    extendable = hp.found and hp.three_prime_paired
    dg = hp.dG37_kcal if hp.found else 0.0
    sev = "WARN" if (extendable and dg < cons.ss_3p_warn_dg37) else "INFO"
    out.append(_rec("SS-02",
                    {"stem_dG37_kcal_mol": round(dg, 3) if hp.found else None,
                     "three_prime_extendable": extendable},
                    {"warn": f"3' uzatilabilir ve dG37 < {cons.ss_3p_warn_dg37}"},
                    sev))
    return out


# ---------------------------------------------------------------------------
# Cift ve amplikon kontrolleri
# ---------------------------------------------------------------------------
def check_pair(fwd: str, rev: str, product_start: int, product_end: int,
               amplicon: str, cond: salt.Conditions | None = None,
               salt_model: str = "owczarzy_2008",
               cons: SeqPairConstraints | None = None) -> list[dict]:
    cons = cons or SeqPairConstraints()
    cond = cond or salt.Conditions()
    out: list[dict] = []
    product_len = product_end - product_start

    # PAIR-02 urun boyu
    lo, hi = cons.product_bp
    span = hi - lo
    edge = cons.product_edge_frac * span
    sev = "INFO"
    if product_len < lo or product_len > hi:
        sev = "FAIL"
    elif product_len < lo + edge or product_len > hi - edge:
        sev = "WARN"
    out.append(_rec("PAIR-02", {"product_bp": product_len},
                    {"warn": f"pencere kenarina < %{int(cons.product_edge_frac*100)}",
                     "fail": f"{lo}-{hi} bp disi"}, sev))

    # PAIR-03 primer cakismasi / arasi mesafe
    inner = product_len - len(fwd) - len(rev)
    sev = "INFO"
    if inner < 0:
        sev = "FAIL"
    elif inner < cons.pair_gap_min_nt:
        sev = "WARN"
    out.append(_rec("PAIR-03", {"inner_gap_nt": inner},
                    {"warn": f"ara < {cons.pair_gap_min_nt} nt", "fail": "ortusme"},
                    sev))

    # PAIR-04 urun Tm'i - primer Tm'i
    prod_tm = concentration.build(amplicon, cond=cond,
                                  salt_model=salt_model).tm_c
    tm_f = concentration.build(fwd, cond=cond, salt_model=salt_model).tm_c
    tm_r = concentration.build(rev, cond=cond, salt_model=salt_model).tm_c
    gap = prod_tm - max(tm_f, tm_r)
    out.append(_rec("PAIR-04",
                    {"product_tm_c": round(prod_tm, 3),
                     "max_primer_tm_c": round(max(tm_f, tm_r), 3),
                     "gap_c": round(gap, 3)},
                    {"warn": f"> {cons.product_primer_tm_gap_warn} C"},
                    "WARN" if gap > cons.product_primer_tm_gap_warn else "INFO"))

    # AMP-01 amplikon GC (global + 50 nt pencere maksimumu)
    gc_global = nn.gc_fraction(amplicon)
    w = cons.gc_window_nt
    win_max = gc_global
    if len(amplicon) >= w:
        win_max = max(nn.gc_fraction(amplicon[i:i + w])
                      for i in range(len(amplicon) - w + 1))
    gw_lo, gw_hi = cons.amplicon_gc_warn
    sev = "INFO"
    if win_max > cons.amplicon_gc_fail_window:
        sev = "FAIL"
    elif gc_global < gw_lo or gc_global > gw_hi or win_max > 0.75:
        sev = "WARN"
    out.append(_rec("AMP-01",
                    {"gc_global": round(gc_global, 4),
                     "gc_window_max": round(win_max, 4), "window_nt": w},
                    {"warn": f"global %{int(gw_lo*100)}-%{int(gw_hi*100)} disi "
                             f"veya pencere > %75",
                     "fail": f"pencere > %{int(cons.amplicon_gc_fail_window*100)}"},
                    sev))

    # AMP-03 amplikonda tekrar / dusuk karmasiklik
    h1 = shannon_entropy(amplicon, 1)
    g4 = has_g_quadruplex(amplicon)
    sev = "WARN" if (h1 < cons.entropy_warn or g4) else "INFO"
    out.append(_rec("AMP-03",
                    {"amplicon_shannon_k1_bit": round(h1, 4),
                     "amplicon_g_quadruplex": g4},
                    {"warn": "dusuk karmasiklik ya da G-quadruplex (amplikonda)"},
                    sev))
    return out


# ---------------------------------------------------------------------------
# Kalip baglami gerektiren kontroller (SEQ-04, SEQ-05)
# ---------------------------------------------------------------------------
def _non_acgt_count(segment: str) -> tuple[int, int]:
    """(N sayisi, soft-mask/kucuk-harf sayisi). Girdi HAM (dogrulanmamis)."""
    n = sum(1 for b in segment if b in "Nn")
    soft = sum(1 for b in segment if b.islower() and b.upper() in "ACGT")
    return n, soft


def check_template(fwd: str, rev: str, product_start: int, product_end: int,
                   raw_target: str) -> list[dict]:
    """SEQ-04 (kalipta belirsiz baz) ve SEQ-05 (koordinat-dizi tutarliligi).

    `raw_target` HAM referans dizisidir (N ve soft-mask KORUNMUS). Primerler
    5'->3' verilir; ileri primer product_start'ta, geri primer product_end'de
    ters-tamamlayici olarak beklenir.
    """
    out: list[dict] = []
    amp_raw = raw_target[product_start:product_end]
    fwd_u, rev_u = fwd.strip().upper(), rev.strip().upper()
    fwd_region = raw_target[product_start:product_start + len(fwd_u)]
    rev_region = raw_target[product_end - len(rev_u):product_end]

    # SEQ-04 kalipta belirsiz baz
    amp_n, amp_soft = _non_acgt_count(amp_raw)
    prim_n = (_non_acgt_count(fwd_region)[0] + _non_acgt_count(rev_region)[0])
    prim_soft = (_non_acgt_count(fwd_region)[1] + _non_acgt_count(rev_region)[1])
    sev = "INFO"
    if prim_n > 0:
        sev = "FAIL"                       # primer baglanma bolgesinde N
    elif amp_n > 0 or amp_soft > 0 or prim_soft > 0:
        sev = "WARN"                       # amplikonda N / soft-mask
    out.append(_rec("SEQ-04",
                    {"amplicon_N": amp_n, "amplicon_softmask": amp_soft,
                     "primer_region_N": prim_n, "primer_region_softmask": prim_soft},
                    {"warn": "amplikonda N / soft-mask > 0",
                     "fail": "primer baglanma bolgesinde N"}, sev))

    # SEQ-05 koordinat-dizi tutarliligi
    ref_fwd = fwd_region.upper()
    ref_rev = reverse_complement_raw(rev_region.upper())
    fwd_ok = ref_fwd == fwd_u
    rev_ok = ref_rev == rev_u
    sev = "INFO" if (fwd_ok and rev_ok) else "FAIL"
    out.append(_rec("SEQ-05",
                    {"forward_matches_reference": fwd_ok,
                     "reverse_matches_reference": rev_ok,
                     "product_start": product_start, "product_end": product_end},
                    {"fail": "koordinattaki referans dizisi primer dizisiyle "
                             "uyusmuyor"}, sev,
                    "" if (fwd_ok and rev_ok) else
                    "tasarim koordinati ile dizi tutarsiz"))
    return out


def reverse_complement_raw(seq: str) -> str:
    """N'yi koruyan ters tamamlayici (SEQ-05 referans karsilastirmasi icin)."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


# ---------------------------------------------------------------------------
# CHIP geni uyarisi (CT-06)
# ---------------------------------------------------------------------------
# Klonal hematopoez (CHIP) ile iliskili yaygin surucu genler.
# Kaynak: Jaiswal 2014 (NEJM 371:2488), Razavi 2019 (Nat Med 25:1928).
# ctDNA panelinde bu genleri hedeflemek, tumor DNA'si yerine kan
# hucrelerinden gelen varyantlari yakalama riski tasir.
CHIP_GENES = frozenset({
    "DNMT3A", "TET2", "ASXL1", "JAK2", "TP53", "PPM1D", "SF3B1", "SRSF2",
    "GNB1", "CBL", "GNAS", "IDH1", "IDH2", "KRAS", "NRAS",
    "CREBBP", "ATM", "CHEK2", "SETD2", "STAG2", "BCOR",
})


def check_chip_gene(gene: str | None, profile_name: str = "ctdna") -> dict:
    """CT-06: hedef gen CHIP suruculerinden biriyse uyari (yalniz ctDNA)."""
    g = (gene or "").strip().upper()
    is_chip = g in CHIP_GENES
    severity = "WARN" if (is_chip and profile_name == "ctdna") else "INFO"
    note = ("Hedef CHIP surucu geni; kan hucresi kokenli varyant riski "
            "(Jaiswal 2014, Razavi 2019)." if is_chip else
            ("gen belirtilmedi" if not g else "CHIP listesinde degil"))
    return _rec("CT-06", {"gene": g or None, "is_chip_gene": is_chip},
                {"warn": "hedef gen CHIP surucu listesinde (ctDNA)"},
                severity, note)


# ---------------------------------------------------------------------------
# Amplikon ikincil yapisi (SS-04) ve primer x prob dimeri (DIM-05)
# ---------------------------------------------------------------------------
def _unimolecular_tm_c(dH_kcal: float, dS_cal: float) -> float:
    """Tek molekullu (hairpin) gecis Tm'i: Tm = dH / dS. Konsantrasyondan
    BAGIMSIZDIR (dimerden farki budur)."""
    if dS_cal >= 0:
        return float("-inf")
    return dH_kcal * 1000.0 / dS_cal - 273.15


def _structure_pairs(structure: str) -> list[tuple[int, int]]:
    """Nokta-parantez yapisindan baz ciftleri."""
    stack, pairs = [], []
    for idx, ch in enumerate(structure):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            pairs.append((stack.pop(), idx))
    return pairs


def _structure_enthalpy(seq: str, pairs: list[tuple[int, int]]) -> float:
    """MFE yapisindaki yigin basamaklarinin toplam dH'si (kcal/mol).

    (i,j) ile (i+1,j-1) ikisi de ciftse bir WC yigin basamagidir; ust
    iplikten okunur. Donguler dH=0 kaynak varsayimi tasir, bu yuzden dH
    yalniz yiginlardan gelir.
    """
    paired = set(pairs)
    nn_tab = nn.WC_SANTALUCIA_1998.require()["nn"]
    dH = 0.0
    for i, j in pairs:
        if (i + 1, j - 1) in paired:
            step = seq[i:i + 2]
            h, _ = nn_tab.get(step) or nn_tab[nn.reverse_complement(step)]
            dH += h
    return dH


def check_amplicon_structure(amplicon: str, anneal_c: float = 60.0,
                             cons: SeqPairConstraints | None = None) -> dict:
    """SS-04: tek iplik amplikon MFE'si ve en kararli yapinin Tm'si.

    Yapi ve MFE `fold.zuker_mfe` (tam DNA dongu modeli) ile bulunur. Yapi
    Tm'i PROJENIN kendi enerji modeliyle TUTARLI turetilir:
      dH   = MFE yapisindaki yigin basamaklarinin toplami (dongu dH=0),
      dG37 = MFE (dongu cezalari dahil),
      dS   = (dH - dG37) * 1000 / 310.15,
      Tm   = dH / dS  (tek molekullu; konsantrasyondan bagimsiz).
    Boylece dongu cezalari dS'e girer ve Tm cipak saptan farkli olarak
    asiri cagirmaz. Esik Ta'ya goredir (qPCR'de uzamayi engelleyen yapi).
    """
    cons = cons or SeqPairConstraints()
    from . import fold          # dongusel import'tan kacin
    s = nn.validate(amplicon)
    mfe, structure = fold.zuker_mfe(s)
    pairs = _structure_pairs(structure)
    dH = _structure_enthalpy(s, pairs)
    if dH < 0:
        dS = (dH - mfe) * 1000.0 / 310.15
        tm = dH * 1000.0 / dS - 273.15 if dS < 0 else float("-inf")
    else:
        tm = float("-inf")
    sev = "INFO"
    if tm > anneal_c:
        sev = "FAIL"
    elif tm > anneal_c - 10.0:
        sev = "WARN"
    return _rec("SS-04",
                {"amplicon_mfe_kcal_mol": round(mfe, 3),
                 "structure_tm_c": round(tm, 3) if tm != float("-inf") else None,
                 "n_base_pairs": len(pairs), "anneal_c": anneal_c},
                {"warn": f"Tm_yapi > {anneal_c - 10:.0f} C",
                 "fail": f"Tm_yapi > {anneal_c:.0f} C"}, sev,
                "Tm, MFE yapisindan tutarli turetildi (dongu cezalari dS'te)")


def check_primer_probe_dimers(fwd: str, rev: str, probe: str,
                              warn_dg37: float = -6.0, fail_dg37: float = -9.0,
                              warn_3p: float = -5.0, fail_3p: float = -7.0) -> dict:
    """DIM-05: primer x prob dimerlerinin tum kombinasyonlari.

    Her primer-prob cifti genel ('any') ve 3' uc ('3prime') modunda
    puanlanir; en kararli dG37 esikle karsilastirilir. Prob 3' ucu
    sondurucuyle blokeli oldugundan 3' esikleri primerdekinden gevsektir
    (PRB-05 ile ayni gerekce).
    """
    probe_u = nn.validate(oligomod.parse(probe).dna)
    worst_any = 0.0
    worst_3p = 0.0
    partner = None
    for name, primer in (("forward", fwd), ("reverse", rev)):
        a = align.best_dimer(nn.validate(primer), probe_u, "any")
        t = align.best_dimer(nn.validate(primer), probe_u, "3prime")
        if a.found and a.dG37_kcal < worst_any:
            worst_any = a.dG37_kcal
        if t.found and t.dG37_kcal < worst_3p:
            worst_3p = t.dG37_kcal
            partner = name
    sev = "INFO"
    if worst_any < fail_dg37 or worst_3p < fail_3p:
        sev = "FAIL"
    elif worst_any < warn_dg37 or worst_3p < warn_3p:
        sev = "WARN"
    return _rec("DIM-05",
                {"worst_any_dG37_kcal_mol": round(worst_any, 3),
                 "worst_3p_dG37_kcal_mol": round(worst_3p, 3),
                 "worst_3p_partner": partner},
                {"warn": f"dG37 < {warn_dg37} (genel) / {warn_3p} (3')",
                 "fail": f"dG37 < {fail_dg37} (genel) / {fail_3p} (3')"}, sev)


# ---------------------------------------------------------------------------
# RT primeri DNA/RNA hibrit Tm (RT-05) — INFO
# ---------------------------------------------------------------------------
def hybrid_duplex_dh_ds(dna_primer: str) -> tuple[float, float]:
    """DNA/RNA hibrit dupleksin (dH, dS)'i, Sugimoto 1995 NN tablosuyla.

    `dna_primer` RT primeridir (DNA, 5'->3'); RNA kalibi ters tamamlayicidir.
    Tablo AKTARIM BEKLIYOR; bos oldugu surece ParameterTableUnavailable.
    """
    from .params import HYBRID_SUGIMOTO_1995, HYBRID_SUGIMOTO_1995_INIT
    table = HYBRID_SUGIMOTO_1995.require()        # bos ise hata verir
    s = nn.validate(dna_primer)
    init_h, init_s = HYBRID_SUGIMOTO_1995_INIT    # baslatma terimi (bir kez)
    dH, dS = float(init_h), float(init_s)
    for i in range(len(s) - 1):
        step = s[i:i + 2]
        h, ds = table[step]
        dH += h
        dS += ds
    return dH, dS


def check_rt05(dna_primer: str, total_strand_M: float = 2.0e-7) -> dict:
    """RT-05: RT primerinin RNA kalibina hibrit Tm'i (INFO).

    Tablo aktarilmamissa sessizce 'pending' raporlar (hata firlatmaz), boylece
    RT profili tablonun gelmesini beklerken calismaya devam eder.
    """
    from .params import HYBRID_SUGIMOTO_1995, ParameterTableUnavailable
    try:
        dH, dS = hybrid_duplex_dh_ds(dna_primer)
    except ParameterTableUnavailable:
        return _rec("RT-05",
                    {"hybrid_tm_c": None, "status": "pending"},
                    {}, "INFO",
                    "DNA/RNA hibrit NN tablosu aktarilmadi "
                    f"({HYBRID_SUGIMOTO_1995.provenance.citation})")
    denom = dS + 1.9872 * math.log(total_strand_M / 4.0)
    tm = dH * 1000.0 / denom - 273.15 if denom < 0 else float("-inf")
    return _rec("RT-05",
                {"hybrid_tm_c": round(tm, 3) if tm != float("-inf") else None,
                 "status": "evaluated"},
                {}, "INFO", "DNA/RNA hibrit Tm (Sugimoto 1995); bilgi kaydi")


# ---------------------------------------------------------------------------
# CT-05: WT blocker / ddPCR prob ayrimi
# ---------------------------------------------------------------------------
def check_ct05(mutant_probe: str, wt_probe: str,
               cond: salt.Conditions | None = None,
               salt_model: str = "owczarzy_2008",
               delta_tm_warn: float = 5.0, delta_tm_fail: float = 3.0) -> dict:
    """CT-05: mutant ve WT prob arasindaki Tm farki (kanal ayrimi).

    Problar genellikle LNA tasir; LNA mutlak Tm'i belirgin yukseltir ve bu
    motorda MODELLENMEZ (ayri bir McTigue 2004 parametre tablosu gerekir).
    Bu yuzden raporlanan Tm'ler DNA ISKELETININ degeridir (gercek Tm'in alt
    siniri); iki probun LNA sayilari da verilir ki fark yorumlanabilsin.
    Kucuk ayrim (ΔTm) iki kanalin birbirine girmesi demektir.

    Girdi modifiye notasyonda olabilir (+LNA, /uc/); iceride ayristirilir.
    """
    from . import oligomod
    cond = cond or salt.Conditions()
    mp = oligomod.parse(mutant_probe)
    wp = oligomod.parse(wt_probe)
    tm_mut = concentration.build(mp.dna, cond=cond, salt_model=salt_model).tm_c
    tm_wt = concentration.build(wp.dna, cond=cond, salt_model=salt_model).tm_c
    delta = abs(tm_mut - tm_wt)
    sev = "INFO"
    if delta < delta_tm_fail:
        sev = "FAIL"
    elif delta < delta_tm_warn:
        sev = "WARN"
    return _rec("CT-05",
                {"mutant_probe_dna": mp.dna, "wt_probe_dna": wp.dna,
                 "mutant_tm_c_dna_backbone": round(tm_mut, 3),
                 "wt_tm_c_dna_backbone": round(tm_wt, 3),
                 "delta_tm_c": round(delta, 3),
                 "mutant_lna_count": mp.lna_count, "wt_lna_count": wp.lna_count,
                 "lna_modelled": False},
                {"warn": f"ΔTm < {delta_tm_warn} C", "fail": f"ΔTm < {delta_tm_fail} C"},
                sev,
                "Tm DNA iskeletinden; LNA stabilizasyonu modellenMEDI (alt "
                "sinir). Mutlak esik LNA-dogru Tm ile kalibre edilmeli (REP-03)")


# ---------------------------------------------------------------------------
# SPC-08: mitokondri / alt contig / rRNA ayrica taranir
# ---------------------------------------------------------------------------
def scan_contig_sites(primer: str, contig_seq: str, max_mismatch: int = 2,
                      anchor: int = 3) -> list[dict]:
    """Primeri bir contig'in iki ipliginde kayar pencereyle tarar.

    Kucuk contig'ler (chrM ~16.5 kb, rDNA) icin saf Python yeterlidir; buyuk
    ALT contig'leri icin C++ FM-index gerekir (WSL scripti ayrica verilir).
    `max_mismatch` ustundeki yerler atilir. 3'-capa = son `anchor` nt tam
    eslesmeli (uzatilabilir off-target).
    """
    s = nn.validate(primer)
    n = len(s)
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    rc = "".join(comp.get(b, "N") for b in reversed(contig_seq.upper()))
    hits: list[dict] = []
    for strand, seqc in (("+", contig_seq.upper()), ("-", rc)):
        for i in range(len(seqc) - n + 1):
            window = seqc[i:i + n]
            mm = 0
            for a, b in zip(s, window):
                if a != b:
                    mm += 1
                    if mm > max_mismatch:
                        break
            if mm <= max_mismatch:
                anchored = all(s[-k] == window[-k] for k in range(1, anchor + 1))
                hits.append({"strand": strand, "pos": i, "mismatches": mm,
                             "three_prime_anchored": anchored})
    return hits


def check_spc08(primer: str, contigs: dict[str, str],
                max_mismatch: int = 2, anchor: int = 3) -> dict:
    """SPC-08: primer mitokondri/alt/rRNA contig'lerinde baglaniyor mu.

    `contigs` {ad: dizi}. 3'-capali site FAIL (uzatilabilir), capasiz
    uyumsuzluklu site WARN. Contig verilmezse INFO/'taranmadi'.
    """
    if not contigs:
        return _rec("SPC-08", {"scanned": False}, {}, "INFO",
                    "mitokondri/alt/rRNA contig verilmedi")
    summary: dict[str, dict] = {}
    anchored_total = 0
    site_total = 0
    for name, seqc in contigs.items():
        hits = scan_contig_sites(primer, seqc, max_mismatch, anchor)
        anch = sum(1 for h in hits if h["three_prime_anchored"])
        summary[name] = {"sites": len(hits), "anchored": anch}
        anchored_total += anch
        site_total += len(hits)
    sev = "FAIL" if anchored_total else ("WARN" if site_total else "INFO")
    return _rec("SPC-08",
                {"scanned": True, "per_contig": summary,
                 "anchored_sites": anchored_total, "total_sites": site_total,
                 "max_mismatch": max_mismatch, "anchor_nt": anchor},
                {"warn": "contig'te uyumsuzluklu baglanma",
                 "fail": "contig'te 3'-capali (uzatilabilir) baglanma"}, sev,
                "kucuk contig'ler Python'da taranir; buyuk ALT icin C++ gerekir")


SEVERITY_ORDER = {"INFO": 0, "WARN": 1, "FAIL": 2}


def check_seq01(raw_seq: str) -> dict:
    """SEQ-01: alfabe + modifikasyon notasyonu.

    Yalniz A/C/G/T beklenir; IUPAC dejenere kodlar yalniz dejenere modda
    kabul edilir; LNA (+N) ve uc modifikasyonlari (/.../) oligomod ile
    ayristirilir ve termodinamik hesap DNA iskeleti uzerinden yapilir.
    LNA / uc modifikasyonu kararlilik katkisi bu motorda MODELLENMEZ ve bu
    durum acikca raporlanir. Taninmayan karakter FAIL dondurur.
    """
    parsed = oligomod.parse(raw_seq)
    core = re.sub(r"/[^/]+/", "", raw_seq.strip()).replace("+", "")
    bases = [c.upper() for c in core if c.isalpha()]
    iupac = sorted({b for b in bases if b in IUPAC and b not in "ACGT"})
    invalid = sorted({b for b in bases if b not in IUPAC})
    has_mod = bool(parsed.lna_positions or parsed.end_mods)
    sev = "FAIL" if invalid else "INFO"
    return _rec(
        "SEQ-01",
        {"length_nt": len(parsed.dna),
         "pure_acgt": not iupac and not invalid and not has_mod,
         "iupac_degenerate": iupac,
         "modifications": {"lna_count": parsed.lna_count,
                           "lna_positions": sorted(parsed.lna_positions),
                           "end_mods": list(parsed.end_mods)},
         "invalid_chars": invalid},
        {"fail": "A/C/G/T, IUPAC ya da taninan modifikasyon disi karakter"},
        sev,
        "LNA/mod ayristirildi; termodinamik DNA iskeleti uzerinden "
        "(LNA/mod kararlilik katkisi modellenmez)")


def check_ss05(probe: str, anneal_c: float = 60.0,
               cons: SeqPairConstraints | None = None) -> dict:
    """SS-05: prob hairpin Tm'i ve yapisi.

    SS-04 ile ayni enerji modeli (fold.zuker_mfe, tam DNA dongu modeli) prob
    dizisine uygulanir. LNA / uc modifikasyonlari oligomod ile ayristirilir
    ve DNA iskeleti kullanilir; LNA hairpin'i KARARLASTIRIR fakat katki
    MODELLENMEZ, bu nedenle bildirilen Tm bir ALT SINIRDIR. Esik Ta'ya
    goredir (uzamayi/sondurucu acilmasini engelleyen prob yapisi).
    """
    cons = cons or SeqPairConstraints()
    from . import fold          # dongusel import'tan kacin
    parsed = oligomod.parse(probe)
    seq = nn.validate(parsed.dna)
    mfe, structure = fold.zuker_mfe(seq)
    pairs = _structure_pairs(structure)
    dH = _structure_enthalpy(seq, pairs)
    if dH < 0:
        dS = (dH - mfe) * 1000.0 / 310.15
        tm = dH * 1000.0 / dS - 273.15 if dS < 0 else float("-inf")
    else:
        tm = float("-inf")
    sev = "INFO"
    if tm > anneal_c:
        sev = "FAIL"
    elif tm > anneal_c - 10.0:
        sev = "WARN"
    note = ("Tm MFE yapisindan tutarli turetildi (dongu cezalari dS'te); "
            "prob 3' ucu genelde sondurucu/blokeli")
    if parsed.lna_count:
        note += f"; {parsed.lna_count} LNA -> bildirilen Tm ALT SINIR"
    return _rec(
        "SS-05",
        {"probe_mfe_kcal_mol": round(mfe, 3),
         "structure_tm_c": round(tm, 3) if tm != float("-inf") else None,
         "n_base_pairs": len(pairs), "anneal_c": anneal_c,
         "lna_count": parsed.lna_count},
        {"warn": f"Tm_yapi > {anneal_c - 10:.0f} C",
         "fail": f"Tm_yapi > {anneal_c:.0f} C"}, sev, note)


def check_dim07(fwd: str, rev: str, fwd_tail: str = "", rev_tail: str = "",
                warn_dg37: float = -6.0, fail_dg37: float = -9.0,
                warn_3p: float = -5.0, fail_3p: float = -7.0,
                cons: SeqPairConstraints | None = None) -> dict:
    """DIM-07: kuyruk-dahil dimer analizi (M13 / NGS adaptor kuyruklari).

    Evrensel kuyruklar primerlerin 5' ucuna eklenir ve self + cross dimerler
    (genel ve 3' uc modunda) kuyruk DAHIL yeniden hesaplanir. Kuyruk dahil en
    kararli dG37 esikle karsilastirilir; kuyruksuz degerle farki da raporlanir
    (kuyruk kaynakli dimer riski). Kuyruk verilmezse kuyruksuz degerlerle ayni
    sonuc doner (kurul izlenebilirligi icin kayit yine uretilir).
    """
    cons = cons or SeqPairConstraints()

    def worst_set(a: str, b: str) -> tuple[float, float, str]:
        wa = w3 = 0.0
        partner = None
        for name, x, y in (("F-F", a, a), ("R-R", b, b), ("F-R", a, b)):
            g = align.best_dimer(x, y, "any")
            t = align.best_dimer(x, y, "3prime")
            if g.found and g.dG37_kcal < wa:
                wa = g.dG37_kcal
            if t.found and t.dG37_kcal < w3:
                w3 = t.dG37_kcal
                partner = name
        return wa, w3, partner

    f_full = nn.validate(fwd_tail + fwd)
    r_full = nn.validate(rev_tail + rev)
    tail_any, tail_3p, partner = worst_set(f_full, r_full)
    base_any, base_3p, _ = worst_set(nn.validate(fwd), nn.validate(rev))
    sev = "INFO"
    if tail_any < fail_dg37 or tail_3p < fail_3p:
        sev = "FAIL"
    elif tail_any < warn_dg37 or tail_3p < warn_3p:
        sev = "WARN"
    return _rec(
        "DIM-07",
        {"tail_included_any_dG37_kcal_mol": round(tail_any, 3),
         "tail_included_3p_dG37_kcal_mol": round(tail_3p, 3),
         "worst_3p_pair": partner,
         "no_tail_any_dG37_kcal_mol": round(base_any, 3),
         "no_tail_3p_dG37_kcal_mol": round(base_3p, 3),
         "delta_any_from_tail_kcal_mol": round(tail_any - base_any, 3),
         "fwd_tail": fwd_tail, "rev_tail": rev_tail},
        {"warn": f"dG37 < {warn_dg37} (genel) / {warn_3p} (3')",
         "fail": f"dG37 < {fail_dg37} (genel) / {fail_3p} (3')"}, sev,
        "M13/NGS adaptor kuyruklari dahil self+cross dimer (NGS-01)")


def check_seq06(header: str, fasta_sha256: str,
                n_bases: int | None = None) -> dict:
    """SEQ-06: referans kimligi. Build + yama + FASTA SHA-256 raporlanir.

    Izlenebilirlik kaydi (REP ailesi): tasarimda kullanilan referans FASTA'nin
    build/yama etiketi (FASTA basligindan ayristirilir) ve tam dosya SHA-256'si
    beyan edilir. Checksum, farkli bir referansla yapilan tasarimin ayni QC
    ciktisini vermeyecegini garanti eder.
    """
    h = header or ""
    build = patch = None
    m = re.search(r"GRCh3\d(?:\.p\d+)?", h)
    if m:
        build = m.group(0).split(".")[0]
        pm = re.search(r"p(\d+)", m.group(0))
        patch = ("p" + pm.group(1)) if pm else None
    sev = "INFO" if fasta_sha256 else "WARN"
    return _rec(
        "SEQ-06",
        {"header": h[:160], "build": build, "patch": patch,
         "fasta_sha256": fasta_sha256, "bases": n_bases},
        {"warn": "FASTA SHA-256 hesaplanamadi"}, sev,
        "referans izlenebilirligi (build + yama + checksum)")


def worst(records: list[dict]) -> str:
    return max((r["severity"] for r in records),
               key=lambda s: SEVERITY_ORDER[s], default="INFO")
