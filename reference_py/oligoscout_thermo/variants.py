"""Varyant katmani: VAR-01..VAR-08.

Primer (ve prob) baglanma yerine dusen bir varyant, o alelde baglanmayi
zayiflatir ya da tamamen engeller: **alel dusmesi (allele dropout)**. Etki
varyantin primer uzerindeki KONUMUNA ve popülasyondaki SIKLIGINA baglidir;
3' uca yakin bir varyant uzamayi dogrudan bloke eder.

Bu modul ag erisimi YAPMAZ; varyantlari disaridan alir (CLI gnomAD v4'ten
ceker ve onbellege yazar), boylece referans katmani standart kutuphane
disina cikmaz ve testler cevrimdisi calisir.

ESIKLERIN KAYNAGI: asagidaki esikler QC katalogundaki VAR maddelerinin
literatur baslangic degerleridir ve KALIBRE DEGILDIR. Katalogun (qc_rules.yaml)
kendi esikleri elinize gectiginde `VariantConstraints` alanlarini degistirmek
yeterlidir; kural kimlikleri ve kayit semasi ayni kalir.

  VAR-01  Baglanma yerinde herhangi bir varyant (sayim + en yuksek AF)
  VAR-02  3' terminal pencerede (son 5 nt) varyant — uzamayi bloke eder
  VAR-03  Baglanma yerinde YAYGIN varyant (AF >= %1)
  VAR-04  Tekli ic SNV icin alternatif eksi referans dG37
  VAR-05  Prob bolgesinde varyant
  VAR-06  Baglanma yeriyle ortusen indel — uzunluk degisimi, dizi kaymasi
  VAR-07  Popülasyona ozgu yuksek frekans (kuresel AF'nin cok ustunde)
  VAR-08  Amplikon icinde (primer disi) uzunluk degistiren yaygin varyant
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

from .engine import ENGINE_VERSION, SEVERITY_ORDER, CheckRecord
from . import nn, params


@dataclass(frozen=True)
class Variant:
    """Tek bir varyant kaydi. Konum 1-tabanli (VCF/gnomAD kurali)."""

    chrom: str
    pos: int
    ref: str
    alt: str
    af: float
    ac: int = 0
    an: int = 0
    rsid: str = ""
    source: str = "gnomad_r4"
    pop_af: dict[str, float] = field(default_factory=dict)

    @property
    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt)

    @property
    def length_change(self) -> int:
        return len(self.alt) - len(self.ref)

    @property
    def span(self) -> tuple[int, int]:
        """Etkilenen genomik aralik, 0-tabanli yari acik."""
        return (self.pos - 1, self.pos - 1 + max(1, len(self.ref)))

    def pop_max(self) -> tuple[str, float]:
        if not self.pop_af:
            return ("", 0.0)
        pid = max(self.pop_af, key=lambda k: self.pop_af[k])
        return (pid, self.pop_af[pid])

    def label(self) -> str:
        return self.rsid or f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


@dataclass(frozen=True)
class Site:
    """Bir primerin (ya da probun) genomdaki baglanma yeri."""

    site_id: str
    chrom: str
    start: int            # 0-tabanli
    end: int              # dislayici
    strand: str           # "+" | "-"
    sequence: str = ""
    # Primer yonundeki referans. None tam eslesmeyi, "" dizinin yoklugunu belirtir.
    template_site: str | None = None

    def three_prime_window(self, nt: int) -> tuple[int, int]:
        """3' uctaki son `nt` bazin genomik araligi (0-tabanli yari acik).

        + iplikte 3' uc sagdadir, - iplikte solda. Yanlis tarafi almak
        VAR-02'yi anlamsiz kilar, bu yuzden ipligi acikca kullaniyoruz.
        """
        if self.strand == "+":
            return (max(self.start, self.end - nt), self.end)
        return (self.start, min(self.end, self.start + nt))


@dataclass(frozen=True)
class VariantConstraints:
    af_report: float = 1e-4      # bunun altindaki varyantlar yalniz sayilir
    af_warn: float = 1e-3
    af_common: float = 0.01      # "yaygin" taniminin esigi
    af_high: float = 0.05
    three_prime_nt: int = 5      # VAR-02 penceresi
    pop_ratio: float = 5.0       # VAR-07: pop AF / kuresel AF


# ---------------------------------------------------------------------------
# Yardimcilar
# ---------------------------------------------------------------------------
def normalize_chrom(name: str) -> str:
    """'chr22' <-> '22' farkini tek bicime indirger."""
    return name[3:] if name.lower().startswith("chr") else name


def overlapping(variants: Iterable[Variant], chrom: str, start: int,
                end: int) -> list[Variant]:
    c = normalize_chrom(chrom)
    out = []
    for v in variants:
        if normalize_chrom(v.chrom) != c:
            continue
        vs, ve = v.span
        if vs < end and ve > start:
            out.append(v)
    return sorted(out, key=lambda v: v.pos)


def distance_from_3p(site: Site, v: Variant) -> int:
    """Varyantin primerin 3' ucuna uzakligi (0 = terminal baz)."""
    vs, _ = v.span
    return (site.end - 1 - vs) if site.strand == "+" else (vs - site.start)


def max_af(variants: Iterable[Variant]) -> float:
    return max((v.af for v in variants), default=0.0)


def _worst(variants: Iterable[Variant]) -> Variant | None:
    vs = list(variants)
    return max(vs, key=lambda v: v.af) if vs else None


def _detail(v: Variant | None, site: Site | None = None) -> dict[str, Any]:
    if v is None:
        return {}
    d = {"variant": v.label(), "af": v.af, "an": v.an, "source": v.source}
    if site is not None:
        d["distance_from_3p_nt"] = distance_from_3p(site, v)
    pid, paf = v.pop_max()
    if pid:
        d["pop_max"] = pid
        d["pop_max_af"] = round(paf, 6)
    return d


def _variant_thermodynamics(site: Site, hits: list[Variant]) -> dict[str, Any]:
    reference = site.sequence if site.template_site is None else site.template_site
    reference = reference.upper()
    if not hits:
        return {"evaluated": False, "reason": "baglanma yerinde varyant yok"}
    if len(reference) != site.end - site.start or site.strand not in ("+", "-"):
        return {"evaluated": False, "reason": "referans baglanma dizisi yok"}
    try:
        ref_dg = nn.duplex_with_mismatches(site.sequence, reference).dG_kcal(37.0)
    except (ValueError, params.ParameterTableUnavailable) as error:
        return {"evaluated": False, "reason": str(error)}

    details = []
    for variant in hits:
        row = _detail(variant, site)
        details.append(row)
        ref, alt = variant.ref.upper(), variant.alt.upper()
        if (len(ref) != 1 or len(alt) != 1 or ref == alt
                or ref not in nn.VALID or alt not in nn.VALID):
            row["reason"] = "yalniz tek baz degisimi desteklenir"
            continue
        index = variant.pos - 1 - site.start
        if site.strand == "-":
            index = len(reference) - 1 - index
            ref, alt = nn.reverse_complement(ref), nn.reverse_complement(alt)
        if not 0 <= index < len(reference) or reference[index] != ref:
            row["reason"] = "varyant REF bazi baglanma dizisiyle eslesmiyor"
            continue
        alternate = reference[:index] + alt + reference[index + 1:]
        try:
            alt_dg = nn.duplex_with_mismatches(site.sequence, alternate).dG_kcal(37.0)
        except (ValueError, params.ParameterTableUnavailable) as error:
            row["reason"] = str(error)
            continue
        row["delta_delta_g37_kcal_mol"] = round(alt_dg - ref_dg, 4)

    scored = sum("delta_delta_g37_kcal_mol" in row for row in details)
    return {"evaluated": bool(scored), "n_scored": scored,
            "n_unscored": len(details) - scored, "variants": details}


# ---------------------------------------------------------------------------
# Kontroller
# ---------------------------------------------------------------------------
def evaluate_site(site: Site, variants: Iterable[Variant],
                  cons: VariantConstraints | None = None,
                  profile_name: str = "", conditions: dict | None = None,
                  probe: bool = False) -> list[CheckRecord]:
    """Bir baglanma yeri icin VAR-01..VAR-07 (prob ise VAR-05) kayitlari."""
    cons = cons or VariantConstraints()
    conditions = conditions or {}
    digest = params.registry_digest()
    hits = overlapping(variants, site.chrom, site.start, site.end)
    checks: list[CheckRecord] = []

    def rec(cid, value, threshold, severity, model, note=""):
        checks.append(CheckRecord(
            id=cid, target=site.site_id, value=value, threshold=threshold,
            severity=severity, model=model, conditions=conditions,
            param_set_sha256=digest, profile=profile_name, note=note))

    source = next(iter({v.source for v in hits}), "gnomad_r4") if hits else "gnomad_r4"
    model = {"variant_source": source,
             "coordinates": "0-tabanli yari acik, GRCh38",
             "threshold_origin": "literatur baslangic degeri, kalibre degil"}

    if probe:
        # --- VAR-05 prob bolgesi ---------------------------------------
        worst = _worst(hits)
        af = max_af(hits)
        sev = ("FAIL" if af >= cons.af_common
               else "WARN" if af >= cons.af_warn else "INFO")
        rec("VAR-05",
            {"n_variants": len(hits), "max_af": round(af, 6), **_detail(worst)},
            {"warn": f"AF >= {cons.af_warn}", "fail": f"AF >= {cons.af_common}"},
            sev, model,
            "prob bolgesindeki varyant sinyali dusurur, uzamayi engellemez")
        return checks

    # --- VAR-01 herhangi bir varyant -----------------------------------
    worst = _worst(hits)
    af = max_af(hits)
    sev = "WARN" if af >= cons.af_warn else "INFO"
    rec("VAR-01",
        {"n_variants": len(hits),
         "n_above_report": sum(1 for v in hits if v.af >= cons.af_report),
         "max_af": round(af, 6), **_detail(worst, site)},
        {"warn": f"AF >= {cons.af_warn}"}, sev, model,
        "baglanma yerinin tamami tarandi")

    # --- VAR-02 3' terminal pencere ------------------------------------
    lo, hi = site.three_prime_window(cons.three_prime_nt)
    tp = overlapping(hits, site.chrom, lo, hi)
    worst_tp = _worst(tp)
    af_tp = max_af(tp)
    sev = ("FAIL" if af_tp >= cons.af_common
           else "WARN" if af_tp >= cons.af_report else "INFO")
    rec("VAR-02",
        {"window_nt": cons.three_prime_nt, "n_variants": len(tp),
         "max_af": round(af_tp, 6), **_detail(worst_tp, site)},
        {"warn": f"AF >= {cons.af_report}", "fail": f"AF >= {cons.af_common}"},
        sev, model,
        "3' uca yakin uyumsuzluk uzamayi dogrudan engeller")

    # --- VAR-03 yaygin varyant -----------------------------------------
    common = [v for v in hits if v.af >= cons.af_common]
    worst_c = _worst(common)
    af_c = max_af(common)
    sev = ("FAIL" if af_c >= cons.af_high
           else "WARN" if common else "INFO")
    rec("VAR-03",
        {"n_common": len(common), "max_af": round(af_c, 6),
         **_detail(worst_c, site)},
        {"warn": f"AF >= {cons.af_common}", "fail": f"AF >= {cons.af_high}"},
        sev, model, "yaygin varyant kohortun buyuk kisminda alel dusmesi demek")

    # --- VAR-04 termodinamik etki --------------------------------------
    mismatch_table = params.MISMATCH_ALLAWI_PEYRET.provenance
    rec("VAR-04", _variant_thermodynamics(site, hits), {}, "INFO",
        {"nn_table": mismatch_table.table_id,
         "nn_table_status": mismatch_table.verification_status,
         "temperature_c": 37.0, "reference_conditions": mismatch_table.conditions},
        "dG37(alternatif) - dG37(referans); kalibre karar esigi yok")

    # --- VAR-06 indel ----------------------------------------------------
    indels = [v for v in hits if v.is_indel]
    scored = [v for v in indels if v.af >= cons.af_report]
    worst_i = _worst(scored) or _worst(indels)
    af_i = max_af(scored)
    # Tek gozlemlik (singleton) indeller her bolgede bulunur; bunlari WARN
    # saymak butun paneli isaretler ve olcumu anlamsiz kilar. Esik af_report.
    sev = ("FAIL" if af_i >= cons.af_common
           else "WARN" if scored else "INFO")
    rec("VAR-06",
        {"n_indels": len(indels), "n_above_report": len(scored),
         "max_af": round(af_i, 6),
         "length_change": worst_i.length_change if worst_i else 0,
         **_detail(worst_i, site)},
        {"warn": f"indel AF >= {cons.af_report}",
         "fail": f"indel AF >= {cons.af_common}"},
        sev, model, "indel dizi kaymasi yaratir; tek baz uyumsuzlugundan agirdir")

    # --- VAR-07 popülasyona ozgu frekans ---------------------------------
    flagged = []
    for v in hits:
        pid, paf = v.pop_max()
        if paf >= cons.af_common and (v.af <= 0 or paf >= cons.pop_ratio * v.af):
            flagged.append((v, pid, paf))
    sev = "WARN" if flagged else "INFO"
    top = max(flagged, key=lambda t: t[2]) if flagged else None
    rec("VAR-07",
        {"n_population_specific": len(flagged),
         "population": top[1] if top else "",
         "pop_af": round(top[2], 6) if top else 0.0,
         "global_af": top[0].af if top else 0.0,
         "variant": top[0].label() if top else ""},
        {"warn": f"pop AF >= {cons.af_common} ve >= {cons.pop_ratio}x kuresel AF"},
        sev, model,
        "belirli atalarda alel dusmesi; kuresel AF bunu gizler")

    return checks


def evaluate_amplicon(amplicon_id: str, chrom: str, start: int, end: int,
                      sites: Iterable[Site], variants: Iterable[Variant],
                      cons: VariantConstraints | None = None,
                      profile_name: str = "",
                      conditions: dict | None = None) -> list[CheckRecord]:
    """VAR-08: amplikon icinde, primerlerin DISINDA kalan yaygin indeller."""
    cons = cons or VariantConstraints()
    inner = overlapping(variants, chrom, start, end)
    site_spans = [(s.start, s.end) for s in sites]

    def in_primer(v: Variant) -> bool:
        vs, ve = v.span
        return any(vs < b and ve > a for a, b in site_spans)

    inner = [v for v in inner if not in_primer(v)]
    indels = [v for v in inner if v.is_indel and v.af >= cons.af_common]
    worst = _worst(indels)
    sev = "WARN" if indels else "INFO"
    return [CheckRecord(
        id="VAR-08", target=amplicon_id,
        value={"amplicon_bp": end - start,
               "n_variants_inside": len(inner),
               "n_common_indels": len(indels),
               "max_length_change": (max((abs(v.length_change) for v in indels),
                                         default=0)),
               **_detail(worst)},
        threshold={"warn": f"primer disi indel AF >= {cons.af_common}"},
        severity=sev,
        model={"variant_source": "gnomad_r4",
               "rule": "primer baglanma yerlerine dusen varyantlar haric"},
        conditions=conditions or {}, param_set_sha256=params.registry_digest(),
        profile=profile_name,
        note="urun boyu degisir; boy ayrimina dayanan okumalar etkilenir")]


VERDICT_FROM_SEVERITY = {"INFO": "PASS", "WARN": "WARN", "FAIL": "FAIL"}
VERDICT_ORDER = {"PASS": 0, "WARN": 1, "FAIL": 2}


def verdict_from_severities(severities: Iterable[str]) -> str:
    """VAR kontrol ciddiyetleri -> cift duzeyi varyant karari."""
    out = "PASS"
    for s in severities:
        cand = VERDICT_FROM_SEVERITY[s]
        if VERDICT_ORDER[cand] > VERDICT_ORDER[out]:
            out = cand
    return out


def combine_verdicts(*verdicts: str) -> str:
    """Katmanlarin EN KOTUSU. Ozgullukte PASS + varyantta FAIL = FAIL.

    Iki katman birbirinin yerine gecmez: biri "baska yerde urun var mi",
    digeri "hedefte baglanma her bireyde tutuyor mu" sorusunu yanitlar.
    """
    out = "PASS"
    for v in verdicts:
        if not v:
            continue
        if VERDICT_ORDER[v] > VERDICT_ORDER[out]:
            out = v
    return out


def worst_severity(checks: Iterable[CheckRecord]) -> str:
    return max((c.severity for c in checks), key=lambda s: SEVERITY_ORDER[s],
               default="INFO")


def summarize(checks: Iterable[CheckRecord]) -> dict[str, Any]:
    checks = list(checks)
    by_id: dict[str, str] = {}
    for c in checks:
        cur = by_id.get(c.id, "INFO")
        by_id[c.id] = max(cur, c.severity, key=lambda s: SEVERITY_ORDER[s])
    return {"worst_severity": worst_severity(checks), "by_check": by_id,
            "engine_version": ENGINE_VERSION}
