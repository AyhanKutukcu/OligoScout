"""Hidroliz probu (TaqMan / ddPCR) kontrolleri: PRB-01..PRB-05.

Kapsam: amplikon icine oturan, 5' florofor + 3' sonlandirici tasiyan hidroliz
probu. Bes kontrol QC katalogundan gelir.

  PRB-01  Prob Tm'i primer Tm'lerinden 6-10 C yuksek olmali.
          Gerekce: prob, primer uzamasi baslamadan once hedefe baglanmis
          olmali; aksi halde polimeraz probu yerinden surer ve sinyal cikmaz.
  PRB-02  5' uctaki baz G olmamali.
          Gerekce: 5' G, sonlandirici kesildikten sonra bile floroforu
          sondurur; raporlanan sinyal bastirilir.
  PRB-03  Probun yaziminda C sayisi G sayisindan fazla olmali.
          Gerekce: G agirlikli iplik hem sondurme hem de G-dortlusu egilimi
          tasir; iki iplikten C agirlikli olani secilir.
  PRB-04  Prob, primer baglanma yerleriyle CAKISMAMALI ve amplikonun icinde
          kalmali. Cakisma, prob ile primerin ayni yeri paylasmasi demektir.
  PRB-05  Prob self-dimeri ve prob x primer dimerleri.
          NOT: probun 3' ucu bloke oldugu icin buradaki risk uzama degil,
          HIBRIDIZASYON REKABETIDIR; dimer, probu hedeften uzak tutar.

PRB-05, `align` modulunden tam eslesen ve tekli ic uyumsuzluklu dimerleri
kullanir. Esikler literatur kilavuz degeridir, KALIBRE DEGILDIR.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from . import align, concentration, nn, params, salt
from .design import longest_homopolymer
from .engine import ENGINE_VERSION, SEVERITY_ORDER, CheckRecord


@dataclass(frozen=True)
class ProbeConstraints:
    """Prob tasarim kisitlari ve PRB esikleri; profil tarafindan ezilebilir."""

    length_min: int = 13
    length_max: int = 30
    tm_delta_min: float = 6.0        # PRB-01 pencere alt siniri
    tm_delta_max: float = 10.0       # PRB-01 pencere ust siniri
    tm_delta_fail_lo: float = 3.0    # bunun altinda FAIL
    tm_delta_fail_hi: float = 13.0   # bunun ustunde FAIL
    gc_min: float = 0.30
    gc_max: float = 0.80
    max_homopolymer: int = 4         # 4 ve uzeri ayni baz elenir
    forbid_g_run: int = 4            # GGGG yasak
    dimer_dg_warn: float = -6.0      # PRB-05
    dimer_dg_fail: float = -9.0
    dimer3p_dg_warn: float = -5.0
    dimer3p_dg_fail: float = -7.0


@dataclass(frozen=True)
class ProbeCandidate:
    sequence: str
    strand: str          # "+" hedefin verilen ipligi, "-" tumleyeni
    start: int           # amplikon icinde 0-tabanli
    end: int             # dislayici
    tm_c: float
    gc: float
    delta_tm_c: float    # prob Tm - en yuksek primer Tm
    g_count: int
    c_count: int
    score: float = 0.0

    def as_row(self) -> list:
        return [self.sequence, self.strand, self.start, self.end,
                len(self.sequence), round(100 * self.gc, 2),
                round(self.tm_c, 3), round(self.delta_tm_c, 3),
                self.g_count, self.c_count, round(self.score, 4)]


# ---------------------------------------------------------------------------
# Kompozisyon yardimcilari
# ---------------------------------------------------------------------------
def gc_counts(seq: str) -> tuple[int, int]:
    s = nn.validate(seq)
    return s.count("G"), s.count("C")


def composition_reject(seq: str, cons: ProbeConstraints) -> tuple[str, str] | None:
    """Tasarim asamasi sert filtreleri. Gecerse None, yoksa (kural, deger)."""
    s = nn.validate(seq)
    if not (cons.length_min <= len(s) <= cons.length_max):
        return ("PRB-00", f"uzunluk {len(s)} nt")
    gc = nn.gc_fraction(s)
    if not (cons.gc_min <= gc <= cons.gc_max):
        return ("CMP-01", f"GC %{100 * gc:.1f}")
    if longest_homopolymer(s) >= cons.max_homopolymer:
        return ("CMP-03", f"homopolimer {longest_homopolymer(s)}")
    if "G" * cons.forbid_g_run in s:
        return ("CMP-03", "poli-G")
    if s[0] == "G":
        return ("PRB-02", "5' uc G")
    g, c = gc_counts(s)
    if c <= g:
        return ("PRB-03", f"G={g} C={c}")
    return None


# ---------------------------------------------------------------------------
# PRB-01..PRB-05 degerlendirmesi
# ---------------------------------------------------------------------------
def evaluate_probe(
    probe_id: str,
    probe: str,
    forward: str,
    reverse: str,
    amplicon: str | None = None,
    cond: salt.Conditions | None = None,
    cons: ProbeConstraints | None = None,
    profile_name: str = "ddpcr",
    salt_model: str = "owczarzy_2008",
) -> dict[str, Any]:
    """Bir probu bes PRB kontrolunden gecirir, REP-02 kanit kaydi uretir."""
    cond = cond or salt.Conditions()
    cons = cons or ProbeConstraints()
    p = nn.validate(probe)
    f = nn.validate(forward)
    r = nn.validate(reverse)
    digest = params.registry_digest()
    checks: list[CheckRecord] = []

    def rec(cid, value, threshold, severity, model, note=""):
        checks.append(CheckRecord(
            id=cid, target=probe_id, value=value, threshold=threshold,
            severity=severity, model=model, conditions=cond.as_dict(),
            param_set_sha256=digest, profile=profile_name, note=note))

    m_p = concentration.build(p, cond=cond, salt_model=salt_model)
    m_f = concentration.build(f, cond=cond, salt_model=salt_model)
    m_r = concentration.build(r, cond=cond, salt_model=salt_model)
    tm_primer_max = max(m_f.tm_c, m_r.tm_c)
    delta = m_p.tm_c - tm_primer_max

    # --- PRB-01  prob Tm'i primerlerin uzerinde mi ---------------------------
    if delta < cons.tm_delta_fail_lo or delta > cons.tm_delta_fail_hi:
        sev = "FAIL"
    elif delta < cons.tm_delta_min or delta > cons.tm_delta_max:
        sev = "WARN"
    else:
        sev = "INFO"
    rec("PRB-01",
        {"tm_probe_c": round(m_p.tm_c, 3),
         "tm_forward_c": round(m_f.tm_c, 3),
         "tm_reverse_c": round(m_r.tm_c, 3),
         "delta_tm_c": round(delta, 3)},
        {"warn": f"{cons.tm_delta_min}-{cons.tm_delta_max} C disi",
         "fail": f"<{cons.tm_delta_fail_lo} veya >{cons.tm_delta_fail_hi} C"},
        sev,
        {"reference": "prob Tm - en yuksek primer Tm",
         "conc_model": m_p.conc_model, "salt_model": m_p.salt_model},
        "prob, uzama baslamadan baglanmis olmali")

    # --- PRB-02  5' uctaki G -------------------------------------------------
    sev = "FAIL" if p[0] == "G" else "INFO"
    rec("PRB-02",
        {"five_prime_base": p[0],
         "five_prime_base_other_strand": nn.reverse_complement(p)[0]},
        {"fail": "5' uc G"},
        sev,
        {"mechanism": "5' G florofor sondurmesi"},
        "G ise tumleyen ipligi kullanin" if sev == "FAIL" else "")

    # --- PRB-03  C > G -------------------------------------------------------
    g, c = gc_counts(p)
    sev = "INFO" if c > g else "WARN"
    rec("PRB-03",
        {"g_count": g, "c_count": c,
         "g_count_other_strand": c, "c_count_other_strand": g},
        {"warn": "C sayisi G sayisindan fazla degil"},
        sev,
        {"rule": "C > G olan iplik secilir"},
        "tumleyen iplikte C=%d G=%d" % (g, c) if sev == "WARN" else "")

    # --- PRB-04  amplikon icinde ve primerlerden ayri -------------------------
    if amplicon is None:
        rec("PRB-04", {"evaluated": False}, {"fail": "primer ile cakisma"},
            "INFO", {"requires": "amplikon dizisi"},
            "amplikon verilmedi; konum denetimi yapilamadi")
    else:
        a = nn.validate(amplicon)
        pos = a.find(p)
        strand = "+"
        if pos < 0:
            pos = a.find(nn.reverse_complement(p))
            strand = "-"
        if pos < 0:
            rec("PRB-04", {"found_in_amplicon": False},
                {"fail": "prob amplikon icinde degil"}, "FAIL",
                {"search": "her iki iplik"},
                "prob dizisi amplikonda bulunamadi")
        else:
            end = pos + len(p)
            f_end = len(f)                 # ileri primerin 3' ucu (dislayici)
            r_start = len(a) - len(r)      # geri primerin baslangici
            overlap_f = max(0, f_end - pos)
            overlap_r = max(0, end - r_start)
            sev = "FAIL" if (overlap_f > 0 or overlap_r > 0) else "INFO"
            rec("PRB-04",
                {"found_in_amplicon": True, "strand": strand,
                 "start": pos, "end": end, "amplicon_bp": len(a),
                 "gap_to_forward_3p_nt": pos - f_end,
                 "gap_to_reverse_5p_nt": r_start - end,
                 "overlap_forward_nt": overlap_f,
                 "overlap_reverse_nt": overlap_r},
                {"fail": "primer baglanma yeriyle cakisma > 0 nt"},
                sev,
                {"coordinates": "amplikon icinde 0-tabanli"},
                "prob ile primer ayni yeri paylasiyor" if sev == "FAIL" else "")

    # --- PRB-05  prob dimerleri ----------------------------------------------
    partners = {"self": p, "forward": f, "reverse": r}
    worst_any = None
    for name, other in partners.items():
        for mode, w_thr, f_thr in (
            ("any", cons.dimer_dg_warn, cons.dimer_dg_fail),
            ("3prime", cons.dimer3p_dg_warn, cons.dimer3p_dg_fail),
        ):
            d = align.best_dimer(p, other, mode)
            cid = "PRB-05"
            tgt = f"{probe_id}:{name}:{mode}"
            if not d.found:
                checks.append(CheckRecord(
                    id=cid, target=tgt, value={"found": False},
                    threshold={"warn": f"dG37 < {w_thr}"}, severity="INFO",
                    model=align.model_metadata(),
                    conditions=cond.as_dict(), param_set_sha256=digest,
                    profile=profile_name, note="puanlanabilir dimer yok"))
                continue
            sev = ("FAIL" if d.dG37_kcal < f_thr
                   else "WARN" if d.dG37_kcal < w_thr else "INFO")
            if mode == "any" and (worst_any is None
                                  or d.dG37_kcal < worst_any[1]):
                worst_any = (name, d.dG37_kcal)
            checks.append(CheckRecord(
                id=cid, target=tgt, value=d.as_dict(cond.primer_M),
                threshold={"warn": f"dG37 < {w_thr}", "fail": f"dG37 < {f_thr}"},
                severity=sev,
                model={**align.model_metadata(),
                       "risk": "probun 3' ucu bloke; risk uzama degil "
                               "hibridizasyon rekabeti"},
                conditions=cond.as_dict(), param_set_sha256=digest,
                profile=profile_name))

    worst = max((ch.severity for ch in checks), key=lambda x: SEVERITY_ORDER[x])
    return {
        "probe_id": probe_id,
        "sequence": p,
        "length_nt": len(p),
        "gc_fraction": round(nn.gc_fraction(p), 4),
        "tm_c": round(m_p.tm_c, 3),
        "delta_tm_c": round(delta, 3),
        "worst_severity": worst,
        "worst_dimer": ({"partner": worst_any[0],
                         "dG37_kcal_mol": round(worst_any[1], 3)}
                        if worst_any else None),
        "checks": [ch.as_dict() for ch in checks],
        "engine_version": ENGINE_VERSION,
        "param_set_sha256": digest,
        "calibration_status": {"thresholds_calibrated": False},
    }


# ---------------------------------------------------------------------------
# Prob secimi
# ---------------------------------------------------------------------------
MAX_DIMER_SCAN = 120


def enumerate_probes(
    amplicon: str,
    forward: str,
    reverse: str,
    cond: salt.Conditions | None = None,
    cons: ProbeConstraints | None = None,
    salt_model: str = "owczarzy_2008",
) -> tuple[list[ProbeCandidate], list[tuple[str, str]]]:
    """Amplikonun primer disi bolgesinde iki iplikte aday prob uretir.

    Sert filtreler (uzunluk, GC, homopolimer, poli-G, 5' G, C>G) tasarim
    asamasinda uygulanir; gecen adaylar PRB-01 penceresine yakinliga gore
    puanlanir. Her elenme gerekcesiyle sayilir.
    """
    cond = cond or salt.Conditions()
    cons = cons or ProbeConstraints()
    a = nn.validate(amplicon)
    f = nn.validate(forward)
    r = nn.validate(reverse)

    inner_lo, inner_hi = len(f), len(a) - len(r)
    if inner_hi - inner_lo < cons.length_min:
        return [], [("PRB-04", "primerler arasi bolge proba yetmiyor")]

    tm_primer_max = max(
        concentration.build(f, cond=cond, salt_model=salt_model).tm_c,
        concentration.build(r, cond=cond, salt_model=salt_model).tm_c,
    )
    tm_mid = (cons.tm_delta_min + cons.tm_delta_max) / 2.0

    accepted: list[ProbeCandidate] = []
    rejections: list[tuple[str, str]] = []
    for start in range(inner_lo, inner_hi):
        for length in range(cons.length_min, cons.length_max + 1):
            end = start + length
            if end > inner_hi:
                break
            window = a[start:end]
            for strand, seq in (("+", window),
                                ("-", nn.reverse_complement(window))):
                bad = composition_reject(seq, cons)
                if bad:
                    rejections.append(bad)
                    continue
                m = concentration.build(seq, cond=cond, salt_model=salt_model)
                delta = m.tm_c - tm_primer_max
                if delta < cons.tm_delta_fail_lo or delta > cons.tm_delta_fail_hi:
                    rejections.append(("PRB-01", f"dTm {delta:.1f} C"))
                    continue
                g, c = gc_counts(seq)
                accepted.append(ProbeCandidate(
                    sequence=seq, strand=strand, start=start, end=end,
                    tm_c=m.tm_c, gc=nn.gc_fraction(seq), delta_tm_c=delta,
                    g_count=g, c_count=c,
                    score=abs(delta - tm_mid),
                ))

    # Dimer denetimi pahali: yalniz PRB-01 penceresine en yakin adaylara.
    accepted.sort(key=lambda p: p.score)
    scored: list[ProbeCandidate] = []
    for cand in accepted[:MAX_DIMER_SCAN]:
        worst = 0.0
        for other in (cand.sequence, f, r):
            d = align.best_dimer(cand.sequence, other, "any")
            if d.found:
                worst = min(worst, d.dG37_kcal)
        if worst < cons.dimer_dg_fail:
            rejections.append(("PRB-05", f"dimer {worst:.1f} kcal/mol"))
            continue
        # ceza: pencere disi Tm + esige yaklasan dimer
        penalty = cand.score + max(0.0, cons.dimer_dg_warn - worst)
        scored.append(ProbeCandidate(
            sequence=cand.sequence, strand=cand.strand, start=cand.start,
            end=cand.end, tm_c=cand.tm_c, gc=cand.gc,
            delta_tm_c=cand.delta_tm_c, g_count=cand.g_count,
            c_count=cand.c_count, score=penalty,
        ))
    for cand in accepted[MAX_DIMER_SCAN:]:
        rejections.append(("LIMIT", f"dimer taramasi siniri {MAX_DIMER_SCAN}"))

    scored.sort(key=lambda p: (p.score, p.start))
    return scored, rejections


def pick_probe(amplicon, forward, reverse, cond=None, cons=None,
               salt_model="owczarzy_2008") -> ProbeCandidate | None:
    cands, _ = enumerate_probes(amplicon, forward, reverse, cond, cons,
                                salt_model)
    return cands[0] if cands else None
