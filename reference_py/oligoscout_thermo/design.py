"""Primer tasarim motoru: aday uretimi, sert filtreler, cift kurgulama, siralama.

QC spesifikasyonu bolum 2, madde 6'nin uygulamasi:
    aday uretimi -> sert filtreler (FAIL) -> cok amacli skor (Pareto)
    -> [OligoScout ozgullugu] -> [varyant/haplotip] -> yeniden siralama
Burada ilk uc asama var. Ozgulluk taramasi ayri bir adimdir: tasarlanan
ciftler `oligoscout_panel` araciyla taranir ve `panel_qc.py` ile yeniden
siralanir.

Tasarim ilkesi: **her elenme gerekcesiyle kaydedilir.** Hicbir aday sessizce
dusmez; `rejections` listesi hangi kuralin hangi degerle tetiklendigini
tasir ve QC katalogundaki madde kimligiyle etiketlenir.

Skor kalibre DEGILDIR ve oyle etiketlenir (REP-03).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from . import align, concentration, nn, salt

MAX_DIMER_SCAN = 4000      # dimer denetimine giren en fazla aday


@dataclass(frozen=True)
class DesignConstraints:
    """Sert filtreler ve pencereler. Hepsi profil tarafindan ezilebilir."""

    length_min: int = 18
    length_max: int = 27
    tm_min: float = 58.0
    tm_max: float = 62.0
    tm_opt: float = 60.0
    gc_min: float = 0.40
    gc_max: float = 0.60
    max_homopolymer: int = 4
    max_repeat_units: int = 4
    gc_clamp_min: int = 1        # son 5 nt icinde G/C sayisi
    gc_clamp_max: int = 3
    self_dimer_dg: float = -6.0
    self_dimer3p_dg: float = -5.0
    hetero_dimer_dg: float = -6.0
    hetero_dimer3p_dg: float = -5.0
    product_min: int = 70
    product_max: int = 200
    max_delta_tm: float = 2.0
    max_pairs: int = 2000
    theta_min: float = 0.0


@dataclass(frozen=True)
class Candidate:
    """Bir aday primer."""

    sequence: str
    start: int                  # hedef dizide 0-tabanli
    end: int                    # dislayici
    orientation: str            # "F" | "R"
    tm_c: float
    gc: float
    theta_ta: float
    self_dimer_dg: float
    self_dimer3p_dg: float

    def as_row(self) -> list:
        return [self.orientation, self.sequence, self.start, self.end,
                len(self.sequence), round(100 * self.gc, 2), round(self.tm_c, 3),
                round(self.theta_ta, 6), round(self.self_dimer_dg, 3),
                round(self.self_dimer3p_dg, 3)]


@dataclass(frozen=True)
class Rejection:
    sequence: str
    start: int
    orientation: str
    rule: str                   # QC katalog kimligi
    reason: str
    value: str

    def as_row(self) -> list:
        return [self.orientation, self.sequence, self.start, self.rule,
                self.reason, self.value]


@dataclass
class PairCandidate:
    left: Candidate
    right: Candidate
    product_start: int
    product_end: int
    product_length: int
    delta_tm: float
    hetero_dg: float
    hetero3p_dg: float
    objectives: tuple[float, ...] = field(default_factory=tuple)
    score: float = 0.0
    pareto_rank: int = 0

    def as_row(self) -> list:
        return [self.left.sequence, self.right.sequence,
                self.product_start, self.product_end, self.product_length,
                round(self.left.tm_c, 3), round(self.right.tm_c, 3),
                round(self.delta_tm, 3), round(self.hetero_dg, 3),
                round(self.hetero3p_dg, 3), round(self.score, 4),
                self.pareto_rank]


# ---------------------------------------------------------------------------
# Kompozisyon denetimleri
# ---------------------------------------------------------------------------
def longest_homopolymer(seq: str) -> int:
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best if seq else 0


def max_repeat_units(seq: str, unit: int) -> int:
    """Ardisik tekrar eden `unit` uzunluklu birimin en fazla tekrar sayisi."""
    best = 1
    for i in range(len(seq) - unit):
        motif = seq[i : i + unit]
        count = 1
        j = i + unit
        while seq[j : j + unit] == motif:
            count += 1
            j += unit
        best = max(best, count)
    return best


def gc_clamp(seq: str, window: int = 5) -> int:
    tail = seq[-window:]
    return tail.count("G") + tail.count("C")


# ---------------------------------------------------------------------------
# Aday uretimi
# ---------------------------------------------------------------------------
def enumerate_candidates(
    target: str,
    orientation: str,
    cons: DesignConstraints,
    cond: salt.Conditions,
    region: tuple[int, int] | None = None,
) -> tuple[list[Candidate], list[Rejection]]:
    """Verilen bolgedeki tum uzunluklarda aday uretir ve sert filtreleri uygular."""
    seq = nn.validate(target)
    lo, hi = region or (0, len(seq))
    lo = max(0, lo)
    hi = min(len(seq), hi)

    accepted: list[Candidate] = []
    rejected: list[Rejection] = []

    for start in range(lo, hi):
        for length in range(cons.length_min, cons.length_max + 1):
            end = start + length
            if end > hi:
                break
            window = seq[start:end]
            primer = window if orientation == "F" else nn.reverse_complement(window)

            def reject(rule: str, reason: str, value) -> None:
                rejected.append(
                    Rejection(primer, start, orientation, rule, reason, str(value))
                )

            gc = nn.gc_fraction(primer)
            if not (cons.gc_min <= gc <= cons.gc_max):
                reject("CMP-01", "GC penceresi disi", f"{100 * gc:.1f}%")
                continue

            clamp = gc_clamp(primer)
            if not (cons.gc_clamp_min <= clamp <= cons.gc_clamp_max):
                reject("CMP-02", "3' GC clamp", clamp)
                continue

            homo = longest_homopolymer(primer)
            if homo >= cons.max_homopolymer:
                reject("CMP-03", "homopolimer", homo)
                continue

            if "GGGG" in primer:
                reject("CMP-03", "poli-G", "GGGG")
                continue

            rep2 = max_repeat_units(primer, 2)
            rep3 = max_repeat_units(primer, 3)
            if max(rep2, rep3) >= cons.max_repeat_units:
                reject("CMP-04", "di/tri-nukleotit tekrar", max(rep2, rep3))
                continue

            model = concentration.build(primer, cond=cond)
            tm = model.tm_c
            if not (cons.tm_min <= tm <= cons.tm_max):
                reject("TH-01", "Tm penceresi disi", f"{tm:.2f} C")
                continue

            theta = model.theta(cond.anneal_c)
            if theta < cons.theta_min:
                reject("TH-03", "hedefe baglanma kesri esik altinda",
                       f"{theta:.6f} < {cons.theta_min:g}")
                continue

            accepted.append(
                Candidate(
                    sequence=primer, start=start, end=end,
                    orientation=orientation, tm_c=tm, gc=gc,
                    theta_ta=theta,
                    self_dimer_dg=0.0, self_dimer3p_dg=0.0,
                )
            )

    # Dimer denetimi pahali: yalniz hayatta kalanlara, Tm'e yakinliga gore
    accepted.sort(key=lambda c: abs(c.tm_c - cons.tm_opt))
    scanned: list[Candidate] = []
    for cand in accepted[:MAX_DIMER_SCAN]:
        d_any = align.self_dimer(cand.sequence, "any")
        d_3p = align.self_dimer(cand.sequence, "3prime")
        g_any = d_any.dG37_kcal if d_any.found else 0.0
        g_3p = d_3p.dG37_kcal if d_3p.found else 0.0
        if g_any < cons.self_dimer_dg:
            rejected.append(Rejection(cand.sequence, cand.start, orientation,
                                      "DIM-01", "self-dimer", f"{g_any:.2f}"))
            continue
        if g_3p < cons.self_dimer3p_dg:
            rejected.append(Rejection(cand.sequence, cand.start, orientation,
                                      "DIM-02", "3'-uclu self-dimer",
                                      f"{g_3p:.2f}"))
            continue
        scanned.append(
            Candidate(
                sequence=cand.sequence, start=cand.start, end=cand.end,
                orientation=cand.orientation, tm_c=cand.tm_c, gc=cand.gc,
                theta_ta=cand.theta_ta, self_dimer_dg=g_any,
                self_dimer3p_dg=g_3p,
            )
        )

    for cand in accepted[MAX_DIMER_SCAN:]:
        rejected.append(Rejection(cand.sequence, cand.start, orientation,
                                  "LIMIT", "dimer taramasi siniri",
                                  str(MAX_DIMER_SCAN)))
    return scanned, rejected


# ---------------------------------------------------------------------------
# Cift kurgulama ve siralama
# ---------------------------------------------------------------------------
def build_pairs(
    lefts: list[Candidate],
    rights: list[Candidate],
    cons: DesignConstraints,
) -> tuple[list[PairCandidate], list[Rejection]]:
    pairs: list[PairCandidate] = []
    rejected: list[Rejection] = []

    rights_sorted = sorted(rights, key=lambda c: c.end)
    for left in lefts:
        for right in rights_sorted:
            if right.start < left.end:
                continue
            length = right.end - left.start
            if length < cons.product_min:
                continue
            if length > cons.product_max:
                break

            dtm = abs(left.tm_c - right.tm_c)
            if dtm > cons.max_delta_tm:
                rejected.append(Rejection(
                    f"{left.sequence}+{right.sequence}", left.start, "P",
                    "PAIR-01", "dTm", f"{dtm:.2f}"))
                continue

            d_any = align.best_dimer(left.sequence, right.sequence, "any")
            d_3p = align.best_dimer(left.sequence, right.sequence, "3prime")
            g_any = d_any.dG37_kcal if d_any.found else 0.0
            g_3p = d_3p.dG37_kcal if d_3p.found else 0.0
            if g_any < cons.hetero_dimer_dg:
                rejected.append(Rejection(
                    f"{left.sequence}+{right.sequence}", left.start, "P",
                    "DIM-03", "hetero-dimer", f"{g_any:.2f}"))
                continue
            if g_3p < cons.hetero_dimer3p_dg:
                rejected.append(Rejection(
                    f"{left.sequence}+{right.sequence}", left.start, "P",
                    "DIM-04", "3'-uclu hetero-dimer", f"{g_3p:.2f}"))
                continue

            pairs.append(PairCandidate(
                left=left, right=right,
                product_start=left.start, product_end=right.end,
                product_length=length, delta_tm=dtm,
                hetero_dg=g_any, hetero3p_dg=g_3p,
            ))
            if len(pairs) >= cons.max_pairs:
                return pairs, rejected
    return pairs, rejected


def _objectives(pair: PairCandidate, cons: DesignConstraints) -> tuple[float, ...]:
    """Hepsi 'kucuk daha iyi'. Kalibre edilmemis ceza terimleri."""
    tm_dev = (abs(pair.left.tm_c - cons.tm_opt)
              + abs(pair.right.tm_c - cons.tm_opt)) / 2.0
    gc_dev = (abs(pair.left.gc - 0.5) + abs(pair.right.gc - 0.5)) / 2.0
    dimer = -min(0.0, min(pair.hetero_dg, pair.left.self_dimer_dg,
                          pair.right.self_dimer_dg))
    dimer3p = -min(0.0, min(pair.hetero3p_dg, pair.left.self_dimer3p_dg,
                            pair.right.self_dimer3p_dg))
    theta_loss = (1.0 - pair.left.theta_ta) + (1.0 - pair.right.theta_ta)
    return (tm_dev, pair.delta_tm, gc_dev, dimer, dimer3p, theta_loss)


def rank_pairs(
    pairs: list[PairCandidate], cons: DesignConstraints
) -> list[PairCandidate]:
    """Pareto katmanlari + beyan edilmis, KALIBRE EDILMEMIS skalar skor."""
    for p in pairs:
        p.objectives = _objectives(p, cons)

    # skalar: her amacin kendi araligina normalize edilmis toplami
    if pairs:
        n_obj = len(pairs[0].objectives)
        lo = [min(p.objectives[i] for p in pairs) for i in range(n_obj)]
        hi = [max(p.objectives[i] for p in pairs) for i in range(n_obj)]
        for p in pairs:
            total = 0.0
            for i in range(n_obj):
                span = hi[i] - lo[i]
                total += 0.0 if span <= 0 else (p.objectives[i] - lo[i]) / span
            p.score = total / n_obj

    # Pareto katmanlari
    remaining = list(pairs)
    rank = 1
    while remaining:
        front = []
        for a in remaining:
            dominated = False
            for b in remaining:
                if b is a:
                    continue
                if (all(x <= y for x, y in zip(b.objectives, a.objectives))
                        and any(x < y for x, y in zip(b.objectives, a.objectives))):
                    dominated = True
                    break
            if not dominated:
                front.append(a)
        for p in front:
            p.pareto_rank = rank
        remaining = [p for p in remaining if p.pareto_rank == 0]
        rank += 1
        if rank > 50:                      # guvenlik siniri
            for p in remaining:
                p.pareto_rank = rank
            break

    pairs.sort(key=lambda p: (p.pareto_rank, p.score))
    return pairs


def design(
    target: str,
    cons: DesignConstraints | None = None,
    cond: salt.Conditions | None = None,
    left_region: tuple[int, int] | None = None,
    right_region: tuple[int, int] | None = None,
) -> dict:
    """Uctan uca tasarim. Ozgulluk taramasi bu adimda YAPILMAZ."""
    cons = cons or DesignConstraints()
    cond = cond or salt.Conditions()

    lefts, rej_l = enumerate_candidates(target, "F", cons, cond, left_region)
    rights, rej_r = enumerate_candidates(target, "R", cons, cond, right_region)
    pairs, rej_p = build_pairs(lefts, rights, cons)
    pairs = rank_pairs(pairs, cons)

    return {
        "target_length": len(nn.validate(target)),
        "left_candidates": lefts,
        "right_candidates": rights,
        "pairs": pairs,
        "rejections": rej_l + rej_r + rej_p,
        "calibration_status": {
            "ranking_score_calibrated": False,
            "specificity_screened": False,
        },
    }
