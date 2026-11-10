"""QC motoru: termodinamik kontroller ve kanit kayitlari.

Her kontrol, QC spesifikasyonundaki REP-02 kanit kaydi semasini uretir:
deger, esik, ciddiyet, model, kosullar, parametre seti checksum'i, motor
surumu, profil, override. Agent katmani bu kayitlarin disinda hicbir sayi
uretemez.

Tek primer kontrolleri termodinamik ve ikincil yapi kayitlari uretir.
Genomik baglanma yerlerine veya varyantlara bagli kontroller kendi QC
asamalarinda degerlendirilir.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from . import align, concentration, fold, nn, params, salt

ENGINE_VERSION = "oligoscout-thermo-ref/0.3.0"

SEVERITY_ORDER = {"INFO": 0, "WARN": 1, "FAIL": 2}


@dataclass(frozen=True)
class Profile:
    """Uygulama profili; her esik profil tarafindan ezilebilir."""

    name: str = "qpcr_sybr"
    tm_window_c: tuple[float, float] = (58.0, 62.0)
    tm_hard_margin_c: float = 3.0
    salt_spread_warn_c: float = 2.0
    salt_spread_fail_c: float = 4.0
    theta_warn: float = 0.95
    theta_fail: float = 0.80
    dg3p_warn_range: tuple[float, float] = (-9.0, -6.0)  # kcal/mol, disi WARN
    dg3p_fail: float = -11.0
    terminal_window_nt: int = 5
    dimer_dg_warn: float = -6.0          # DIM-01 / DIM-03
    dimer_dg_fail: float = -9.0
    dimer3p_dg_warn: float = -5.0        # DIM-02 / DIM-04
    dimer3p_dg_fail: float = -7.0
    delta_tm_warn: float = 2.0           # PAIR-01
    delta_tm_fail: float = 5.0

    def as_dict(self) -> dict[str, Any]:
        return {
            "profile": self.name,
            "tm_window_c": list(self.tm_window_c),
            "theta_warn": self.theta_warn,
            "theta_fail": self.theta_fail,
        }


@dataclass
class CheckRecord:
    """REP-02 kanit kaydi."""

    id: str
    target: str
    value: dict[str, Any]
    threshold: dict[str, Any]
    severity: str
    model: dict[str, Any]
    conditions: dict[str, Any]
    param_set_sha256: str
    engine_version: str = ENGINE_VERSION
    profile: str = ""
    override: Any = None
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "target": self.target,
            "value": self.value,
            "threshold": self.threshold,
            "severity": self.severity,
            "model": self.model,
            "conditions": self.conditions,
            "param_set_sha256": self.param_set_sha256,
            "engine_version": self.engine_version,
            "profile": self.profile,
            "override": self.override,
            "note": self.note,
        }


def _terminal_dG37(seq: str, window: int) -> float:
    """3' uctaki `window` bazin NN yigin dG37 toplami (baslatma terimi yok).

    TH-04 icin kullanilir. Baslatma teriminin disarida birakildigi acikca
    beyan edilir; literaturdeki '3' end stability' tanimi budur.
    """
    s = nn.validate(seq)[-window:]
    table = params.WC_SANTALUCIA_1998.require()["nn"]
    dH = dS = 0.0
    for i in range(len(s) - 1):
        step = s[i : i + 2]
        h, sv = table.get(step) or table[nn.reverse_complement(step)]
        dH += h
        dS += sv
    return dH - 310.15 * dS / 1000.0


class ThermoEngine:
    def __init__(
        self,
        cond: salt.Conditions | None = None,
        profile: Profile | None = None,
        salt_model: str = "owczarzy_2008",
    ) -> None:
        self.cond = cond or salt.Conditions()
        self.profile = profile or Profile()
        self.salt_model = salt_model
        self.param_digest = params.registry_digest()

    def evaluate_ct04(self, target: str, primer: str, matched_plus: str,
                      mismatched_plus: str, strand: str) -> CheckRecord:
        """Compare allele binding using 5′→3′ genomic plus-strand windows."""
        if self.profile.name != "ctdna" or strand not in {"+", "-"}:
            raise ValueError("CT-04 ctDNA profili ve + veya - iplik bilgisi gerektirir")
        p, match, mismatch = map(nn.validate, (primer, matched_plus, mismatched_plus))
        if strand == "-":
            match, mismatch = map(nn.reverse_complement, (match, mismatch))
        if p != match:
            raise ValueError("CT-04: eşleşen alel, primer yönünde tam eşleşmelidir")
        on = concentration.build(p, self.cond, "primer_excess", self.salt_model)
        off = concentration.build(p, self.cond, "primer_excess", self.salt_model,
                                  template_site=mismatch, arms_mismatch=True)
        position = nn.mismatch_positions(p, mismatch)[0]
        ta = self.cond.anneal_c
        theta_on, theta_off = on.theta(ta), off.theta(ta)
        margin = theta_on - theta_off
        delta_dg = off.duplex.dG_kcal() - on.duplex.dG_kcal()
        # Literature motivates positive allele discrimination. These numeric
        # starting points are project heuristics and require experimental data.
        warn_theta, warn_dg = 0.20, 0.50
        severity = "WARN" if margin < warn_theta or delta_dg < warn_dg else "INFO"
        terminal = position == len(p) - 1
        table = (params.TERMINAL_MISMATCH_SANTALUCIA_PEYRET if terminal
                 else params.MISMATCH_ALLAWI_PEYRET).provenance
        return CheckRecord(
            id="CT-04", target=target, severity=severity, profile=self.profile.name,
            conditions=self.cond.as_dict(), param_set_sha256=self.param_digest,
            value={"strand": strand, "position_from_3prime": len(p) - position,
                   "mismatch_kind": "terminal" if terminal else "penultimate",
                   "theta_matched": theta_on, "theta_mismatched": theta_off, "theta_margin": margin,
                   "dG37_matched_kcal_mol": on.duplex.dG_kcal(),
                   "dG37_mismatched_kcal_mol": off.duplex.dG_kcal(), "delta_delta_G37_kcal_mol": delta_dg},
            threshold={"theta_margin_warn_below": warn_theta, "delta_delta_G37_warn_below": warn_dg,
                       "status": "REP-03 / KALIBRE DEGIL"},
            model={"mismatch_nn": table.table_id, "verification_status": table.verification_status,
                   "salt_model": self.salt_model, "conc_model": "primer_excess",
                   "reference_conditions": "37 C, 1 M NaCl (dG37)"},
            note="Newton 1989 ve Huang 1992 temelli karar desteği. Sayısal eşikler proje başlangıç değerleridir. "
                 "Polimeraz uzaması bu kayıtta modellenmez.",
        )

    # -- tek primer ---------------------------------------------------------
    def evaluate(self, primer_id: str, seq: str) -> dict[str, Any]:
        s = nn.validate(seq)
        model = concentration.build(
            s, cond=self.cond, salt_model=self.salt_model
        )
        ta = self.cond.anneal_c
        checks: list[CheckRecord] = []

        def rec(
            cid: str, value: dict, threshold: dict, severity: str, extra_model: dict,
            note: str = "", *, conditions: dict | None = None,
        ) -> None:
            checks.append(
                CheckRecord(
                    id=cid,
                    target=primer_id,
                    value=value,
                    threshold=threshold,
                    severity=severity,
                    model={
                        "nn_table": params.WC_SANTALUCIA_1998.provenance.table_id,
                        "nn_table_status":
                            params.WC_SANTALUCIA_1998.provenance.verification_status,
                        **extra_model,
                    },
                    conditions=self.cond.as_dict() if conditions is None else conditions,
                    param_set_sha256=self.param_digest,
                    profile=self.profile.name,
                    note=note,
                )
            )

        # --- TH-01  Tm (NN) ------------------------------------------------
        lo, hi = self.profile.tm_window_c
        m = self.profile.tm_hard_margin_c
        tm = model.tm_c
        sev = "INFO"
        if tm < lo - m or tm > hi + m:
            sev = "FAIL"
        elif tm < lo or tm > hi:
            sev = "WARN"
        rec(
            "TH-01",
            {"tm_c": round(tm, 3), "tm_1M_c": round(model.tm_1m_K - params.KELVIN, 3),
             "dH_kcal_mol": round(model.duplex.dH_kcal, 3),
             "dS_cal_mol_K": round(model.duplex.dS_cal, 3)},
            {"warn": f"{lo}-{hi} C disi", "fail": f"{lo - m}-{hi + m} C disi"},
            sev,
            {"conc_model": model.conc_model, "salt_model": model.salt_model,
             "salt_regime": model.salt_detail.get("regime")},
        )

        # --- TH-02  tuz modeli duyarliligi ---------------------------------
        sp = salt.spread(model.tm_1m_K, model.f_gc, model.duplex.n_bp, self.cond)
        spread_c = sp["spread_c"]
        sev = ("FAIL" if spread_c > self.profile.salt_spread_fail_c
               else "WARN" if spread_c > self.profile.salt_spread_warn_c else "INFO")
        rec(
            "TH-02",
            {k: round(v, 3) for k, v in sp.items()},
            {"warn": f">{self.profile.salt_spread_warn_c} C",
             "fail": f">{self.profile.salt_spread_fail_c} C"},
            sev,
            {"models_reported": list(salt.MODELS),
             "models_compared": list(salt.applicable_models(self.cond))},
            "model belirsizliginin raporlanmasi (REP-03); yayilim yalniz "
            "beyan edilen kosullara uygulanabilir modeller uzerinden",
        )

        # --- TH-03  hedefe bagli fraksiyon theta(Ta) ------------------------
        theta = model.theta(ta)
        sev = ("FAIL" if theta < self.profile.theta_fail
               else "WARN" if theta < self.profile.theta_warn else "INFO")
        rec(
            "TH-03",
            {"theta_at_ta": round(theta, 6), "ta_c": ta,
             "dG_ta_kcal_mol": round(model.dG_kcal(ta), 4)},
            {"warn": f"<{self.profile.theta_warn}", "fail": f"<{self.profile.theta_fail}"},
            sev,
            {"conc_model": model.conc_model},
        )

        # --- TH-04  3' uc stabilitesi --------------------------------------
        w = self.profile.terminal_window_nt
        dg3 = _terminal_dG37(s, w)
        lo_g, hi_g = self.profile.dg3p_warn_range
        if dg3 < self.profile.dg3p_fail:
            sev, side = "FAIL", "asiri_kararli"
        elif dg3 < lo_g:
            sev, side = "WARN", "fazla_kararli"
        elif dg3 > hi_g:
            sev, side = "WARN", "zayif"
        else:
            sev, side = "INFO", "pencere_ici"
        rec(
            "TH-04",
            {"dG37_last_%d_nt_kcal_mol" % w: round(dg3, 4), "side": side},
            {"warn": f">{hi_g} veya <{lo_g} kcal/mol",
             "fail": f"<{self.profile.dg3p_fail} kcal/mol"},
            sev,
            {"definition": "son %d bazin NN yigin toplami, baslatma terimi haric" % w},
            "esik literatur kilavuz degeridir; islak-lab kalibrasyonu bekler",
        )

        # --- TH-05  erime egrisi ic tutarliligi -----------------------------
        theta_at_tm = model.theta(model.tm_c)
        err = abs(theta_at_tm - 0.5)
        sev = "FAIL" if err > 1e-9 else "INFO"
        rec(
            "TH-05",
            {"theta_at_tm": round(theta_at_tm, 12), "abs_error": err},
            {"fail": "|theta(Tm) - 0.5| > 1e-9"},
            sev,
            {"identity": "theta(Tm) = 0.5"},
            "analitik ozdeslik testi",
        )

        # --- DIM-01 / DIM-02  self-dimer -----------------------------------
        for cid, mode, w_thr, f_thr in (
            ("DIM-01", "any", self.profile.dimer_dg_warn, self.profile.dimer_dg_fail),
            ("DIM-02", "3prime", self.profile.dimer3p_dg_warn,
             self.profile.dimer3p_dg_fail),
        ):
            d = align.self_dimer(s, mode)
            if not d.found:
                rec(cid, {"found": False}, {"warn": f"dG37 < {w_thr}"}, "INFO",
                    align.model_metadata(), "puanlanabilir dimer yok")
                continue
            sev = ("FAIL" if d.dG37_kcal < f_thr
                   else "WARN" if d.dG37_kcal < w_thr else "INFO")
            rec(cid, d.as_dict(self.cond.primer_M),
                {"warn": f"dG37 < {w_thr}", "fail": f"dG37 < {f_thr}"}, sev,
                align.model_metadata())

        # --- SS-01 / SS-03: DNA ikincil yapi, referans kosullarinda ---------
        loop_table = params.DNA_LOOP_SANTALUCIA_HICKS_2004.provenance
        fold_model = {
            "model": fold.DnaLoopEnergies().name,
            "loop_table": loop_table.table_id,
            "loop_table_status": loop_table.verification_status,
            "scope": "WC yiginlari, uzunluga dayali donguler ve secili kapanis terimleri",
            "excluded": ["cok dalli donguler", "1x1 ic donguler",
                         "diger uc terimleri ve kisa sac tokasi bonuslari"],
        }
        fold_conditions = {"temperature_c": 37.0, "monovalent_mM": 1000.0}
        mfe, structure = fold.zuker_mfe(s)
        rec("SS-01", {"mfe_kcal_mol": round(mfe, 4), "dot_bracket": structure},
            {}, "INFO", fold_model, "kalibre karar esigi yok", conditions=fold_conditions)
        unpaired = fold.unpaired_probability(s)
        terminal = unpaired[-min(len(s), max(1, self.profile.terminal_window_nt)):]
        rec("SS-03", {"p_unpaired_by_base": [round(p, 6) for p in unpaired],
                      "mean_p_unpaired_3p": round(sum(terminal) / len(terminal), 6),
                      "window_3p_nt": len(terminal)},
            {}, "INFO", fold_model, "her baz icin marjinal eslesmeme olasiligi",
            conditions=fold_conditions)

        # --- degerlendirilemeyenler ----------------------------------------
        not_eval = {
            "SPC-04": "genomik baglanma kayitlari panel QC asamasinda kullanilir",
            "VAR-04": "referans ve alternatif aleller varyant QC asamasinda kullanilir",
            "TH-06": "katki maddesi katsayilari aktarilmadi",
            "TH-07": "urun Tm'i gerekir (amplikon modulu)",
        }

        worst = max((c.severity for c in checks), key=lambda x: SEVERITY_ORDER[x])
        return {
            "primer_id": primer_id,
            "sequence": s,
            "length_nt": len(s),
            "gc_fraction": round(model.f_gc, 4),
            "thermo": model.as_dict(),
            "worst_severity": worst,
            "checks": [c.as_dict() for c in checks],
            "not_evaluated": not_eval,
            "engine_version": ENGINE_VERSION,
            "param_set_sha256": self.param_digest,
            "calibration_status": {
                "ranking_score_calibrated": False,
                "parameter_tables_verified": all(
                    t.provenance.verification_status == "verified"
                    for t in params.REGISTRY.values()
                ),
            },
        }

    # -- primer cifti -------------------------------------------------------
    def evaluate_pair(self, pair_id: str, fwd: str, rev: str) -> dict[str, Any]:
        """DIM-03, DIM-04 (hetero-dimer) ve PAIR-01 (dTm)."""
        f, r = nn.validate(fwd), nn.validate(rev)
        mf = concentration.build(f, cond=self.cond, salt_model=self.salt_model)
        mr = concentration.build(r, cond=self.cond, salt_model=self.salt_model)
        checks: list[CheckRecord] = []

        def rec(cid, value, threshold, severity, model, note=""):
            checks.append(CheckRecord(
                id=cid, target=pair_id, value=value, threshold=threshold,
                severity=severity, model=model, conditions=self.cond.as_dict(),
                param_set_sha256=self.param_digest, profile=self.profile.name,
                note=note))

        for cid, mode, w_thr, f_thr in (
            ("DIM-03", "any", self.profile.dimer_dg_warn, self.profile.dimer_dg_fail),
            ("DIM-04", "3prime", self.profile.dimer3p_dg_warn,
             self.profile.dimer3p_dg_fail),
        ):
            d = align.best_dimer(f, r, mode)
            if not d.found:
                rec(cid, {"found": False}, {"warn": f"dG37 < {w_thr}"}, "INFO",
                    align.model_metadata(), "puanlanabilir dimer yok")
                continue
            sev = ("FAIL" if d.dG37_kcal < f_thr
                   else "WARN" if d.dG37_kcal < w_thr else "INFO")
            rec(cid, d.as_dict(self.cond.primer_M),
                {"warn": f"dG37 < {w_thr}", "fail": f"dG37 < {f_thr}"}, sev,
                align.model_metadata())

        dtm = abs(mf.tm_c - mr.tm_c)
        sev = ("FAIL" if dtm > self.profile.delta_tm_fail
               else "WARN" if dtm > self.profile.delta_tm_warn else "INFO")
        rec("PAIR-01",
            {"delta_tm_c": round(dtm, 3), "tm_f_c": round(mf.tm_c, 3),
             "tm_r_c": round(mr.tm_c, 3)},
            {"warn": f">{self.profile.delta_tm_warn} C",
             "fail": f">{self.profile.delta_tm_fail} C"},
            sev, {"salt_model": self.salt_model})

        worst = max((c.severity for c in checks), key=lambda x: SEVERITY_ORDER[x])
        return {
            "pair_id": pair_id, "forward": f, "reverse": r,
            "worst_severity": worst,
            "checks": [c.as_dict() for c in checks],
            "engine_version": ENGINE_VERSION,
            "param_set_sha256": self.param_digest,
        }

    def melting_curve(self, seq: str, t_min: float = 20.0, t_max: float = 95.0,
                      step: float = 0.25) -> list[tuple[float, float]]:
        model = concentration.build(nn.validate(seq), cond=self.cond,
                                    salt_model=self.salt_model)
        return model.melting_curve(t_min, t_max, step)
