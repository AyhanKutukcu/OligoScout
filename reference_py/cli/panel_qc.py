#!/usr/bin/env python3
"""oligoscout_panel kayitlarini termodinamik QC motoruyla birlestirir.

Girdi : <prefix>.bindings.tsv, <prefix>.products.tsv, <prefix>.contract.tsv
Cikti : panel_qc.tsv          cift basina ozgulluk + termodinamik ozet
        offtarget_products.tsv  hedef disi urunler, en kotuden siralanmis
        panel_qc_records.jsonl  her kontrol icin REP-02 kanit kaydi

C++ baglanma yerlerini arar. Python katmani SPC-04 icin desteklenen tekli
ic uyumsuzluklari, kayittaki referans dizisini kullanarak puanlar.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import (annotations, concentration, engine, nn,  # noqa: E402
                               params, salt)

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def three_prime_anchor_intact(mask: str, anchor: int) -> bool:
    """3' uctaki `anchor` nt icinde uyumsuzluk var mi (maske 5'->3')."""
    return "X" not in mask[-anchor:] if mask else True


_TOLERANT_MM = frozenset({frozenset("GT")})   # G.T / T.G sallanma


def _classify_terminal_mismatch(binding: dict, tally, window: int = 3) -> None:
    """SPC-05: 3' son `window` nt'deki uyumsuzluklari tolerant/disruptive say.

    Uyumsuzluk cifti = (primer bazi, karsisindaki kaliptaki baz). Maske
    primer 5'->3' hizalidir ve primerle AYNI ANLAMLI iplige gore X isaretler
    (+ : site; - : revcomp(site)). Primerin 3' ucu `mask[-1]`dir.
    """
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    mask = binding["mask"]
    if not mask or "X" not in mask[-window:]:
        return
    primer = binding["primer"]
    site = binding["site_plus_strand"]
    same_sense = (site if binding["strand"] == "+"
                  else "".join(comp.get(b, "N") for b in reversed(site)))
    for i in range(len(mask) - window, len(mask)):
        if i < 0 or mask[i] != "X":
            continue
        opp = comp.get(same_sense[i], "N")       # primerin karsisindaki baz
        pair = frozenset({primer[i], opp})
        tally["tolerant" if pair in _TOLERANT_MM else "disruptive"] += 1


def score_mismatched_site(binding: dict, cond: salt.Conditions, salt_model: str) -> dict:
    if int(binding["mismatches"]) != 1:
        raise ValueError("yalniz tek uyumsuzluk desteklenir")
    target = binding.get("site_plus_strand", "")
    strand = binding.get("strand")
    if not target or strand not in ("+", "-"):
        raise ValueError("referans baglanma dizisi veya iplik bilgisi yok")
    if strand == "-":
        target = nn.reverse_complement(target)
    primer = binding["primer"]
    if len(nn.mismatch_positions(primer, target)) != 1:
        raise ValueError("referans dizisi arama kaydiyla uyusmuyor")
    model = concentration.build(primer, cond=cond, conc_model="primer_excess",
                                salt_model=salt_model, template_site=target)
    return {"chrom": binding["chrom"], "start": int(binding["start"]),
            "strand": strand, "theta_at_ta": model.theta(cond.anneal_c),
            "dG37_1M_kcal_mol": round(model.duplex.dG_kcal(), 4)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True,
                    help="oligoscout_panel cikti oneki")
    ap.add_argument("--out", default="out/panel")
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    ap.add_argument("--tm-window", default="58:62")
    ap.add_argument("--anchor", type=int, default=12)
    ap.add_argument("--panel", help="beklenen koordinatlari iceren panel TSV'si; "
                                    "verilirse hedef/hedef disi ayrimi yapilir")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    lo, hi = (float(x) for x in args.tm_window.split(":"))

    cond = salt.Conditions(
        monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
        dntp_M=args.dntp / 1e3, primer_M=args.primer_nm / 1e9,
        anneal_c=args.ta,
    )
    eng = engine.ThermoEngine(cond=cond, profile=engine.Profile(tm_window_c=(lo, hi)))

    bindings = read_tsv(args.prefix + ".bindings.tsv")
    products = read_tsv(args.prefix + ".products.tsv")
    contract = {r["field"]: r["value"]
                for r in read_tsv(args.prefix + ".contract.tsv")}

    # --- primer duzeyi termodinamik -----------------------------------------
    primers: dict[tuple[str, str], str] = {}
    for b in bindings:
        primers[(b["pair_id"], b["role"])] = b["primer"]

    thermo: dict[str, dict] = {}
    for (pid, role), seq in primers.items():
        key = f"{pid}_{role}"
        thermo[key] = eng.evaluate(key, seq)

    # --- ozgulluk sayimlari -------------------------------------------------
    sites = collections.defaultdict(lambda: collections.Counter())
    anchored = collections.defaultdict(int)
    exact = collections.defaultdict(int)
    # SPC-05: hedef disi yerlerde 3' terminal (son 3 nt) uyumsuzluk tipi.
    # tolerant (G.T / T.G) 3' uc uyumsuzlugu uzamaya izin verebilir -> risk.
    term3 = collections.defaultdict(lambda: collections.Counter())
    # SPC-06 / SPC-07: baglanma yerinin psodogen / segmental duplikasyon
    # ortusmesi. Anotasyon dosyalari yoksa indexler bostur ve kurallar
    # 'anotasyon yok' raporlar.
    ann_paths = annotations.default_paths(_REPO_ROOT)
    pseudo_idx = annotations.load_intervals(ann_paths["pseudogene"])
    segdup_idx = annotations.load_intervals(ann_paths["segdup"])
    pseudo_hits = collections.defaultdict(set)
    segdup_hits = collections.defaultdict(int)
    scored = collections.Counter()
    unscored = collections.defaultdict(collections.Counter)
    strongest: dict[tuple[str, str], dict] = {}
    for b in bindings:
        key = (b["pair_id"], b["role"])
        mm = int(b["mismatches"])
        sites[key][mm] += 1
        if pseudo_idx.total or segdup_idx.total:
            bs, be = int(b["start"]), int(b["end_exclusive"])
            for name in pseudo_idx.overlaps(b["chrom"], bs, be):
                pseudo_hits[key].add(name)
            if segdup_idx.overlaps_any(b["chrom"], bs, be):
                segdup_hits[key] += 1
        if mm == 0:
            exact[key] += 1
            continue
        if three_prime_anchor_intact(b["mask"], args.anchor):
            anchored[key] += 1
        _classify_terminal_mismatch(b, term3[key])
        try:
            result = score_mismatched_site(b, cond, eng.salt_model)
        except (ValueError, params.ParameterTableUnavailable) as error:
            unscored[key][str(error)] += 1
            continue
        scored[key] += 1
        if key not in strongest or result["theta_at_ta"] > strongest[key]["theta_at_ta"]:
            strongest[key] = result

    mismatch_table = params.MISMATCH_ALLAWI_PEYRET.provenance
    for (pid, role) in primers:
        key = (pid, role)
        rec = thermo[f"{pid}_{role}"]
        best = strongest.get(key)
        if best:
            best = {**best, "theta_at_ta": round(best["theta_at_ta"], 6)}
            rec["not_evaluated"].pop("SPC-04", None)
        else:
            rec["not_evaluated"]["SPC-04"] = ("desteklenen uyumsuzluk yok"
                                             if unscored[key] else "uyumsuzluklu yer yok")
        rec["checks"].append(engine.CheckRecord(
            id="SPC-04", target=f"{pid}_{role}",
            value={"evaluated": bool(best), "n_scored": scored[key],
                   "n_unscored": sum(unscored[key].values()),
                   "unscored_reasons": dict(unscored[key]), "strongest_site": best},
            threshold={}, severity="INFO",
            model={"nn_table": mismatch_table.table_id,
                   "nn_table_status": mismatch_table.verification_status,
                   "salt_model": eng.salt_model,
                   "salt_correction": "primer GC oranina dayali yaklasim"},
            conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
            profile=eng.profile.name,
            note="tekli ic uyumsuzluk icin theta(Ta); kalibre karar esigi yok",
        ).as_dict())

        # SPC-01: hedef tekilligi — tam eslesmeli genom sitesi sayisi.
        n_exact = exact[key]
        sev = "INFO"
        if n_exact > 1:
            sev = "FAIL"          # birden fazla tam eslesme -> ozgul degil
        elif n_exact == 0:
            sev = "WARN"          # tam eslesme yok (referans/tasarim uyumu?)
        rec["checks"].append(engine.CheckRecord(
            id="SPC-01", target=f"{pid}_{role}",
            value={"exact_match_sites": n_exact},
            threshold={"warn": "tam eslesmeli site yok",
                       "fail": "birden fazla tam eslesmeli site (ozgul degil)"},
            severity=sev, model={"source": "genome_scan_bindings"},
            conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
            profile=eng.profile.name,
            note="primer duzeyi genom tekilligi (tam eslesme sayimi)",
        ).as_dict())

        # SPC-02: hedef disi baglanma yerleri, uyumsuzluk sinifi x 3'-capa.
        hist = sites[key]
        mm_hist = {str(k): int(v) for k, v in sorted(hist.items())}
        anchored_low = sum(hist[m] for m in (1, 2) if m in hist) and anchored[key]
        extra_exact = max(0, exact[key] - 1)   # 1 tam eslesme hedef beklenir
        sev = "INFO"
        if extra_exact > 0:
            sev = "FAIL"             # fazladan tam eslesmeli hedef disi yer
        elif anchored[key] > 0:
            sev = "WARN"             # 3'-capali, uyumsuzluklu hedef disi yer
        rec["checks"].append(engine.CheckRecord(
            id="SPC-02", target=f"{pid}_{role}",
            value={"mismatch_histogram": mm_hist, "exact_sites": exact[key],
                   "anchored_mismatch_sites": anchored[key],
                   "anchor_nt": args.anchor},
            threshold={"warn": "3'-capali uyumsuzluklu hedef disi yer var",
                       "fail": "fazladan tam eslesmeli hedef disi yer"},
            severity=sev, model={"source": "genome_scan_bindings"},
            conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
            profile=eng.profile.name,
            note="genom taramasi baglanma kayitlarindan; kalibre esik yok",
        ).as_dict())

        # SPC-05: 3' terminal (son 3 nt) uyumsuzluk tipi.
        tol = term3[key].get("tolerant", 0)
        dis = term3[key].get("disruptive", 0)
        sev = "WARN" if tol > 0 else "INFO"
        rec["checks"].append(engine.CheckRecord(
            id="SPC-05", target=f"{pid}_{role}",
            value={"terminal3_tolerant": tol, "terminal3_disruptive": dis,
                   "window_nt": 3},
            threshold={"warn": "hedef disi yerde tolerant (G.T/T.G) 3' uc "
                               "uyumsuzlugu (uzamaya izin verir)"},
            severity=sev,
            model={"source": "genome_scan_bindings",
                   "tolerant_pairs": "G.T, T.G"},
            conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
            profile=eng.profile.name,
            note="tolerant 3' uc uyumsuzlugu ozgulluk riski (Kwok 1990)",
        ).as_dict())

        # SPC-06: psodogen / paralog ortusmesi.
        if pseudo_idx.total:
            names = sorted(pseudo_hits[key])
            rec["checks"].append(engine.CheckRecord(
                id="SPC-06", target=f"{pid}_{role}",
                value={"pseudogene_overlap_sites": len(names),
                       "pseudogenes": names[:10]},
                threshold={"warn": "baglanma yeri psodogenle ortusuyor"},
                severity="WARN" if names else "INFO",
                model={"source": "GENCODE gene_type *_pseudogene"},
                conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
                profile=eng.profile.name,
                note="psodogen uzerindeki baglanma hedef disi urun riski "
                     "(paralog tespiti BioMart gerektirir, kapsam disi)",
            ).as_dict())
        else:
            rec["not_evaluated"]["SPC-06"] = "psodogen anotasyonu yok"

        # SPC-07: segmental duplikasyon ortusmesi.
        if segdup_idx.total:
            rec["checks"].append(engine.CheckRecord(
                id="SPC-07", target=f"{pid}_{role}",
                value={"segdup_overlap_sites": segdup_hits[key]},
                threshold={"warn": "baglanma yeri segmental duplikasyonda"},
                severity="WARN" if segdup_hits[key] else "INFO",
                model={"source": "UCSC genomicSuperDups"},
                conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
                profile=eng.profile.name,
                note="segmental duplikasyonda baglanma kopya ayrimini zorlastirir",
            ).as_dict())
        else:
            rec["not_evaluated"]["SPC-07"] = "segdup anotasyonu yok"

    # --- hedef / hedef disi ayrimi ------------------------------------------
    # Beklenen koordinat verilmisse bir urun ancak baslangici ve uzunlugu
    # panelde beyan edilenle ayni ise HEDEF sayilir. Aksi halde, uyumsuzlugu
    # sifir olsa bile hedef disidir — tam eslesmeli hedef disi urun en kotu
    # siniftir ve "exact" etiketiyle temiz gorunmemelidir.
    expected: dict[str, tuple[int, int]] = {}
    if args.panel:
        for row in read_tsv(args.panel):
            low = {k.lower(): v for k, v in row.items()}
            pid = low.get("pair_id")
            start = low.get("left_start")
            length = low.get("expected_amplicon_length")
            if pid and start and length:
                expected[pid] = (int(start), int(length))

    def is_intended(p: dict) -> bool:
        exp = expected.get(p["pair_id"])
        if exp is None:
            return int(p["total_mismatches"]) == 0
        return (int(p["amplicon_start"]) == exp[0]
                and int(p["amplicon_length"]) == exp[1])

    prod_by_pair = collections.defaultdict(list)
    for p in products:
        prod_by_pair[p["pair_id"]].append(p)

    # --- raporlar ------------------------------------------------------------
    qc_path = os.path.join(args.out, "panel_qc.tsv")
    off_path = os.path.join(args.out, "offtarget_products.tsv")
    rec_path = os.path.join(args.out, "panel_qc_records.jsonl")

    pair_ids = sorted({b["pair_id"] for b in bindings},
                      key=lambda x: (len(x), x))
    n_pairs_with_off = 0

    with open(qc_path, "w", newline="", encoding="utf-8") as fq, \
         open(off_path, "w", newline="", encoding="utf-8") as fo, \
         open(rec_path, "w", encoding="utf-8") as fj:

        qw = csv.writer(fq, delimiter="\t", lineterminator="\n")
        qw.writerow([
            "pair_id", "tm_f_c", "tm_r_c", "theta_f_at_ta", "theta_r_at_ta",
            "sites_f", "sites_r", "exact_f", "exact_r",
            "anchored_mm_f", "anchored_mm_r",
            "products_total", "products_intended", "offtarget_total",
            "offtarget_exact", "offtarget_anchored", "min_offtarget_mismatches",
            "spc01", "spc03", "worst_primer_severity", "spc04_status",
        ])
        ow = csv.writer(fo, delimiter="\t", lineterminator="\n")
        ow.writerow(["pair_id", "chrom", "amplicon_start", "amplicon_length",
                     "total_mismatches", "class", "left_mask", "right_mask",
                     "left_anchor_intact", "right_anchor_intact"])

        for pid in pair_ids:
            f_key, r_key = (pid, "F"), (pid, "R")
            tf, tr = thermo.get(f"{pid}_F"), thermo.get(f"{pid}_R")
            prods = prod_by_pair.get(pid, [])
            intended = [p for p in prods if is_intended(p)]
            off = [p for p in prods if not is_intended(p)]
            off_exact = [p for p in off if int(p["total_mismatches"]) == 0]
            off_anchored = [
                p for p in off
                if three_prime_anchor_intact(p["left_mask"], args.anchor)
                and three_prime_anchor_intact(p["right_mask"], args.anchor)
            ]
            if off:
                n_pairs_with_off += 1

            # SPC-01: hedef tekilligi. SPC-03: pencere icinde hedef disi urun.
            spc01 = ("FAIL" if len(intended) != 1 else "INFO")
            spc03 = ("FAIL" if off_exact else
                     "WARN" if off else "INFO")

            def theta(rec):
                if not rec:
                    return ""
                c = [x for x in rec["checks"] if x["id"] == "TH-03"][0]
                return c["value"]["theta_at_ta"]

            sev = max(
                (r["worst_severity"] for r in (tf, tr) if r),
                key=lambda x: engine.SEVERITY_ORDER[x], default="INFO",
            )
            n_scored = scored[f_key] + scored[r_key]
            n_unscored = sum(unscored[f_key].values()) + sum(unscored[r_key].values())
            if n_scored:
                spc04 = f"{'partial' if n_unscored else 'evaluated'}: {n_scored} yer"
            else:
                spc04 = ("not_evaluated: desteklenen uyumsuzluk yok" if n_unscored
                         else "not_applicable: uyumsuzluklu yer yok")
            qw.writerow([
                pid,
                tf["thermo"]["tm_c"] if tf else "",
                tr["thermo"]["tm_c"] if tr else "",
                theta(tf), theta(tr),
                sum(sites[f_key].values()), sum(sites[r_key].values()),
                exact[f_key], exact[r_key],
                anchored[f_key], anchored[r_key],
                len(prods), len(intended), len(off), len(off_exact),
                len(off_anchored),
                min((int(p["total_mismatches"]) for p in off), default=""),
                spc01, spc03, sev, spc04,
            ])

            for p in sorted(off, key=lambda x: int(x["total_mismatches"]))[:50]:
                mm = int(p["total_mismatches"])
                anchored_both = (
                    three_prime_anchor_intact(p["left_mask"], args.anchor)
                    and three_prime_anchor_intact(p["right_mask"], args.anchor)
                )
                klass = ("exact_offtarget" if mm == 0
                         else "anchored_offtarget" if anchored_both
                         else "mismatched_offtarget")
                ow.writerow([
                    p["pair_id"], p["chrom"], p["amplicon_start"],
                    p["amplicon_length"], mm, klass,
                    p["left_mask"], p["right_mask"],
                    three_prime_anchor_intact(p["left_mask"], args.anchor),
                    three_prime_anchor_intact(p["right_mask"], args.anchor),
                ])

            for rec in (tf, tr):
                if rec:
                    fj.write(json.dumps(rec, ensure_ascii=False,
                                        sort_keys=True) + "\n")

            # SPC-03: cift duzeyi hedef disi urun kaydi.
            spc03_rec = engine.CheckRecord(
                id="SPC-03", target=pid,
                value={"products_total": len(prods),
                       "products_intended": len(intended),
                       "offtarget_total": len(off),
                       "offtarget_exact": len(off_exact),
                       "offtarget_anchored": len(off_anchored),
                       "min_offtarget_mismatches":
                           min((int(p["total_mismatches"]) for p in off),
                               default=None)},
                threshold={"warn": "pencere icinde hedef disi urun",
                           "fail": "tam eslesmeli hedef disi urun"},
                severity=spc03,
                model={"source": "genome_scan_products",
                       "note": "F-R urunler sozlesme araligindan; F-F/R-R ve "
                               "genisletilmis aralik ayri C++ taramasi ister"},
                conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
                profile=eng.profile.name,
                note="cift duzeyi ozgulluk; kalibre esik yok")
            fj.write(json.dumps(spc03_rec.as_dict(), ensure_ascii=False,
                                sort_keys=True) + "\n")

        # --- panel duzeyi beyan kayitlari (SPC-09, COND-01/02) -------------
        import hashlib

        def _win(path: str) -> str:
            m = re.match(r"/mnt/([a-z])/(.*)", path or "")
            return f"{m.group(1).upper()}:\\" + m.group(2).replace("/", "\\") if m else path

        # SPC-09: arama sozlesmesi beyani + indeks checksum'u.
        index_digest = None
        manifest_src = contract.get("manifest", "")
        manifest_path = _win(manifest_src)
        if manifest_path and os.path.exists(manifest_path):
            shard_sha = []
            for mrow in read_tsv(manifest_path):
                if mrow.get("sha256"):
                    shard_sha.append(f"{mrow.get('chromosome','')}:{mrow['sha256']}")
            if shard_sha:
                index_digest = hashlib.sha256(
                    "\n".join(sorted(shard_sha)).encode()).hexdigest()
        contract_ok = all(contract.get(k) for k in
                          ("anchor_length", "max_mismatches",
                           "min_amplicon_length", "max_amplicon_length"))
        spc09 = engine.CheckRecord(
            id="SPC-09", target="__panel__",
            value={"anchor_length": contract.get("anchor_length"),
                   "max_mismatches": contract.get("max_mismatches"),
                   "product_range_bp": [contract.get("min_amplicon_length"),
                                        contract.get("max_amplicon_length")],
                   "shards": contract.get("shards"),
                   "mask_convention": contract.get("mask_convention"),
                   "site_strand": contract.get("site_strand"),
                   "index_sha256": index_digest,
                   "index_manifest_found": bool(index_digest),
                   "param_set_sha256": params.registry_digest()},
            threshold={"fail": "sozlesme alanlari eksik"},
            severity="INFO" if contract_ok else "FAIL",
            model={"source": "panel.contract.tsv + manifest.sha256"},
            conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
            profile=eng.profile.name,
            note=("arama sozlesmesi beyani: capa/uyumsuzluk/urun araligi + "
                  "indeks checksum (manifest per-shard sha256 uzerinden)"))
        fj.write(json.dumps(spc09.as_dict(), ensure_ascii=False,
                            sort_keys=True) + "\n")

        # COND-01: reaksiyon kosullari beyani.
        cond01 = engine.CheckRecord(
            id="COND-01", target="__panel__", value=cond.as_dict(),
            threshold={}, severity="INFO",
            model={"source": "salt.Conditions"},
            conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
            profile=eng.profile.name,
            note="Na/K + Mg + dNTP + primer + Ta reaksiyon kosullari")
        fj.write(json.dumps(cond01.as_dict(), ensure_ascii=False,
                            sort_keys=True) + "\n")

        # COND-02: serbest Mg2+ beyani ([Mg]free = [Mg] - [dNTP]).
        cond02 = engine.CheckRecord(
            id="COND-02", target="__panel__",
            value={"magnesium_mM": round(cond.magnesium_M * 1e3, 4),
                   "dntp_mM": round(cond.dntp_M * 1e3, 4),
                   "free_magnesium_mM": round(cond.free_magnesium_M * 1e3, 4)},
            threshold={"warn": "serbest Mg2+ <= 0 (dNTP tum Mg'yi sekestre eder)"},
            severity="WARN" if cond.free_magnesium_M <= 0 else "INFO",
            model={"source": "[Mg]free = [Mg] - [dNTP]"},
            conditions=cond.as_dict(), param_set_sha256=eng.param_digest,
            profile=eng.profile.name,
            note="serbest Mg2+ Tm/salt duzeltmesinin girdisidir")
        fj.write(json.dumps(cond02.as_dict(), ensure_ascii=False,
                            sort_keys=True) + "\n")

    print(f"cift sayisi            : {len(pair_ids)}")
    print(f"baglanma kaydi         : {len(bindings)}")
    print(f"urun kaydi             : {len(products)}")
    print(f"hedef disi urunu olan  : {n_pairs_with_off} cift")
    if not expected:
        print("  NOT: beklenen koordinat verilmedi; hedef ayrimi yalniz "
              "uyumsuzluk sayisina gore yapildi. --panel ile kesinlestirin.")
    print(f"arama sozlesmesi       : capa={contract.get('anchor_length')}, "
          f"maks uyumsuzluk={contract.get('max_mismatches')}, "
          f"urun={contract.get('min_amplicon_length')}-"
          f"{contract.get('max_amplicon_length')} bp")
    print(f"parametre digest       : {params.registry_digest()[:16]}...")
    print("ciktilar:")
    for p in (qc_path, off_path, rec_path):
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
