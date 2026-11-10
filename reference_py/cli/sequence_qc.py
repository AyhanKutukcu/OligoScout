#!/usr/bin/env python3
"""Dizi/kompozisyon/cift/amplikon kontrolleri (SEQ, CMP, SS-02, PAIR, AMP).

Zincirdeki yeri:
    design_primers.py  -> cift adaylari (bu betigin girdisi)
    sequence_qc.py     -> kendine yeten dizi ve cift kontrolleri   <-- burasi
    panel_qc.py        -> ozgulluk taramasi (ayri adim)

Girdi : --design klasoru (design_pairs.tsv) ve --target FASTA.
Cikti : seqpair_qc.tsv (satir basina bir kontrol),
        seqpair_panel.tsv (cift basina ozet karar),
        seqpair_qc_records.jsonl (REP-02 kanit kaydi).

Esikler literatur kilavuz degeridir, KALIBRE DEGILDIR (REP-03).
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import (bisulfite, nn, params, profiles,  # noqa: E402
                               salt, seqpair)


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path: str) -> str:
    parts = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if parts:
                    break
                continue
            parts.append(line.strip())
    return nn.validate("".join(parts))


def read_fasta_raw(path: str) -> tuple[str, str]:
    """(ham dizi, baslik). SEQ-04 icin N ve soft-mask KORUNUR."""
    parts, header = [], ""
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if parts:
                    break
                header = line[1:].strip()
                continue
            parts.append(line.strip())
    return "".join(parts), header


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True,
                    help="design_pairs.tsv iceren klasor")
    ap.add_argument("--target", required=True, help="tasarimda kullanilan FASTA")
    ap.add_argument("--out", help="cikti klasoru (varsayilan: --design)")
    ap.add_argument("--profile", default="ctdna",
                    choices=sorted(profiles.PROFILES))
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=1.5)
    ap.add_argument("--dntp", type=float, default=0.6)
    ap.add_argument("--primer-nm", type=float, default=200.0)
    ap.add_argument("--ta", type=float, default=60.0)
    ap.add_argument("--salt-model", default="owczarzy_2008",
                    choices=list(salt.MODELS))
    ap.add_argument("--product-bp", default="",
                    help="urun boyu penceresi min:max (varsayilan: profil)")
    ap.add_argument("--gene", default="",
                    help="hedef gen adi (CT-06 CHIP uyarisi icin)")
    ap.add_argument("--msp", action="store_true",
                    help="bisulfit: MSP (metilasyona-ozgu) primer; CpG serbest "
                         "(BIS-01 CpG'yi hata saymaz)")
    ap.add_argument("--probes", default="",
                    help="DIM-05 icin prob paneli: pair_id<TAB>probe TSV'si "
                         "(probe_qc.py ciktisi probe_panel.tsv da olur)")
    ap.add_argument("--contigs", default="",
                    help="SPC-08 icin kucuk contig FASTA'lari (virgulle; "
                         "ornegin chrM.fa). Buyuk ALT icin C++ gerekir.")
    ap.add_argument("--tails", default="",
                    help="DIM-07 kuyruk-dahil analiz: fwd_tail:rev_tail "
                         "(M13/NGS adaptor kuyruklari, 5' ucuna eklenir)")
    args = ap.parse_args()

    out_dir = args.out or args.design
    os.makedirs(out_dir, exist_ok=True)
    target = read_fasta(args.target)
    raw_target, header = read_fasta_raw(args.target)
    pairs = read_tsv(os.path.join(args.design, "design_pairs.tsv"))

    prof = profiles.get(args.profile)
    pc = getattr(prof, "constraints", None)
    pbp = ((pc.product_min, pc.product_max)
           if pc and hasattr(pc, "product_min") else (60, 200))
    if args.product_bp:
        lo, hi = (int(x) for x in args.product_bp.split(":"))
        pbp = (lo, hi)
    cons = seqpair.SeqPairConstraints(product_bp=tuple(pbp))
    cond = salt.Conditions(
        monovalent_M=args.na / 1e3, magnesium_M=args.mg / 1e3,
        dntp_M=args.dntp / 1e3, primer_M=args.primer_nm / 1e9,
        anneal_c=args.ta,
    )
    digest = params.registry_digest()
    cond_dict = cond.as_dict()

    # SEQ-06 icin referans FASTA tam dosya SHA-256'si (izlenebilirlik).
    import hashlib
    with open(args.target, "rb") as _fh:
        fasta_sha256 = hashlib.sha256(_fh.read()).hexdigest()

    # DIM-07 icin opsiyonel kuyruklar (M13/NGS adaptorleri).
    fwd_tail, rev_tail = "", ""
    if args.tails:
        parts = args.tails.split(":")
        fwd_tail = parts[0].strip().upper()
        rev_tail = (parts[1].strip().upper() if len(parts) > 1 else "")

    # SPC-08 icin opsiyonel kucuk contig'ler (chrM vb.).
    contigs: dict[str, str] = {}
    for path in [c for c in args.contigs.split(",") if c.strip()]:
        name = os.path.splitext(os.path.basename(path))[0]
        contigs[name] = read_fasta_raw(path)[0].upper()

    # DIM-05 icin opsiyonel prob haritasi (pair_id -> prob dizisi).
    probe_map: dict[str, str] = {}
    if args.probes:
        for row in read_tsv(args.probes):
            seq = row.get("probe") or ""
            if seq:
                probe_map[row.get("pair_id") or row.get("id") or ""] = seq

    qc_rows: list[list] = []
    panel_rows: list[list] = []
    records: list[dict] = []
    verdicts: collections.Counter = collections.Counter()
    by_check: collections.Counter = collections.Counter()

    for idx, row in enumerate(pairs):
        pid = row.get("pair_id") or row.get("id") or f"PAIR_{idx:03d}"
        fwd = nn.validate(row["left_primer"])
        rev = nn.validate(row["right_primer"])
        p_start = int(row["product_start"])
        p_end = int(row["product_end"])
        amplicon = target[p_start:p_end]

        checks: list[dict] = []
        for role, raw_seq in (("forward", row["left_primer"]),
                              ("reverse", row["right_primer"])):
            # SEQ-01 ham dizi uzerinden (alfabe + modifikasyon notasyonu).
            checks.append({**seqpair.check_seq01(raw_seq), "role": role})
        for role, seq in (("forward", fwd), ("reverse", rev)):
            for ch in seqpair.check_primer(seq, cons):
                ch = {**ch, "role": role}
                checks.append(ch)
        for ch in seqpair.check_pair(fwd, rev, p_start, p_end, amplicon,
                                     cond=cond, salt_model=args.salt_model,
                                     cons=cons):
            checks.append({**ch, "role": "pair"})
        for ch in seqpair.check_template(fwd, rev, p_start, p_end, raw_target):
            checks.append({**ch, "role": "template"})
        checks.append({**seqpair.check_amplicon_structure(amplicon, args.ta, cons),
                       "role": "amplicon"})
        if args.profile == "rt_qpcr":
            for role, seq in (("forward", fwd), ("reverse", rev)):
                checks.append({**seqpair.check_rt05(seq), "role": role})
        if args.profile == "bisulfite":
            # BIS-01/03 donusmus primer uzerinden; BIS-02 ORIJINAL hedef
            # footprint'inden (--target donusmemis referans olmali).
            msp = bool(args.msp)
            fwd_fp = raw_target[p_start:p_start + len(fwd)].upper()
            rev_fp = seqpair.reverse_complement_raw(
                raw_target[p_end - len(rev):p_end]).upper()
            for role, seq, fp in (("forward", fwd, fwd_fp),
                                  ("reverse", rev, rev_fp)):
                checks.append({**bisulfite.check_bis01(seq, msp=msp), "role": role})
                checks.append({**bisulfite.check_bis02(fp), "role": role})
                checks.append({**bisulfite.check_bis03(seq), "role": role})
        if contigs:
            for role, seq in (("forward", fwd), ("reverse", rev)):
                checks.append({**seqpair.check_spc08(seq, contigs), "role": role})
        if fwd_tail or rev_tail:
            checks.append({**seqpair.check_dim07(
                fwd, rev, fwd_tail, rev_tail), "role": "pair"})
        if pid in probe_map:
            checks.append({**seqpair.check_primer_probe_dimers(
                fwd, rev, probe_map[pid]), "role": "probe"})
            # SS-05: prob hairpin Tm (prob verildiginde).
            checks.append({**seqpair.check_ss05(probe_map[pid], args.ta, cons),
                           "role": "probe"})

        verdict = seqpair.worst(checks)
        verdicts[verdict] += 1
        panel_rows.append([pid, fwd, rev, len(amplicon), verdict])
        for ch in checks:
            if ch["severity"] != "INFO":
                by_check[(ch["id"], ch["severity"])] += 1
            qc_rows.append([
                pid, ch.get("role", ""), ch["id"], ch["severity"],
                json.dumps(ch["value"], ensure_ascii=False, separators=(",", ":")),
                json.dumps(ch["threshold"], ensure_ascii=False,
                           separators=(",", ":")),
                ch.get("note", ""),
            ])
        records.append({
            "pair_id": pid, "forward": fwd, "reverse": rev,
            "amplicon_bp": len(amplicon), "worst_severity": verdict,
            "checks": checks, "conditions": cond_dict,
            "param_set_sha256": digest, "profile": args.profile,
            "engine_version": "oligoscout-thermo-ref/seqpair",
            "calibration_status": "REP-03 / KALIBRE DEGIL",
        })

    # CT-06 panel/hedef duzeyinde: her cifte degil, bir kez uygulanir.
    chip = seqpair.check_chip_gene(args.gene or None, args.profile)
    if chip["severity"] != "INFO":
        by_check[(chip["id"], chip["severity"])] += 1
    qc_rows.append(["__panel__", "target", chip["id"], chip["severity"],
                    json.dumps(chip["value"], ensure_ascii=False,
                               separators=(",", ":")),
                    json.dumps(chip["threshold"], ensure_ascii=False,
                               separators=(",", ":")),
                    chip.get("note", "")])
    records.append({"scope": "panel", "check": chip, "target_header": header,
                    "param_set_sha256": digest, "profile": args.profile})

    # SEQ-06 panel/hedef duzeyinde: referans kimligi (build + yama + sha256).
    seq06 = seqpair.check_seq06(header, fasta_sha256, len(raw_target))
    if seq06["severity"] != "INFO":
        by_check[(seq06["id"], seq06["severity"])] += 1
    qc_rows.append(["__panel__", "reference", seq06["id"], seq06["severity"],
                    json.dumps(seq06["value"], ensure_ascii=False,
                               separators=(",", ":")),
                    json.dumps(seq06["threshold"], ensure_ascii=False,
                               separators=(",", ":")),
                    seq06.get("note", "")])
    records.append({"scope": "panel", "check": seq06,
                    "param_set_sha256": digest, "profile": args.profile})

    panel_path = os.path.join(out_dir, "seqpair_panel.tsv")
    with open(panel_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pair_id", "forward_primer", "reverse_primer",
                    "amplicon_bp", "verdict"])
        w.writerows(panel_rows)

    qc_path = os.path.join(out_dir, "seqpair_qc.tsv")
    with open(qc_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pair_id", "role", "check", "severity", "value",
                    "threshold", "note"])
        w.writerows(qc_rows)

    rec_path = os.path.join(out_dir, "seqpair_qc_records.jsonl")
    with open(rec_path, "w", encoding="utf-8", newline="\n") as fh:
        for r in records:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")

    print(f"cift sayisi    : {len(pairs)}")
    print("karar          : " +
          ", ".join(f"{k}={v}" for k, v in sorted(verdicts.items())))
    if by_check:
        print("WARN/FAIL      : " +
              ", ".join(f"{c}:{s}={n}" for (c, s), n in sorted(by_check.items())))
    print("esikler kalibre degildir (REP-03)")
    print("ciktilar:")
    for p in (panel_path, qc_path, rec_path):
        print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
