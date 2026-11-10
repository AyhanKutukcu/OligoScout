#!/usr/bin/env python3
"""Tuz duzeltmesi katmanini primer3'un `oligotm` ikilisine karsi denetler.

NEDEN AYRI BIR BETIK. `param_crosscheck.py` tablo tablo (sayi sayi)
karsilastirir; `salt.py` ise tablo degil FORMUL tasir (Owczarzy 2004 /
Owczarzy 2008 / von Ahsen 2001 katsayilari). primer3 bu katsayilari ayri
bir dosyada yayinlamaz, `oligotm.c` icinde kodlar. Bu yuzden denetim
SAYISAL yapilir: ayni dizi ve ayni kosullar iki uygulamadan gecirilir.

IKI ASAMA. Karsilastirma tek basina yeterli degildir, cunku Tm iki
katmandan gelir: (1) NN + konsantrasyon sozlesmesi, (2) tuz duzeltmesi.

  Asama A  1 M Na+, Mg = 0, dNTP = 0 -> tuz duzeltmesi etkisizdir.
           Burada cikan fark KONSANTRASYON SOZLESMESI farkidir, tuz
           katsayisi degil.
  Asama B  gercekci tuz izgarasi. Asama A tutuyorsa buradaki fark
           dogrudan tuz katmanina aittir.

Asama A tutmuyorsa betik bunu acikca yazar ve B'yi "kesin degil" olarak
isaretler; sessizce gecmez.

Model eslesmesi: `oligotm -tp 1 -sc 2` = SantaLucia 1998 birlesik tablosu +
Owczarzy tuz duzeltmesi. Mg yokken Python karsiligi `owczarzy_2004`,
Mg varken `owczarzy_2008`'dir; betik ikisini de raporlar.

Kullanim:
    python3 cli/salt_crosscheck.py [--oligotm YOL] [--out out/salt_crosscheck]

primer3 ikilisi yoksa cikis kodu 0 ve "ATLANDI" raporu (cevrimdisi
calisabilirlik korunur).
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from oligoscout_thermo import concentration, nn, salt  # noqa: E402
from thermo_report import DEMO_PANEL  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(HERE)
DEFAULT_OLIGOTM = os.path.join(REPO_ROOT, "OligoScout", "third_party",
                               "primer3", "src", "oligotm")

# (etiket, Na mM, Mg mM, dNTP mM)
STAGE_A = (("1M_Na", 1000.0, 0.0, 0.0),)
STAGE_B = (
    ("Na50", 50.0, 0.0, 0.0),
    ("Na50_Mg1.5_dNTP0.6", 50.0, 1.5, 0.6),
    ("Na50_Mg3_dNTP0.8", 50.0, 3.0, 0.8),
    ("Na100_Mg2_dNTP0.2", 100.0, 2.0, 0.2),
    ("Na20_Mg1_dNTP0.2", 20.0, 1.0, 0.2),
)
PRIMER_NM = 200.0
TOL_C = 0.05          # 0.05 C: oligotm iki ondalik basar


def oligotm_tm(binary: str, seq: str, na_mM: float, mg_mM: float,
               dntp_mM: float, primer_nM: float) -> float | None:
    cmd = [binary, "-mv", f"{na_mM}", "-dv", f"{mg_mM}", "-n", f"{dntp_mM}",
           "-d", f"{primer_nM}", "-tp", "1", "-sc", "2", seq]
    try:
        done = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        print(f"  oligotm calistirilamadi: {exc}", file=sys.stderr)
        return None
    if done.returncode != 0:
        return None
    text = done.stdout.strip().split()
    for token in text:
        try:
            return float(token)
        except ValueError:
            continue
    return None


def python_tm(seq: str, na_mM: float, mg_mM: float, dntp_mM: float,
              primer_nM: float, model: str) -> float:
    cond = salt.Conditions(monovalent_M=na_mM / 1e3, magnesium_M=mg_mM / 1e3,
                           dntp_M=dntp_mM / 1e3, primer_M=primer_nM / 1e9,
                           anneal_c=60.0)
    return concentration.build(seq, cond=cond, salt_model=model).tm_c


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--oligotm", default=DEFAULT_OLIGOTM)
    ap.add_argument("--out", default=os.path.join(HERE, "out", "salt_crosscheck"))
    ap.add_argument("--primer-nm", type=float, default=PRIMER_NM)
    ap.add_argument("--tol", type=float, default=TOL_C)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    md_path = os.path.join(args.out, "salt_crosscheck_summary.md")
    tsv_path = os.path.join(args.out, "salt_crosscheck.tsv")

    if not (os.path.isfile(args.oligotm) and os.access(args.oligotm, os.X_OK)):
        message = (f"primer3 `oligotm` ikilisi yok ya da calistirilabilir "
                   f"degil: {args.oligotm}")
        with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# Tuz katmani capraz denetimi\n\nATLANDI: " + message + "\n")
            fh.write("\nIkiliyi uretmek icin: `RUN_SETUP.cmd --only 4`\n")
        print("ATLANDI: " + message)
        print(f"  {md_path}")
        return 0

    sequences = [(pid, nn.validate(seq)) for pid, seq in DEMO_PANEL
                 if len(seq) >= 8]
    rows: list[dict] = []
    for stage, grid in (("A", STAGE_A), ("B", STAGE_B)):
        for label, na, mg, dntp in grid:
            model = "owczarzy_2004" if mg <= 0 else "owczarzy_2008"
            for pid, seq in sequences:
                theirs = oligotm_tm(args.oligotm, seq, na, mg, dntp,
                                    args.primer_nm)
                if theirs is None:
                    rows.append({"stage": stage, "kosul": label, "id": pid,
                                 "model": model, "python_tm_c": "",
                                 "oligotm_tm_c": "", "delta_c": "",
                                 "status": "ORACLE_HATASI"})
                    continue
                ours = python_tm(seq, na, mg, dntp, args.primer_nm, model)
                delta = ours - theirs
                rows.append({
                    "stage": stage, "kosul": label, "id": pid, "model": model,
                    "python_tm_c": f"{ours:.4f}",
                    "oligotm_tm_c": f"{theirs:.4f}",
                    "delta_c": f"{delta:+.4f}",
                    "status": "MATCH" if abs(delta) <= args.tol else "DIFF",
                })

    with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    def stats(stage: str) -> tuple[int, int, float, float]:
        group = [r for r in rows if r["stage"] == stage and r["delta_c"]]
        deltas = [float(r["delta_c"]) for r in group]
        diff = sum(1 for r in group if r["status"] == "DIFF")
        if not deltas:
            return 0, diff, 0.0, 0.0
        mean = sum(deltas) / len(deltas)
        spread = max(abs(d - mean) for d in deltas)
        return len(group), diff, mean, spread

    a_n, a_diff, a_mean, a_spread = stats("A")
    b_n, b_diff, b_mean, b_spread = stats("B")
    a_clean = a_diff == 0
    failed = sum(1 for r in rows if r["status"] == "ORACLE_HATASI")
    if a_n + b_n == 0:
        # Ikili bulundu ama hic calismadi (ornegin WSL disinda bir ELF).
        # Bunu "sapan yok" diye raporlamak yanlis olurdu.
        with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# Tuz katmani capraz denetimi\n\n")
            fh.write(f"**ORACLE CALISMADI.** `{args.oligotm}` bulundu ama "
                     f"{failed} cagrinin hicbiri sonuc vermedi. Betigi WSL "
                     "icinden calistirin (ikili Linux ELF'tir).\n")
        print(f"ORACLE CALISMADI: {failed} cagri basarisiz. "
              "Betigi WSL icinden calistirin.")
        print(f"  {md_path}")
        return 2

    with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# Tuz katmani capraz denetimi (`salt.py` <-> primer3 oligotm)\n\n")
        fh.write(f"- oracle: `{args.oligotm}` (`-tp 1 -sc 2`)\n")
        fh.write(f"- tolerans: {args.tol} C (oligotm iki ondalik basar)\n")
        fh.write(f"- dizi: {len(sequences)}, primer {args.primer_nm} nM\n\n")
        fh.write("| Asama | Ne olcer | Karsilastirma | Sapan | Ortalama fark | "
                 "Ortalamadan sapma |\n|---|---|---:|---:|---:|---:|\n")
        fh.write(f"| A | NN + konsantrasyon sozlesmesi (1 M Na) | {a_n} | "
                 f"{a_diff} | {a_mean:+.4f} C | {a_spread:.4f} C |\n")
        fh.write(f"| B | tuz duzeltmesi | {b_n} | {b_diff} | {b_mean:+.4f} C | "
                 f"{b_spread:.4f} C |\n\n")
        if a_clean:
            fh.write("Asama A temiz: iki uygulama 1 M Na'da ayni Tm'i veriyor, "
                     "yani konsantrasyon sozlesmesi ayni. Asama B'deki fark "
                     "dogrudan TUZ KATMANINA aittir.\n")
        else:
            fh.write("**Asama A tutmuyor.** Iki uygulama 1 M Na'da bile "
                     "ayrisiyor; bu bir tuz katsayisi farki DEGIL, "
                     "konsantrasyon sozlesmesi farkidir. Asama B sonucu bu "
                     "yuzden tuz katmani icin KESIN DEGILDIR; once A'daki "
                     "sabit fark aciklanmalidir.\n")
        worst = sorted((r for r in rows if r["status"] == "DIFF"),
                       key=lambda r: -abs(float(r["delta_c"])))[:15]
        fh.write("\n## En buyuk farklar\n\n")
        if worst:
            fh.write("| Asama | Kosul | id | Model | Python | oligotm | Fark |\n")
            fh.write("|---|---|---|---|---:|---:|---:|\n")
            for r in worst:
                fh.write(f"| {r['stage']} | {r['kosul']} | {r['id']} | "
                         f"{r['model']} | {r['python_tm_c']} | "
                         f"{r['oligotm_tm_c']} | {r['delta_c']} |\n")
        else:
            fh.write("Yok.\n")
        fh.write("\nBu denetim `salt.py` kunyesini 'verified' YAPMAZ; "
                 "katsayilarin birincil makaleden ikinci aktarimi ayri bir "
                 "istir.\n")

    print(f"asama A (1 M Na) : {a_n} karsilastirma, sapan {a_diff}, "
          f"ortalama fark {a_mean:+.4f} C")
    print(f"asama B (tuz)    : {b_n} karsilastirma, sapan {b_diff}, "
          f"ortalama fark {b_mean:+.4f} C")
    if not a_clean:
        print("UYARI: asama A tutmuyor; B tuz katmani icin kesin degil.")
    print("ciktilar:")
    for path in (tsv_path, md_path):
        print(f"  {path}  ({os.path.getsize(path)} bayt)")
    return 1 if (a_clean and b_diff) else 0


if __name__ == "__main__":
    raise SystemExit(main())
