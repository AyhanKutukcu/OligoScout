#!/usr/bin/env python3
"""Varyant katmani CLI: gnomAD v4 -> VAR-01..VAR-08.

Zincirdeki yeri:
    design_primers.py -> probe_qc.py -> oligoscout_panel -> panel_qc.py
    variant_qc.py                                              <-- burasi

Ne yapar: tasarlanan her ciftin HEDEF baglanma yerlerini (ve istenirse hedef
disi urun veren yerleri) gnomAD v4 ile tarar, varyantlari primer uzerindeki
konumlariyla birlikte degerlendirir ve kanit kaydi yazar.

Cevrimdisi calisma: cekilen her bolge `variant_cache.jsonl` dosyasina yazilir;
sonraki kosular once onbellege bakar. `--offline` verilirse ag'a hic cikilmaz.

Esikler literatur baslangic degeridir, KALIBRE DEGILDIR (REP-03).
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from oligoscout_thermo import variants as var  # noqa: E402

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
QUERY = """
query Region($chrom: String!, $start: Int!, $stop: Int!, $ds: DatasetId!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: $ds) {
      variant_id rsids pos ref alt
      genome { af ac an populations { id ac an } }
      exome  { af ac an populations { id ac an } }
    }
  }
}
"""


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# gnomAD erisimi ve onbellek
# ---------------------------------------------------------------------------
def fetch_region(chrom, start, stop, dataset, retries=3):
    """gnomAD v4 bolge sorgusu. Konumlar 1-tabanli, kapsayici."""
    payload = {"query": QUERY,
               "variables": {"chrom": var.normalize_chrom(chrom),
                             "start": int(start), "stop": int(stop),
                             "ds": dataset}}
    req = urllib.request.Request(
        GNOMAD_URL, data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"})
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                body = json.load(r)
            if body.get("errors"):
                raise RuntimeError(body["errors"][0].get("message", "gnomAD hatasi"))
            return body["data"]["region"]["variants"]
        except Exception as exc:                      # ag/limit hatalari
            last = exc
            time.sleep(2 * (attempt + 1))
    raise SystemExit(f"gnomAD sorgusu basarisiz ({chrom}:{start}-{stop}): {last}\n"
                     "Cevrimdisi calismak icin --offline verin (onbellek kullanilir).")


def parse_variants(rows, chrom):
    out = []
    for r in rows:
        src = r.get("genome") or r.get("exome")
        if not src or src.get("af") is None:
            continue
        which = "genome" if r.get("genome") else "exome"
        pops = {}
        for p in (src.get("populations") or []):
            an = p.get("an") or 0
            if an:
                pops[p["id"]] = (p.get("ac") or 0) / an
        out.append(var.Variant(
            chrom=var.normalize_chrom(chrom), pos=int(r["pos"]),
            ref=r["ref"], alt=r["alt"], af=float(src["af"]),
            ac=int(src.get("ac") or 0), an=int(src.get("an") or 0),
            rsid=(r.get("rsids") or [""])[0] or "",
            source=f"gnomad_r4_{which}", pop_af=pops))
    return out


def load_cache(path):
    regions, variants = {}, []
    if not os.path.exists(path):
        return regions, variants
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("_region"):
                regions[rec["_region"]] = True
                continue
            variants.append(var.Variant(
                chrom=rec["chrom"], pos=rec["pos"], ref=rec["ref"],
                alt=rec["alt"], af=rec["af"], ac=rec.get("ac", 0),
                an=rec.get("an", 0), rsid=rec.get("rsid", ""),
                source=rec.get("source", "gnomad_r4"),
                pop_af=rec.get("pop_af", {})))
    return regions, variants


def append_cache(path, key, variants):
    with open(path, "a", encoding="utf-8", newline="\n") as fh:
        for v in variants:
            fh.write(json.dumps({
                "chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": v.alt,
                "af": v.af, "ac": v.ac, "an": v.an, "rsid": v.rsid,
                "source": v.source, "pop_af": v.pop_af}, ensure_ascii=False) + "\n")
        fh.write(json.dumps({"_region": key}) + "\n")


def merge_intervals(intervals, pad=0):
    """[(chrom, start, end)] -> ortusenleri birlestirilmis liste."""
    by_chrom = collections.defaultdict(list)
    for c, s, e in intervals:
        by_chrom[var.normalize_chrom(c)].append((max(0, s - pad), e + pad))
    merged = []
    for c, spans in by_chrom.items():
        spans.sort()
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((c, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((c, cur_s, cur_e))
    return merged


# ---------------------------------------------------------------------------
def build_sites(design, products, probes, scope):
    """Tasarlanan ciftlerin hedef (ve istenirse hedef disi) baglanma yerleri."""
    by_coord = collections.defaultdict(list)
    for p in products:
        by_coord[(p["pair_id"], int(p["amplicon_start"]),
                  int(p["amplicon_end_exclusive"]))].append(p)

    sites, amplicons, missing = [], [], []
    for row in design:
        pid = row["pair_id"]
        f, r = row["primer1"], row["primer2"]
        if "left_start" not in row:
            missing.append(pid)
            continue
        key = (pid, int(row["left_start"]), int(row["right_start"]))
        hits = by_coord.get(key)
        if not hits:
            missing.append(pid)
            continue
        p = hits[0]
        chrom = p["chrom"]
        lf = int(p["left_position"]); rt = int(p["right_position"])
        f_ref = f if p.get("left_mismatches") == "0" else ""
        r_ref = r if p.get("right_mismatches") == "0" else ""
        sites.append((pid, "F", var.Site(f"{pid}_F", chrom, lf, lf + len(f), "+", f, f_ref)))
        sites.append((pid, "R", var.Site(f"{pid}_R", chrom, rt, rt + len(r), "-", r, r_ref)))
        amplicons.append((pid, chrom, int(p["amplicon_start"]),
                          int(p["amplicon_end_exclusive"])))
        if pid in probes:
            s, e = probes[pid]
            gs = int(row["left_start"]) + s
            sites.append((pid, "P", var.Site(f"{pid}_P", chrom, gs,
                                             gs + (e - s), "+", "")))

    if scope == "products":
        designed = {(pid, s, e) for pid, s, e in
                    ((p[0], p[2], p[3]) for p in amplicons)}
        for p in products:
            key = (p["pair_id"], int(p["amplicon_start"]),
                   int(p["amplicon_end_exclusive"]))
            if key in designed:
                continue
            pid = p["pair_id"]
            row = next((d for d in design if d["pair_id"] == pid), None)
            if row is None:
                continue
            f, r = row["primer1"], row["primer2"]
            lf = int(p["left_position"]); rt = int(p["right_position"])
            tag = f'{p["chrom"]}:{lf}'
            sites.append((pid, "F_off", var.Site(f"{pid}_F@{tag}", p["chrom"],
                                                 lf, lf + len(f), "+", f, "")))
            sites.append((pid, "R_off", var.Site(f"{pid}_R@{tag}", p["chrom"],
                                                 rt, rt + len(r), "-", r, "")))
    return sites, amplicons, missing


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True, help="RUN_PANEL cikti klasoru")
    ap.add_argument("--design", required=True, help="design_panel.tsv klasoru")
    ap.add_argument("--out", help="varsayilan: --panel")
    ap.add_argument("--scope", choices=("intended", "products"),
                    default="intended",
                    help="intended: yalniz hedef yerler; products: hedef disi "
                         "urun veren yerler de taranir")
    ap.add_argument("--dataset", default="gnomad_r4")
    ap.add_argument("--cache", help="varsayilan: <out>/variant_cache.jsonl")
    ap.add_argument("--offline", action="store_true",
                    help="ag'a cikma, yalniz onbellegi kullan")
    ap.add_argument("--af-common", type=float, default=0.01)
    ap.add_argument("--af-warn", type=float, default=1e-3)
    ap.add_argument("--three-prime-nt", type=int, default=5)
    args = ap.parse_args()

    out_dir = args.out or args.panel
    os.makedirs(out_dir, exist_ok=True)
    cache_path = args.cache or os.path.join(out_dir, "variant_cache.jsonl")

    design = read_tsv(os.path.join(args.design, "design_panel.tsv"))
    products = read_tsv(os.path.join(args.panel, "panel.products.tsv"))
    probes = {}
    pp = os.path.join(args.design, "probe_panel.tsv")
    if os.path.exists(pp):
        for r in read_tsv(pp):
            if r.get("probe") and r.get("start_in_amplicon"):
                probes[r["pair_id"]] = (int(r["start_in_amplicon"]),
                                        int(r["end_in_amplicon"]))

    sites, amplicons, missing = build_sites(design, products, probes, args.scope)
    if not sites:
        raise SystemExit("hedef baglanma yeri bulunamadi; design_panel.tsv "
                         "genom koordinati tasiyor mu? (FASTA basligi chrN:a-b)")

    cons = var.VariantConstraints(af_common=args.af_common,
                                  af_warn=args.af_warn,
                                  three_prime_nt=args.three_prime_nt)

    cached_regions, pool = load_cache(cache_path)
    regions = merge_intervals([(s.chrom, s.start, s.end) for _, _, s in sites]
                              + [(c, a, b) for _, c, a, b in amplicons], pad=5)
    fetched = 0
    for chrom, start, end in regions:
        key = f"{chrom}:{start}-{end}"
        if key in cached_regions:
            continue
        if args.offline:
            print(f"  onbellekte yok, atlandi (offline): {key}")
            continue
        rows = fetch_region(chrom, start + 1, end, args.dataset)
        new = parse_variants(rows, chrom)
        append_cache(cache_path, key, new)
        pool.extend(new)
        fetched += len(new)
        time.sleep(0.5)

    conditions = {"dataset": args.dataset, "scope": args.scope,
                  "reference_genome": "GRCh38"}
    rows_qc, records, site_rows = [], [], []
    for pid, role, site in sites:
        checks = var.evaluate_site(site, pool, cons, conditions=conditions,
                                   probe=(role == "P"))
        hits = var.overlapping(pool, site.chrom, site.start, site.end)
        site_rows.append([pid, role, site.site_id, site.chrom, site.start,
                          site.end, site.strand, len(hits),
                          round(var.max_af(hits), 6),
                          var.worst_severity(checks)])
        for c in checks:
            rows_qc.append([pid, role, c.id, c.target, c.severity,
                            json.dumps(c.value, ensure_ascii=False,
                                       separators=(",", ":")),
                            json.dumps(c.threshold, ensure_ascii=False,
                                       separators=(",", ":")), c.note])
            records.append(c.as_dict())

    for pid, chrom, a, b in amplicons:
        own = [s for p, _, s in sites if p == pid]
        checks = var.evaluate_amplicon(f"{pid}_amplicon", chrom, a, b, own, pool,
                                       cons, conditions=conditions)
        for c in checks:
            rows_qc.append([pid, "amplicon", c.id, c.target, c.severity,
                            json.dumps(c.value, ensure_ascii=False,
                                       separators=(",", ":")),
                            json.dumps(c.threshold, ensure_ascii=False,
                                       separators=(",", ":")), c.note])
            records.append(c.as_dict())

    qc_path = os.path.join(out_dir, "variant_qc.tsv")
    with open(qc_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pair_id", "role", "check", "target", "severity", "value",
                    "threshold", "note"])
        w.writerows(rows_qc)

    sites_path = os.path.join(out_dir, "variant_sites.tsv")
    with open(sites_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pair_id", "role", "site_id", "chrom", "start", "end",
                    "strand", "n_variants", "max_af", "worst_severity"])
        w.writerows(site_rows)

    rec_path = os.path.join(out_dir, "variant_qc_records.jsonl")
    with open(rec_path, "w", encoding="utf-8", newline="\n") as fh:
        for r in records:
            fh.write(json.dumps(r, ensure_ascii=False) + "\n")

    by_check = collections.Counter((r[2], r[4]) for r in rows_qc)
    worst_pair = collections.Counter()
    for r in rows_qc:
        if r[4] != "INFO":
            worst_pair[r[0]] += 1
    print(f"taranan bolge        : {len(regions)} "
          f"({'onbellek' if not fetched else str(fetched) + ' yeni varyant'})")
    print(f"baglanma yeri        : {len(sites)} ({args.scope})")
    print(f"varyant havuzu       : {len(pool)}")
    if missing:
        print(f"koordinatsiz cift    : {len(missing)} ({', '.join(missing[:5])})")
    sev_rows = sorted((k for k in by_check if k[1] != "INFO"))
    print("WARN/FAIL dagilimi   : " +
          (", ".join(f"{c}:{s}={by_check[(c, s)]}" for c, s in sev_rows)
           if sev_rows else "yok"))
    print(f"etkilenen cift       : {len(worst_pair)} / {len(design)}")
    var04_scored = sum(r["value"].get("n_scored", 0) for r in records if r["id"] == "VAR-04")
    print(f"VAR-04 hesaplanan    : {var04_scored} varyant-baglanma yeri")
    print("esikler kalibre degildir (REP-03)")
    print("ciktilar:")
    for p in (sites_path, qc_path, rec_path, cache_path):
        if os.path.exists(p):
            print(f"  {p}  ({os.path.getsize(p)} bayt)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
