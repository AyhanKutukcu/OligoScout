"""Cross-check catalog IDs and record evidence, without treating coverage as PASS."""
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "reference_py"
OUT = REF / "out"
AUDIT = OUT / "qc_catalog_audit"
RUN = AUDIT / "regresyon"

# Some reports use rule-specific columns instead of an ID in each record.
MEASURES = {
    "AMP-02": ("melt/amplicon_melt.tsv", "amp02"),
    "QPCR-02": ("melt/amplicon_melt.tsv", "qpcr02"),
    "CT-01": ("design_ctdna/design_pairs.tsv", "product_length"),
    "CT-02": ("design_ctdna/ctdna_detectability.tsv", "p_intact_fragment"),
    "CT-03": ("design_ctdna/ctdna_detectability.tsv", "p_detect"),
    "PAIR-01": ("design_ctdna/design_pairs.tsv", "delta_tm_c"),
}
KEYS = ("id", "check", "rule", "kural", "check_id")


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def nested_records(value):
    if isinstance(value, dict):
        yield value
        for item in value.values():
            yield from nested_records(item)
    elif isinstance(value, list):
        for item in value:
            yield from nested_records(item)


def records(path):
    if path.suffix == ".tsv":
        yield from read_tsv(path)
    elif path.suffix == ".jsonl":
        with path.open(encoding="utf-8-sig") as stream:
            for line in stream:
                if line.strip():
                    yield from nested_records(json.loads(line))


def index(paths, ids, base):
    by_rule = defaultdict(set)
    for path in sorted(paths):
        if path.suffix not in {".tsv", ".jsonl"}:
            continue
        for row in records(path):
            for key in KEYS:
                rid = row.get(key)
                if isinstance(rid, str) and rid in ids:
                    by_rule[rid].add(path.relative_to(base).as_posix())
    return by_rule


def references(files, ids):
    by_rule = defaultdict(set)
    lookup = {rid.replace("-", "").lower(): rid for rid in ids}
    pattern = re.compile(r"(?<![a-z0-9])([a-z]+)-?(\d{2})(?![a-z0-9])", re.I)
    for file in sorted(files):
        for no, line in enumerate(file.read_text(encoding="utf-8-sig").splitlines(), 1):
            for fam, num in pattern.findall(line):
                rid = lookup.get(fam.lower() + num)
                if rid and not any(ref.startswith(file.relative_to(ROOT).as_posix() + ":") for ref in by_rule[rid]):
                    by_rule[rid].add(f"{file.relative_to(ROOT).as_posix()}:{no}")
    return by_rule


def main():
    inventory = list(read_tsv(OUT / "durum" / "qc_katalog_kapsami.tsv"))
    guide = {r["id"]: r for r in read_tsv(ROOT / "web/src/qc_rules.tsv")}
    bundle = json.loads((ROOT / "web/oligoscout_data.json").read_text(encoding="utf-8"))
    ids = {r["kural"] for r in inventory}
    assert len(ids) == len(inventory) == len(guide) == len(bundle["catalog"]) == 91
    assert ids == set(guide) == {r["id"] for r in bundle["catalog"]}
    assert len({r["aile"] for r in inventory}) == 19
    bundled_status = {r["id"]: r["status"] for r in bundle["catalog"]}
    assert all(bundled_status[r["kural"]] == r["durum"] for r in inventory)

    fresh = index(RUN.rglob("*"), ids, RUN)
    for rid, (file, column) in MEASURES.items():
        if any(row.get(column) not in (None, "") for row in read_tsv(RUN / file)):
            fresh[rid].add(f"{file} [{column}]")

    threshold_note = RUN / "kalibrasyon" / "esik_karar_notu.md"
    if threshold_note.is_file() and "REP-03" in threshold_note.read_text(encoding="utf-8"):
        fresh["REP-03"].add(threshold_note.relative_to(RUN).as_posix())

    # Exclude the inventory and this audit's own files as sources of evidence.
    archive = index((p for p in OUT.glob("*/*") if p.is_file() and
                     p.parent.name not in {"durum", AUDIT.name}), ids, OUT)
    sources = references([*(REF / "oligoscout_thermo").glob("*.py"),
                          *(p for p in (REF / "cli").glob("*.py") if p != Path(__file__).resolve()),
                          ROOT / "web/local_engine.py", ROOT / "web/src/core.js", ROOT / "web/src/views2.js"], ids)
    tests = references((REF / "tests").glob("test_*.py"), ids)
    assert all(sources[rid] for rid in ids), sorted(rid for rid in ids if not sources[rid])

    with (AUDIT / "qc_kural_dogrulama.tsv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("kural", "aile", "sade_ad", "katalog_durumu", "yeni_regresyon_kaydi",
                         "arsiv_kaydi", "testte_kimlik_veya_sembol", "kodda_kimlik_veya_sembol"))
        for r in inventory:
            rid = r["kural"]
            writer.writerow((rid, r["aile"], guide[rid]["name_tr"], r["durum"],
                             "; ".join(sorted(fresh[rid])), "; ".join(sorted(archive[rid])),
                             "; ".join(sorted(tests[rid])), "; ".join(sorted(sources[rid]))))

    groups = Counter("regresyon" if fresh[rid] else "arşiv" if archive[rid] else
                     "test atfı" if tests[rid] else "kod atfı" for rid in ids)
    only_code = sorted(rid for rid in ids if not (fresh[rid] or archive[rid] or tests[rid]))
    steps = Counter(row["status"] for row in read_tsv(RUN / "regresyon.tsv"))
    print(f"91 benzersiz kural / 19 aile; envanter, veri paketi ve iki dilli kılavuz eşleşiyor.")
    print("Kural başına en güçlü bulunan kayıt:", dict(groups))
    print("Yalnız kod/profil atfı:", ", ".join(only_code))
    print("Regresyon:", dict(steps))


if __name__ == "__main__":
    main()
