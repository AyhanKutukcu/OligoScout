
"""Build the self-contained OligoScout page from the UI and bundled evidence."""
import base64
import csv
import json
from pathlib import Path

WEB = Path(__file__).resolve().parent
SOURCE = WEB / "src"


def read(name):
    return (SOURCE / name).read_text(encoding="utf-8")


def asset_url(name, mime):
    encoded = base64.b64encode((SOURCE / "assets" / name).read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def rule_guide(catalog):
    fields = ("name", "check", "why")
    guide, translations = {}, {}
    rows = csv.DictReader(read("qc_rules.tsv").splitlines(), delimiter="\t")
    for number, row in enumerate(rows, 2):
        rid = row["id"]
        if rid in guide or None in row or not all(row.get(f"{key}_{lang}", "").strip()
                                                 for key in fields for lang in ("tr", "en")):
            raise ValueError(f"qc_rules.tsv:{number}: duplicate ID or incomplete translation: {rid}")
        guide[rid] = {key: row[f"{key}_tr"] for key in fields}
        for key in fields:
            translations[row[f"{key}_tr"]] = row[f"{key}_en"]
    ids = [rule["id"] for rule in catalog]
    if len(ids) != len(set(ids)) or set(guide) != set(ids):
        raise ValueError(f"QC catalog and guide IDs differ: missing {set(ids) - guide.keys()}, "
                         f"extra {guide.keys() - set(ids)}")
    return guide, translations


def script_json(value):
    return json.dumps(value, ensure_ascii=False).replace("</", r"<\/")


def english_catalog(guide_translations):
    messages = {}
    for number, line in enumerate(read("english.tsv").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        source, separator, english = line.partition("\t")
        if not separator or not source or not english:
            raise ValueError(f"english.tsv:{number}: expected a source and a translation")
        if source in messages:
            raise ValueError(f"english.tsv:{number}: duplicate message {source!r}")
        messages[source] = english
    return script_json(messages | guide_translations)


shell = read("shell.html")
for theme in ("light", "dark"):
    shell = shell.replace(
        f"__LOGO_{theme.upper()}__", asset_url(f"oligoscout-logo-{theme}.png", "image/png")
    )

css = read("styles.css")
for name in ("geist", "instrument-sans"):
    for subset in ("latin", "latin-ext"):
        token = f"__{name}_{subset}__".upper().replace("-", "_")
        css = css.replace(token, asset_url(f"{name}-{subset}.woff2", "font/woff2"))

scripts = (
    "engine.js", "i18n.js", "core.js", "viz.js", "views1.js", "views2.js",
    "views3.js", "catalog.js", "live.js", "multiplex.js", "docs.js", "init.js",
)
raw = (WEB / "oligoscout_data.json").read_text(encoding="utf-8")
guide, translations = rule_guide(json.loads(raw)["catalog"])
raw = raw.replace("</", r"<\/")
js = '"use strict";\nconst ENGLISH_MESSAGES = ' + english_catalog(translations) + ';\n'
js += 'const QC_RULE_GUIDE = ' + script_json(guide) + ';\n'
js += "\n".join(read(name) for name in scripts)
js = js.replace("__OVERVIEW_ART__", asset_url("overview-helix.png", "image/png"))
html = (
    '<!doctype html>\n<html lang="tr">\n<head>\n<meta charset="utf-8">\n'
    '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
    '<title>OligoScout</title>\n<style>' + css + '</style>\n</head>\n<body>\n'
    + shell + '\n<script type="application/json" id="osdata">' + raw + '</script>\n'
    '<script>\n' + js + '\n</script>\n</body>\n</html>\n'
)
output = WEB / "index.html"
output.write_text(html, encoding="utf-8")
print("web/index.html", output.stat().st_size, "bytes")
