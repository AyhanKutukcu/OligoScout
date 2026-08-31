#!/usr/bin/env python3
"""Generate GitHub-safe SVG assets for the OligoScout README.

GitHub sanitises SVG embedded in Markdown, so these files use presentation
attributes only: no <style> blocks, no external fonts, no scripts. Light and
dark variants are emitted so the README can switch them with <picture>.
"""
import os

FONT = ("ui-sans-serif, -apple-system, BlinkMacSystemFont, 'Segoe UI', "
        "Helvetica, Arial, sans-serif")
MONO = ("ui-monospace, SFMono-Regular, 'SF Mono', Menlo, Consolas, "
        "'Liberation Mono', monospace")

LIGHT = dict(
    bg="#ffffff", panel="#f6f8fa", line="#d0d7de", ink="#1f2328",
    muted="#59636e", brand="#1f4e79", brandsoft="#b9cde2",
    exh="#2e7d5b", pcr="#7b4ea3", stat="#c2701c", heur="#a32e38",
    grid="#e4e8ec",
)
DARK = dict(
    bg="#0d1117", panel="#151b23", line="#2f3742", ink="#e6edf3",
    muted="#9198a1", brand="#5aa2e0", brandsoft="#1e3a55",
    exh="#4cbf8b", pcr="#b28ad6", stat="#e0a24a", heur="#e06c78",
    grid="#232a33",
)


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def text(x, y, s, size=14, fill="#000", weight="400", anchor="start",
         family=FONT, opacity=None, spacing=None):
    a = (f'<text x="{x}" y="{y}" font-family="{family}" font-size="{size}" '
         f'fill="{fill}" font-weight="{weight}" text-anchor="{anchor}"')
    if opacity is not None:
        a += f' opacity="{opacity}"'
    if spacing is not None:
        a += f' letter-spacing="{spacing}"'
    return a + f'>{esc(s)}</text>'


def rect(x, y, w, h, fill, rx=0, stroke=None, sw=1, opacity=None):
    a = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}"'
    if stroke:
        a += f' stroke="{stroke}" stroke-width="{sw}"'
    if opacity is not None:
        a += f' opacity="{opacity}"'
    return a + '/>'


# =====================================================================
# Hero banner
# =====================================================================
def hero(C):
    W, H = 1200, 300
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
         f'width="{W}" height="{H}" role="img" '
         f'aria-label="OligoScout: anchor-aware exhaustive search for primer '
         f'binding sites and off-target amplicons">']
    o.append(rect(0, 0, W, H, C["bg"]))
    o.append(rect(0, 0, 6, H, C["brand"]))

    # ---- wordmark and tagline
    o.append(text(54, 92, "OligoScout", 52, C["ink"], "700", spacing="-0.5"))
    o.append(text(54, 126, "ANCHOR-AWARE EXHAUSTIVE SEARCH FOR PRIMER PAIRS",
                  13, C["brand"], "600", spacing="1.6"))
    o.append(text(54, 162,
                  "Every binding site and every amplifiable product across the",
                  16, C["muted"]))
    o.append(text(54, 184,
                  "human reference — under one explicit molecular contract.",
                  16, C["muted"]))

    # ---- contract chips
    chips = ["24-nt primers", "exact 12-nt 3′ anchor",
             "≤ 3 mismatches", "50–3,000 bp products"]
    x = 54
    for c in chips:
        w = 11 + len(c) * 7.4
        o.append(rect(x, 210, w, 28, C["panel"], rx=14, stroke=C["line"]))
        o.append(text(x + w / 2, 229, c, 12.5, C["muted"], "500",
                      anchor="middle", family=MONO))
        x += w + 10

    # ---- schematic: inward-facing pair on a reference line
    gx, gy = 690, 132
    gw = 452

    o.append(rect(gx, gy, gw, 11, C["panel"], stroke=C["line"]))
    o.append(rect(gx, gy + 26, gw, 11, C["panel"], stroke=C["line"]))
    o.append(text(gx - 12, gy + 9, "+", 13, C["muted"], "600", anchor="middle"))
    o.append(text(gx - 12, gy + 35, "−", 13, C["muted"], "600", anchor="middle"))

    cw, ch = 7.4, 17

    def primer(x0, y0, flip, mism):
        parts = []
        for i in range(24):
            anchor = (i >= 12) if not flip else (i < 12)
            parts.append(rect(x0 + i * cw, y0, cw - 1.1, ch,
                              C["brand"] if anchor else C["brandsoft"], rx=1))
        for m in mism:
            idx = m if not flip else 23 - m
            cx = x0 + idx * cw + (cw - 1.1) / 2
            parts.append(f'<path d="M{cx - 3.4},{y0 - 8} L{cx + 3.4},{y0 - 8} '
                         f'L{cx},{y0 - 2.4} Z" fill="{C["heur"]}"/>')
        return parts

    o += primer(gx + 8, gy - 30, False, (1, 6))
    o.append(text(gx + 8, gy - 36, "5′", 10.5, C["muted"]))
    o.append(text(gx + 8 + 24 * cw, gy - 36, "3′", 10.5, C["muted"],
                  anchor="end"))
    o.append(f'<path d="M{gx + 12 + 24 * cw},{gy + 5.5} h44" stroke="{C["brand"]}" '
             f'stroke-width="1.6" marker-end="url(#ah)"/>')

    rx0 = gx + gw - 8 - 24 * cw
    o += primer(rx0, gy + 50, True, (0, 9))
    o.append(text(rx0, gy + 84, "3′", 10.5, C["muted"]))
    o.append(text(rx0 + 24 * cw, gy + 84, "5′", 10.5, C["muted"], anchor="end"))
    o.append(f'<path d="M{rx0 - 12},{gy + 31.5} h-44" stroke="{C["brand"]}" '
             f'stroke-width="1.6" marker-end="url(#ah)"/>')

    o.append(f'<path d="M{gx + 8},{gy - 52} H{gx + gw - 8}" stroke="{C["ink"]}" '
             f'stroke-width="1.1" marker-start="url(#ah2)" marker-end="url(#ah2)" '
             f'opacity="0.75"/>')
    o.append(text(gx + gw / 2, gy - 58, "amplicon", 11.5, C["ink"], "500",
                  anchor="middle", opacity="0.8"))

    o.append(rect(gx + 8, gy + 96, 13, 9, C["brandsoft"], rx=1))
    o.append(rect(gx + 23, gy + 96, 13, 9, C["brand"], rx=1))
    o.append(text(gx + 42, gy + 104, "5′ region / exact anchor", 11, C["muted"]))
    cxl = gx + 214
    o.append(f'<path d="M{cxl - 3.4},{gy + 96} L{cxl + 3.4},{gy + 96} '
             f'L{cxl},{gy + 105} Z" fill="{C["heur"]}"/>')
    o.append(text(cxl + 12, gy + 104, "mismatch (never inside the anchor)", 11,
                  C["muted"]))

    o.append(f'<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" '
             f'markerWidth="5" markerHeight="5" orient="auto-start-reverse">'
             f'<path d="M0,1 L9,5 L0,9 z" fill="{C["brand"]}"/></marker>'
             f'<marker id="ah2" viewBox="0 0 10 10" refX="8" refY="5" '
             f'markerWidth="5" markerHeight="5" orient="auto-start-reverse">'
             f'<path d="M0,1 L9,5 L0,9 z" fill="{C["ink"]}" opacity="0.75"/>'
             f'</marker></defs>')
    o.append('</svg>')
    return "\n".join(o)


# =====================================================================
# Benchmark overview
# =====================================================================
ROWS32 = [
    ("OligoScout", 65467, "brand", True),
    ("Bowtie 1", 65467, "exh", True),
    ("BWA-aln", 65467, "exh", True),
    ("SeqKit", 65467, "exh", True),
    ("BLAST+", 65454, "stat", False),
    ("minimap2", 6398, "heur", False),
    ("HISAT2", 5894, "heur", False),
]
UNIVERSE = 65467
T81 = [("OligoScout", 111.80, "brand"), ("BWA-aln", 141.14, "exh"),
       ("Bowtie 1", 148.64, "exh"), ("oracle", 194.11, "muted"),
       ("BLAST+", 3065.39, "stat")]


def benchmark(C):
    W, H = 1200, 520
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
         f'width="{W}" height="{H}" role="img" '
         f'aria-label="Amplicon recovery on the 32-pair panel and search wall '
         f'time on the 1,024-pair panel">']
    o.append(rect(0, 0, W, H, C["bg"]))

    # ---------- panel A
    o.append(text(40, 44, "Completed searches: 32-pair panel", 17,
                  C["ink"], "700"))
    o.append(text(40, 66, "fraction of the 65,467-amplicon exhaustive set "
                  "recovered under one shared contract", 12.5, C["muted"]))

    x0, y0, bw, rowh = 152, 92, 462, 34
    for gx in (0, 0.25, 0.5, 0.75, 1.0):
        px = x0 + gx * bw
        o.append(f'<path d="M{px},{y0 - 6} V{y0 + len(ROWS32) * rowh - 6}" '
                 f'stroke="{C["grid"]}" stroke-width="1"/>')
        o.append(text(px, y0 + len(ROWS32) * rowh + 12,
                      f"{gx:.2f}", 10.5, C["muted"], anchor="middle"))

    for i, (name, amp, col, eq) in enumerate(ROWS32):
        y = y0 + i * rowh
        frac = amp / UNIVERSE
        o.append(text(x0 - 12, y + 15, name, 12.5, C["ink"], "500",
                      anchor="end"))
        o.append(rect(x0, y, bw, 20, C["panel"], rx=3))
        o.append(rect(x0, y, max(2, bw * frac), 20, C[col], rx=3))
        lbl = f"{frac:.4f}"
        if frac > 0.55:
            o.append(text(x0 + bw * frac - 10, y + 15, lbl, 11.5,
                          "#ffffff", "700", anchor="end"))
        else:
            o.append(text(x0 + bw * frac + 9, y + 15, lbl, 11.5,
                          C["ink"], "600"))
        o.append(text(x0 + bw + 76, y + 15, f"{amp:,}", 11, C["muted"],
                      anchor="end", family=MONO))
        if eq:
            o.append(rect(x0 + bw + 86, y + 3, 46, 15, C["exh"], rx=7,
                          opacity=0.14))
            o.append(text(x0 + bw + 109, y + 14, "exact", 9.5, C["exh"],
                          "700", anchor="middle", spacing="0.4"))

    o.append(text(40, 386, "Bowtie 2: partial (31/32); details in README.",
                  12, C["muted"]))
    o.append(text(40, 407, "MFEprimer: different panel; comparison not verified.",
                  12, C["muted"]))

    # ---------- panel B
    bx = 800
    o.append(f'<path d="M{bx - 40},80 V470" stroke="{C["line"]}" stroke-width="1"/>')
    o.append(text(bx, 44, "Publication-scale panel", 17, C["ink"], "700"))
    o.append(text(bx, 66, "1,024 primer pairs vs. an independent oracle",
                  12.5, C["muted"]))

    cards = [("8,262,803", "binding sites"), ("425,055", "amplicons"),
             ("1.0000000000", "precision / recall / F1 / Jaccard")]
    cy = 88
    for big, small in cards:
        o.append(rect(bx, cy, 360, 52, C["panel"], rx=8, stroke=C["line"]))
        o.append(text(bx + 16, cy + 26, big, 20, C["brand"], "700",
                      family=MONO))
        o.append(text(bx + 16, cy + 43, small, 11.5, C["muted"]))
        cy += 60

    o.append(text(bx, cy + 22, "Reported search-stage wall time", 12.5,
                  C["ink"], "600"))
    ty0 = cy + 36
    tmax = max(v for _, v, _ in T81)
    for i, (name, val, col) in enumerate(T81):
        y = ty0 + i * 30
        o.append(text(bx + 74, y + 14, name, 11.5, C["ink"], "500",
                      anchor="end"))
        w = max(3, 160 * (val ** 0.5) / (tmax ** 0.5))
        o.append(rect(bx + 84, y, w, 18, C[col] if col != "muted" else C["muted"],
                      rx=3))
        o.append(text(bx + 84 + w + 8, y + 14, f"{val:,.2f} s", 11,
                      C["muted"], family=MONO))
    o.append(text(bx, ty0 + len(T81) * 30 + 16,
                  "bars ∝ √time; scope details in README",
                  10.5, C["muted"]))

    o.append(text(40, H - 16,
                  "GRCh38.p14 Primary 24  ·  24-nt primers  ·  exact 12-nt "
                  "3′ anchor  ·  ≤ 3 mismatches  ·  50–3,000 bp products",
                  11.5, C["muted"], family=MONO))
    o.append('</svg>')
    return "\n".join(o)


if __name__ == "__main__":
    os.makedirs("docs/assets", exist_ok=True)
    for name, fn in (("oligoscout-hero", hero), ("oligoscout-benchmark", benchmark)):
        for suffix, pal in (("", LIGHT), ("-dark", DARK)):
            path = f"docs/assets/{name}{suffix}.svg"
            open(path, "w", encoding="utf-8").write(fn(pal))
            print("wrote", path)
