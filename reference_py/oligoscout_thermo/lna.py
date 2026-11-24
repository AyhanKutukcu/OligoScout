"""LNA (locked nucleic acid) katkisi: McTigue, Peterson & Kahn 2004.

Model (makale, denklem 12 ve "Proposed LNA Incorporation Algorithm"):
    dH(LNA+DNA:DNA) = dH(DNA:DNA) + sum_LNA [ddH(5'-M X^L) + ddH(5'-X^L N)]
    dS(LNA+DNA:DNA) = dS(DNA:DNA) + sum_LNA [ddS(5'-M X^L) + ddS(5'-X^L N)]
DNA dupleksi birlesik NN parametreleriyle (SantaLucia 1998) hesaplanir; tuz
ve konsantrasyon DNA ile ayni yoldan uygulanir (makale: tuz etkisi dS'e
yedirilir). Tm denklemi degismez.

Gecerlilik alani (makalenin olcum tasarimi): tek LNA, iki ucdan en az 2 nt
iceride, tam eslesen DNA tamamlayicisi. Alan disindaki durumlar
(uca yakin LNA, bitisik LNA'lar, LNA'ya komsu uyumsuzluk, birden cok LNA)
yine hesaplanir; ancak `in_domain=False` ve gerekceleriyle isaretlenir.
Birden cok aralikli LNA icin katkilar toplanir (makale: "should apply to
multiple incorporations if spaced far enough apart").
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .nn import DuplexThermo
from .params import LNA_MCTIGUE_2004

TABLE_ID = LNA_MCTIGUE_2004.provenance.table_id
MIN_END_DISTANCE = 2   # makalede LNA her zaman uclardan en az 2 nt iceride


@dataclass(frozen=True)
class LnaIncrement:
    ddH_kcal: float
    ddS_cal: float
    ddG37_kcal: float
    per_position: tuple[dict, ...]
    in_domain: bool
    out_of_domain: tuple[str, ...] = field(default_factory=tuple)

    def as_dict(self) -> dict:
        return {
            "model": TABLE_ID,
            "model_status": LNA_MCTIGUE_2004.provenance.verification_status,
            "lna_count": len(self.per_position),
            "ddH_kcal_mol": round(self.ddH_kcal, 4),
            "ddS_cal_mol_K": round(self.ddS_cal, 4),
            "ddG37_kcal_mol": round(self.ddG37_kcal, 4),
            "per_position": list(self.per_position),
            "in_domain": self.in_domain,
            "out_of_domain": list(self.out_of_domain),
        }


def increments(dna: str, positions, template_site: str | None = None) -> LnaIncrement:
    """LNA konumlari icin toplam ddH/ddS.

    `dna` 5'->3' DNA iskeleti; `positions` 0-tabanli LNA konumlari.
    `template_site` verilirse (primer ile ayni yonde, ayni uzunlukta hedef
    dizisi) LNA'nin kendisi ya da komsusu uyumsuz oldugunda o dinukleotidin
    artisi uygulanmaz ve alan disi olarak isaretlenir.
    """
    table = LNA_MCTIGUE_2004.require()
    t5, t3 = table["lna_5prime"], table["lna_3prime"]
    s = dna.upper()
    n = len(s)
    pos = sorted(set(int(p) for p in positions))
    for p in pos:
        if not 0 <= p < n:
            raise ValueError(f"LNA konumu dizinin disinda: {p}")
    mism = set()
    if template_site is not None:
        site = template_site.upper()
        if len(site) != n:
            raise ValueError("template_site uzunlugu dizi ile ayni olmali")
        mism = {i for i, (a, b) in enumerate(zip(s, site)) if a != b}
    lset = set(pos)
    reasons: list[str] = []
    if len(pos) > 1:
        reasons.append("birden cok LNA (katkilar toplanarak; makale tek LNA olctu)")
    dH = dS = dG = 0.0
    rows = []
    for p in pos:
        row = {"position": p, "base": s[p], "terms": []}
        if p < MIN_END_DISTANCE or p > n - 1 - MIN_END_DISTANCE:
            reasons.append(f"LNA {p}: uca {MIN_END_DISTANCE} nt'den yakin")
        if (p - 1) in lset or (p + 1) in lset:
            reasons.append(f"LNA {p}: bitisik LNA (LNA:LNA adimi olculmedi)")
        if p in mism or (p - 1) in mism or (p + 1) in mism:
            reasons.append(f"LNA {p}: uyumsuzluk LNA'da ya da komsusunda (model tam eslesme icin)")
        terms = []
        if p > 0:                       # 5'-M X^L-3'
            terms.append((s[p - 1:p + 1], t3, f"{s[p - 1]}{s[p]}L", {p - 1, p}))
        if p < n - 1:                   # 5'-X^L N-3'
            terms.append((s[p:p + 2], t5, f"{s[p]}L{s[p + 1]}", {p, p + 1}))
        for key, tab, label, members in terms:
            if members & mism:
                row["terms"].append({"step": label, "applied": False})
                continue
            h, sv, g = tab[key]
            dH += h
            dS += sv
            dG += g
            row["terms"].append({"step": label, "ddH": h, "ddS": sv, "ddG37": g, "applied": True})
        rows.append(row)
    # tekrar eden gerekceleri sirasini koruyarak birak
    seen: list[str] = []
    for r in reasons:
        if r not in seen:
            seen.append(r)
    return LnaIncrement(dH, dS, dG, tuple(rows), not seen, tuple(seen))


def apply(duplex: DuplexThermo, inc: LnaIncrement) -> DuplexThermo:
    """DNA dupleksine LNA artisini ekler."""
    return DuplexThermo(dH_kcal=duplex.dH_kcal + inc.ddH_kcal,
                        dS_cal=duplex.dS_cal + inc.ddS_cal,
                        self_complementary=duplex.self_complementary,
                        n_bp=duplex.n_bp)
