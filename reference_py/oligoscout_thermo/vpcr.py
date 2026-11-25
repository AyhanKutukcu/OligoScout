"""Sanal (in silico) PCR: kullanici kaliplari, dejenere oligolar, dairesel kalip.

Kalendar et al. 2024 (Front. Bioinform. 4:1464197, virtualPCR) ile ayni
gorev kapsami:
  * herhangi bir FASTA (cok kayitli; dogrusal ya da dairesel) uzerinde arama
  * primer VE kalipta IUPAC dejenere bazlar (kume kesisimi eslesme sayilir)
  * 3' uc kurali: son `terminal_window` (7) bazda en fazla `terminal_max` (1)
    uyumsuzluk, toplamda en fazla `max_mismatch`
  * uca eklenmis tamamlayici olmayan diziler (kucuk harf = kuyruk)
  * urunler: tasarlanmis ciftler (F-R, R-F, F-F, R-R) ya da TUM oligolar
    birbiriyle (parmak izi: IRAP / ISSR / RAPD / SINE-PCR, tek primer dahil)
  * dairesel kalipta orijinden gecen urunler ve ters (inverse) PCR
  * prob arama, urun icindeki problar ve baglantili dizi arama
  * urun basina Ta = en dusuk primer-yer Tm'i + ln(urun uzunlugu)
  * amplikon dizisi (primer kuyruklari dahil)

Arama ya C++ `oligoscout_vpcr_scan` (hizli, eksiksiz) ya da bu modulun saf
Python tarayicisi (kucuk kaliplar, ayni kural; testlerde oracle) ile yapilir.
"""

from __future__ import annotations

import bisect
import collections
import math
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Iterable

from . import align, nn, salt
from .params import KELVIN, R_CAL

IUPAC_SETS = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
_MASK = {k: sum(1 << "ACGT".index(b) for b in v) for k, v in IUPAC_SETS.items()}
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
         "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N", "U": "A"}


def revcomp_iupac(s: str) -> str:
    return "".join(_COMP[c] for c in reversed(s.upper()))


# ---------------------------------------------------------------------------
# Oligo tanimi
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Oligo:
    id: str
    core: str                 # baglanan kisim, IUPAC, 5'->3'
    tail5: str = ""           # 5' kuyruk (baglanmaz; urune girer)
    tail3: str = ""           # 3' kuyruk (baglanmaz; primerde uzamayi engeller)
    role: str = "primer"      # primer | probe | link
    pair: str = ""            # tasarlanmis cift kimligi (pairs modu)
    lna: tuple[int, ...] = () # cekirdekteki LNA konumlari (yalniz termodinamik)

    @property
    def full(self) -> str:
        return self.tail5 + self.core + self.tail3

    @property
    def degeneracy(self) -> int:
        out = 1
        for c in self.core:
            out *= len(IUPAC_SETS[c])
        return out


_SPEC = re.compile(r"^([a-z]*)([A-Z+]+)([a-z]*)$")


def parse_oligo(oid: str, raw: str, role: str = "primer", pair: str = "",
                tail5: str = "", tail3: str = "") -> Oligo:
    """Oligo yazimini ayristirir.

    Kucuk harfli bas/son = baglanmayan kuyruk (ornek: `gtaaaacgacggccagtACGT..`).
    `+X` LNA isaretidir (termodinamik icin). Bosluk ve /uc etiketleri/ atilir.
    """
    s = re.sub(r"/[^/]*/", "", raw or "")
    s = re.sub(r"[\s\-_]", "", s)
    if not s:
        raise ValueError(f"{oid}: bos oligo")
    m = _SPEC.match(s)
    if not m:
        if any(c.islower() for c in s.replace("+", "")) and any(c.isupper() for c in s):
            raise ValueError(f"{oid}: kucuk harf yalniz 5'/3' kuyrukta olabilir")
        m = _SPEC.match(s.upper())
    t5, core, t3 = m.group(1), m.group(2), m.group(3)
    lna, plain = [], []
    nxt = False
    for c in core:
        if c == "+":
            nxt = True
            continue
        if nxt:
            lna.append(len(plain))
        plain.append(c.upper())
        nxt = False
    core = "".join(plain).replace("U", "T")
    bad = sorted({c for c in core if c not in IUPAC_SETS})
    if bad:
        raise ValueError(f"{oid}: tanimsiz baz {bad} (IUPAC: ACGT RYSWKM BDHV N)")
    tail5 = (tail5 or t5).upper()
    tail3 = (tail3 or t3).upper()
    if not 8 <= len(core) <= 512:
        raise ValueError(f"{oid}: baglanan kisim 8-512 nt olmali (su an {len(core)})")
    return Oligo(oid, core, tail5, tail3, role, pair, tuple(lna))


# ---------------------------------------------------------------------------
# Kaliplar
# ---------------------------------------------------------------------------
@dataclass
class Record:
    name: str
    seq: str
    circular: bool = False

    @property
    def length(self) -> int:
        return len(self.seq)


def read_fasta(path: str, circular: Iterable[str] | str = ()) -> list[Record]:
    circ_all = circular == "all"
    circ = set() if circ_all else set(circular)
    out: list[Record] = []
    name, parts = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    out.append(Record(name, "".join(parts).upper(), circ_all or name in circ))
                name = line[1:].split()[0] if line[1:].split() else f"seq{len(out) + 1}"
                parts = []
            elif name is not None:
                parts.append(line.strip())
    if name is not None:
        out.append(Record(name, "".join(parts).upper(), circ_all or name in circ))
    return out


# ---------------------------------------------------------------------------
# Baglanma yerleri
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Site:
    oligo: str
    record: str
    start: int                # arti iplikte 0-tabanli
    end: int                  # dislayici; dairesel kayitta kayit uzunlugunu asabilir
    strand: str               # "+" oligo arti iplikte okunur (sag uca dogru uzar)
    mismatches: int
    terminal_mismatches: int
    mask: str                 # oligonun 5' ucundan; X = uyumsuz, g = hos gorulen G-ciftli uyumsuzluk
    site: str                 # arti iplik dizisi
    wraps: bool = False
    g_mismatches: int = 0     # G.T / G.G / G.A (yalniz g_tolerant kuralinda sayilmaz)


@dataclass(frozen=True)
class ScanRule:
    max_mismatch: int = 3
    terminal_window: int = 7
    terminal_max: int = 1
    template_n: str = "reject"     # reject | mismatch | match
    g_tolerant: bool = False       # G.T, G.G, G.A uyumsuzluklari sinira sayilmaz ...
    g_max: int = 2                 # ... en fazla bu kadari (fazlasi siradan sayilir)


_CB = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _stable_g_pair(x: str, y: str) -> bool:
    """x.y uyumsuz cifti G.T, G.G ya da G.A mi (Kalendar et al. 2024)."""
    return (x == "G" and y in "AGT") or (y == "G" and x in "AGT")


def _g_tolerated(oligo_set: str, tmpl_base: str, strand: str) -> bool:
    """Uyumsuzluk G-ciftli mi. Pencere harfi kalibin arti ipligidir.

    Arti iplikteki oligo (pencere j'de oligo bazi x) comp(b) ile eslesir;
    eksi iplikteki oligo pencerede revcomp olarak yazilir: harf q icin oligo
    bazi comp(q) olup b ile eslesir.
    """
    if tmpl_base not in "ACGT":
        return False
    if strand == "+":
        return any(_stable_g_pair(x, _CB[tmpl_base]) for x in IUPAC_SETS[oligo_set])
    return any(_stable_g_pair(_CB[q], tmpl_base) for q in IUPAC_SETS[oligo_set])


def _site_from_window(o: Oligo, rec: Record, start: int, window: str, strand: str,
                      rule: ScanRule) -> Site | None:
    L = len(o.core)
    pat = o.core if strand == "+" else revcomp_iupac(o.core)
    mm = tm = gm = gt = 0
    mask = ["."] * L
    for j in range(L):
        t = window[j]
        if t == "N" and rule.template_n == "reject":
            return None
        tm_ = 15 if (t == "N" and rule.template_n == "match") else (0 if t == "N" else _MASK.get(t, 0))
        if not (tm_ & _MASK[pat[j]]):
            p = j if strand == "+" else L - 1 - j          # oligo konumu (5' kokenli)
            terminal = p >= L - rule.terminal_window
            if rule.g_tolerant and _g_tolerated(pat[j], t, strand):
                mask[p] = "g"
                gm += 1
                gt += terminal
                continue
            mask[p] = "X"
            mm += 1
            tm += terminal
            if mm > rule.max_mismatch or tm > rule.terminal_max:
                return None
    if rule.g_tolerant and gm > rule.g_max:
        tol_t = min(gt, rule.g_max)
        tol_o = min(gm - gt, rule.g_max - tol_t)
        mm += (gt - tol_t) + (gm - gt - tol_o)
        tm += gt - tol_t
        if mm > rule.max_mismatch or tm > rule.terminal_max:
            return None
    total = sum(1 for c in mask if c != ".")
    return Site(o.id, rec.name, start, start + L, strand, total, tm, "".join(mask), window,
                start + L > rec.length, gm)


def scan_python(records: list[Record], oligos: list[Oligo], rule: ScanRule) -> list[Site]:
    """Saf Python tarama: her konum, iki iplik, her oligo (kucuk kaliplar icin)."""
    out: list[Site] = []
    for rec in records:
        n = rec.length
        maxlen = max(len(o.core) for o in oligos)
        work = rec.seq + (rec.seq[: maxlen - 1] if rec.circular else "")
        for o in oligos:
            L = len(o.core)
            for i in range(n):
                if i + L > len(work):
                    break
                win = work[i:i + L]
                for strand in "+-":
                    s = _site_from_window(o, rec, i, win, strand, rule)
                    if s:
                        out.append(s)
    return out


def default_scan_binary() -> str | None:
    env = os.environ.get("OLIGOSCOUT_VPCR_BINARY")
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    for p in (env, os.path.join(root, "OligoScout", "build", "release", "oligoscout_vpcr_scan")):
        if p and os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return None


def scan_cpp(fasta: str, oligos: list[Oligo], rule: ScanRule, circular: Iterable[str] | str = (),
             binary: str | None = None, threads: int = 2, primary_only: bool = False
             ) -> tuple[list[Site], dict[str, tuple[int, bool]]]:
    binary = binary or default_scan_binary()
    if not binary:
        raise FileNotFoundError("oligoscout_vpcr_scan bulunamadi (OLIGOSCOUT_VPCR_BINARY)")
    with tempfile.TemporaryDirectory(prefix="oligoscout-vpcr-") as tmp:
        ol = os.path.join(tmp, "oligos.tsv")
        out = os.path.join(tmp, "sites.tsv")
        with open(ol, "w", encoding="utf-8") as fh:
            fh.write("oligo_id\tsequence\n")
            for o in oligos:
                fh.write(f"{o.id}\t{o.core}\n")
        cmd = [binary, fasta, ol, out, "--max-mismatch", str(rule.max_mismatch),
               "--terminal-window", str(rule.terminal_window), "--terminal-max", str(rule.terminal_max),
               "--template-n", rule.template_n, "--threads", str(threads)]
        if rule.g_tolerant:
            cmd += ["--g-tolerant", str(rule.g_max)]
        if circular == "all":
            cmd += ["--circular", "all"]
        elif circular:
            cmd += ["--circular", ",".join(circular)]
        if primary_only:
            cmd.append("--primary-only")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        sites = []
        with open(out, encoding="utf-8") as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                sites.append(Site(f[0], f[1], int(f[2]), int(f[3]), f[4], int(f[5]), int(f[6]),
                                  f[7], f[8], f[9] == "1", int(f[10])))
        recs = {}
        with open(out + ".records", encoding="utf-8") as fh:
            next(fh)
            for line in fh:
                name, length, topo = line.rstrip("\n").split("\t")
                recs[name] = (int(length), topo == "circular")
    return sites, recs


# ---------------------------------------------------------------------------
# Termodinamik: yer Tm'i ve Ta
# ---------------------------------------------------------------------------
def resolve_core(core: str, site_sense: str) -> str:
    """Dejenere primeri, yere uyan bazlarla cozer (uyumsuz konumda ilk baz)."""
    out = []
    for p, t in zip(core, site_sense):
        allowed = IUPAC_SETS[p]
        out.append(t if t in allowed else allowed[0])
    return "".join(out)


def site_sense(site: Site) -> str:
    """Yer dizisi, oligo ile ayni yonde (tam eslesmede == oligo)."""
    return site.site if site.strand == "+" else revcomp_iupac(site.site)


def mismatched_duplex(primer: str, template: str) -> tuple[float, float, int, str, list[str]]:
    """Birden cok ic uyumsuzluklu dupleks (Allawi & SantaLucia NN, toplanir).

    `template` primerle ayni yonde yazilmis yer dizisidir (tam eslesmede ==).
    Uclara 2 nt'den yakin uyumsuzluklar, ucu o uyumsuzlukla birlikte kirparak
    cikarilir (o kisim eslesmemis sarkan uc sayilir). Kalan ic uyumsuzluklar
    yalitilmissa (aralarinda en az bir WC cifti) Allawi parametreleri
    toplanir. Bitisik (tandem) uyumsuzluk icin parametre yoktur: None doner.
    Donus: dH, dS, cift sayisi, eslesen cekirdek, notlar.
    """
    lo, hi = 0, len(primer)
    notes: list[str] = []
    while True:
        mm = [i for i in range(lo, hi) if primer[i] != template[i]]
        edge = [i for i in mm if i < lo + 2 or i >= hi - 2]
        if not edge:
            break
        for i in edge:
            if i < lo + 2:
                lo = max(lo, i + 1)
            else:
                hi = min(hi, i)
        if hi - lo < 6:
            raise ValueError("eslesen cekirdek 6 nt'den kisa")
    if lo or hi != len(primer):
        notes.append(f"uc uyumsuzluklari kirpildi: cekirdek {lo}-{hi}")
    p, t = primer[lo:hi], template[lo:hi]
    mm = [i for i in range(len(p)) if p[i] != t[i]]
    if any(b - a == 1 for a, b in zip(mm, mm[1:])):
        raise ValueError("bitisik (tandem) uyumsuzluk")
    base = nn.watson_crick(p)
    h, sv = base.dH_kcal, base.dS_cal
    bottom = t.translate(str.maketrans("ACGT", "TGCA"))
    for i in mm:
        dh, ds = nn._mismatch_stack_adjustment(p, bottom, i)
        h += dh
        sv += ds
    if mm:
        notes.append(f"{len(mm)} ic uyumsuzluk (Allawi-Peyret NN)")
    return h, sv, len(p), p, notes


def site_tm(o: Oligo, site: Site, cond: salt.Conditions, salt_model: str = "owczarzy_2008") -> dict:
    """Oligo ile kalip yeri arasindaki dupleksin Tm'i (primer fazlasi; tuz duzeltmeli).

    Tam eslesme: WC NN. Yalitilmis ic uyumsuzluklar: Allawi-Peyret NN toplami.
    Uca yakin uyumsuzluklar kirpilir. Tandem uyumsuzlukta en kararli tek
    kusurlu sarmal (align.best_dimer) kullanilir ve not edilir.
    """
    sense = "".join(c if c in "ACGT" else "N" for c in site_sense(site))
    primer = resolve_core(o.core, sense)
    template = "".join(c if c in "ACGT" else primer[i] for i, c in enumerate(sense))
    model = "nn_multi_mismatch"
    try:
        h, s, n_bp, core, notes = mismatched_duplex(primer, template)
        gc = nn.gc_fraction(core)
        three = primer[-1] == template[-1]
    except ValueError as exc:
        d = align.best_dimer(primer, nn.reverse_complement(template), "any")
        if not d.found or d.helix_len < 2:
            return {"tm_c": None, "dG37_kcal_mol": None, "helix_len": 0, "three_prime_paired": False,
                    "model": "none", "notes": [str(exc)]}
        h, s, n_bp = d.dH_kcal, d.dS_cal, d.helix_len
        seg = d.top.replace("-", "").strip()
        gc = nn.gc_fraction(seg) if seg else 0.5
        three = d.three_prime_paired
        model, notes = "best_single_defect_helix", [str(exc)]
    denom = s + R_CAL * math.log(cond.primer_M)
    out = {"dG37_kcal_mol": round(h - 310.15 * s / 1000, 3), "helix_len": n_bp,
           "three_prime_paired": three, "model": model, "notes": notes}
    if denom >= 0:
        return {**out, "tm_c": None}
    tm1 = h * 1000.0 / denom
    tmK, _ = salt.correct(tm1, gc, max(n_bp, 2), cond, salt_model)
    return {**out, "tm_c": round(tmK - KELVIN, 2)}


def annealing_temperature(tm_left: float | None, tm_right: float | None, length: int) -> float | None:
    """Kalendar et al. 2024: Ta = en dusuk primer Tm'i + ln(urun uzunlugu)."""
    tms = [t for t in (tm_left, tm_right) if t is not None]
    if not tms or length <= 0:
        return None
    return round(min(tms) + math.log(length), 2)


# ---------------------------------------------------------------------------
# Urunler
# ---------------------------------------------------------------------------
@dataclass
class Product:
    record: str
    circular: bool
    left: Site                 # arti iplikte uzayan (sag uca dogru)
    right: Site                # eksi iplikte uzayan
    start: int
    end: int                   # dislayici; dairesel urunde kayit uzunlugunu asabilir
    template_length: int       # kalipta kapsanan uzunluk (kuyruklar haric)
    length: int                # kuyruklar dahil amplikon uzunlugu
    orientation: str           # ornek: "F-R", "RT1-RT1"
    wraps: bool
    record_length: int = 0
    extra: dict = field(default_factory=dict)


def _allowed(left: Oligo, right: Oligo, mode: str) -> tuple[bool, str]:
    if left.role != "primer" or right.role != "primer":
        return False, ""
    if left.tail3 or right.tail3:
        return False, ""           # 3' kuyruklu oligo uzayamaz
    if mode == "all":
        return True, f"{left.id}|{right.id}"
    if left.pair and left.pair == right.pair:
        return True, f"{left.id}|{right.id}"
    return False, ""


def assemble(sites: list[Site], oligos: dict[str, Oligo], records: dict[str, tuple[int, bool]],
             min_len: int = 50, max_len: int = 5000, mode: str = "pairs") -> list[Product]:
    """Iki yerden urun kurar.

    Sol yer arti iplikte (oligo okunusu sag uca dogru), sag yer eksi iplikte;
    sag yerin baslangici sol yerin sonundan once olamaz (ortusme yok).
    Uzunluk kalip uzerindeki kapsama + iki 5' kuyruktur ve [min, max]
    araliginda olmalidir. Dairesel kayitta mesafe halka uzerinde olculur;
    boylece orijinden gecen urunler ve ters (inverse) PCR urunleri de bulunur.
    """
    by_rec: dict[str, tuple[list[Site], list[Site]]] = collections.defaultdict(lambda: ([], []))
    for s in sites:
        o = oligos.get(s.oligo)
        if o is None or o.role != "primer":
            continue
        by_rec[s.record][0 if s.strand == "+" else 1].append(s)
    products: list[Product] = []
    for rec, (plus, minus) in by_rec.items():
        n, circ = records.get(rec, (None, False))
        minus.sort(key=lambda s: s.start)
        starts = [s.start for s in minus]
        for p in plus:
            lo_ = oligos[p.oligo]
            span_lo = p.start + len(lo_.core)                     # sag yerin en erken baslangici
            span_hi = p.start + max_len                           # kabaca ust sinir
            windows = [(span_lo, span_hi, 0)]
            if circ and n:
                windows = [(span_lo, span_hi, 0), (span_lo - n, span_hi - n, n)]
            seen = set()
            for a, b, shift in windows:
                i = bisect.bisect_left(starts, a)
                while i < len(minus) and starts[i] <= b:
                    m = minus[i]
                    i += 1
                    if (m.start, m.oligo, shift) in seen:
                        continue
                    seen.add((m.start, m.oligo, shift))
                    ro = oligos[m.oligo]
                    ok, label = _allowed(lo_, ro, mode)
                    if not ok:
                        continue
                    m_start = m.start + shift
                    if m_start < p.start + len(lo_.core):
                        continue
                    tlen = m_start + len(ro.core) - p.start
                    if circ and n and tlen > n + len(ro.core):
                        continue
                    total = tlen + len(lo_.tail5) + len(ro.tail5)
                    if not min_len <= total <= max_len:
                        continue
                    products.append(Product(rec, circ, p, m, p.start, p.start + tlen, tlen, total,
                                            label, bool(circ and n and p.start + tlen > n), n or 0))
    products.sort(key=lambda x: (x.record, x.start, x.length, x.orientation))
    return products


def template_slice(rec: Record, start: int, end: int) -> str:
    n = rec.length
    if end <= n:
        return rec.seq[start:end]
    return rec.seq[start:] + rec.seq[: end - n]


def amplicon_sequence(prod: Product, rec: Record, oligos: dict[str, Oligo]) -> dict:
    """Kalip dizisi ve gercek amplikon (primer dizileri + kuyruklar)."""
    lo_, ro = oligos[prod.left.oligo], oligos[prod.right.oligo]
    tmpl = template_slice(rec, prod.start, prod.end)
    inner = tmpl[len(lo_.core): len(tmpl) - len(ro.core)]
    amp = lo_.tail5 + lo_.core + inner + revcomp_iupac(ro.core) + revcomp_iupac(ro.tail5)
    return {"template_sequence": tmpl, "amplicon_sequence": amp}


def probes_inside(prod: Product, probe_sites: list[Site]) -> list[str]:
    out = []
    for s in probe_sites:
        if s.record != prod.record:
            continue
        st = s.start
        if prod.wraps and st < prod.start:
            st += prod.record_length          # halkanin ikinci turu
        if prod.start <= st and st + (s.end - s.start) <= prod.end:
            out.append(f"{s.oligo}{s.strand}@{s.start}")
    return out


@dataclass
class Link:
    record: str
    a: Site
    b: Site
    strand: str
    gap: int
    insert_start: int
    insert_end: int


def find_links(sites: list[Site], a_id: str, b_id: str, min_gap: int, max_gap: int) -> list[Link]:
    """Iki dizinin ayni iplikte, A once B sonra, aralarinda [min,max] baz ile yer almasi.

    Ornek: bir plazmit ekinin iki yanindaki diziler (Kalendar et al. 2024,
    "linked sequences" / "probe search").
    """
    A = [s for s in sites if s.oligo == a_id]
    B = collections.defaultdict(list)
    for s in sites:
        if s.oligo == b_id:
            B[(s.record, s.strand)].append(s)
    out = []
    for a in A:
        for b in B.get((a.record, a.strand), ()):
            if a.strand == "+":
                gap = b.start - a.end
                ins = (a.end, b.start)
            else:
                gap = a.start - b.end
                ins = (b.end, a.start)
            if min_gap <= gap <= max_gap:
                out.append(Link(a.record, a, b, a.strand, gap, ins[0], ins[1]))
    return out
