"""Parametre tablolarinin bagimsiz kaynaklara karsi diferansiyel denetimi.

NE YAPAR. `params.py` icindeki tablolari, ayni birincil kaynaklardan
BAGIMSIZ olarak aktarilmis iki kaynakla karsilastirir:

  * `primer3_thal`    -- primer3'un `primer3_config/` dosyalari
                         (stack, stackmm, dangle, loops).
  * `primer3_oligotm` -- primer3'un `src/oligotm.c` dosyasindaki
                         SantaLucia_1998_dH/dS sabitleri, baslatma ve
                         simetri terimleri.

NE YAPMAZ. Bu modul hicbir degeri `params.py`'ye KOPYALAMAZ. Projenin
kurali degismeden durur: parametreler birincil makalelerden aktarilir,
baska yazilimlarin parametre dosyalari urune alinmaz. Buradaki kaynaklar
yalnizca ORACLE'dir; `tests/test_thermo.py` icindeki Biopython oracle'i ile
ayni roldedir ve urun hesaplarina girmez.

NEYIN YERINE GECMEZ. Bir tablonun `verification_status` alani ancak
birincil PDF'ten IKINCI BIR KISI tarafindan bagimsiz aktarim ve diff ile
"verified" olur. Diferansiyel denetim transkripsiyon hatasini yakalar,
o insan adiminin yerine gecmez. Sonuc, kunyedeki ayri bir `cross_check`
alanina yazilir.

ANAHTAR SOZLESMELERI (kaynaktan dogrulanmistir, bkz. thal_params.c):
  * Baz sirasi A=0, C=1, G=2, T=3 (`str2int`).
  * stack / stackmm: `XY_ZW` -> ust iplik 5'-XY-3', alt iplik 3'-ZW-5'.
    `params.py` uyumsuzluk anahtari "XY/ZW" ile birebir ayni sozlesme.
  * dangle 3' bolumu: `AK_J` -> i=A (sarkan iplikte eslesen baz),
    k=K (sarkan baz), j=J (karsi baz) -> params anahtari "D3:K:A-J".
  * dangle 5' bolumu: `I_JK` -> params anahtari "D5:K:J-I".
  * loops: satir = uzunluk, sutunlar INTERNAL, BULGE, HAIRPIN.
  * Birimler: primer3 dosyalarinda dH cal/mol, dS cal/(mol K);
    oligotm.c'de dH = -100*v cal/mol, dS = -0.1*v cal/(mol K).

KAPSAM DISI. `salt.py` katsayilari (Owczarzy 2004/2008, von Ahsen 2001)
primer3'te ayni biciminde bulunmaz; denetlenmez. `params.py`'deki
baslatma ve simetri terimleri yalniz `primer3_oligotm` ile denetlenir:
primer3'un thal motoru farkli bir baslatma sozlesmesi kullanir
(uc basina AT cezasi 2.2 kcal / 6.9 cal), bu yuzden thal bu karsilastirmaya
katilmaz.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, replace

from .params import (
    DANGLING_BOMMARITO_2000,
    DNA_LOOP_SANTALUCIA_HICKS_2004,
    MISMATCH_ALLAWI_PEYRET,
    TERMINAL_MISMATCH_SANTALUCIA_PEYRET,
    WC_SANTALUCIA_1998,
)

T37 = 310.15
BASES = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
WC_PAIRS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}

# Aktarim hatasi esikleri. Her iki kaynak da ayni ondalik cozunurlukte
# yayinlandigi icin gercek fark sifir olmalidir; tolerans yalnizca
# yuvarlama ve dG37 turetmesi icindir.
TOL_DH_KCAL = 0.051
TOL_DS_CAL = 0.051
TOL_DG37_KCAL = 0.051

ORACLE_THAL = "primer3_thal"
ORACLE_OLIGOTM = "primer3_oligotm"
ORACLE_TERMINAL = "primer3_tstack2"


class OracleUnavailable(RuntimeError):
    """Bagimsiz kaynak bulunamadi ya da beklenen bicimde degil."""


# ---------------------------------------------------------------------------
# Kaynak bulma
# ---------------------------------------------------------------------------
def default_primer3_src() -> str | None:
    """primer3 kaynak klasorunu bul: once ortam degiskeni, sonra depo yolu."""
    env = os.environ.get("OLIGOSCOUT_PRIMER3_SRC")
    if env:
        return env if os.path.isdir(env) else None
    here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(here))
    guess = os.path.join(
        repo_root, "OligoScout", "third_party", "primer3", "src"
    )
    return guess if os.path.isdir(guess) else None


def oracle_sources(primer3_src: str | None = None) -> dict[str, str]:
    """Mevcut oracle dosyalarinin yollari. Eksikse anahtar bulunmaz."""
    src = primer3_src or default_primer3_src()
    if not src:
        return {}
    cfg = os.path.join(src, "primer3_config")
    found: dict[str, str] = {}
    needed = (
        "stack.dh", "stack.ds", "stackmm.dh", "stackmm.ds",
        "dangle.dh", "dangle.ds", "loops.dh", "loops.ds",
    )
    if all(os.path.isfile(os.path.join(cfg, f)) for f in needed):
        found["config_dir"] = cfg
    if all(os.path.isfile(os.path.join(cfg, f)) for f in ("tstack2.dh", "tstack2.ds")):
        found["terminal_dir"] = cfg
    oligotm = os.path.join(src, "oligotm.c")
    if os.path.isfile(oligotm):
        found["oligotm_c"] = oligotm
    return found


# ---------------------------------------------------------------------------
# primer3_config ayristirma
# ---------------------------------------------------------------------------
def _data_rows(path: str) -> list[list[str]]:
    rows = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.split("#")[0].strip()
            if line:
                rows.append(line.split())
    return rows


def _value(token: str) -> float:
    return float("inf") if token == "inf" else float(token)


def _read_quad(path: str) -> dict[tuple[str, str, str, str], float]:
    """stack / stackmm bicimi: 256 satir, okuma sirasi i, ii, j, jj."""
    rows = _data_rows(path)
    if len(rows) != 256:
        raise OracleUnavailable(f"{path}: 256 satir bekleniyordu, {len(rows)} bulundu")
    table: dict[tuple[str, str, str, str], float] = {}
    n = 0
    for i in BASES:
        for ii in BASES:
            for j in BASES:
                for jj in BASES:
                    label, token = rows[n][0], rows[n][1]
                    n += 1
                    expected = f"{i}{ii}_{j}{jj}"
                    if label != expected:
                        raise OracleUnavailable(
                            f"{path}: satir {n} etiketi {label}, beklenen {expected}"
                        )
                    table[(i, ii, j, jj)] = _value(token)
    return table


def _read_dangle(path: str) -> dict[str, float]:
    """dangle bicimi: 64 satir 3' bolumu, 64 satir 5' bolumu.

    Donen anahtarlar dogrudan `params.py` sozlesmesindedir (D3:/D5:).
    Yalnizca Watson-Crick kapanis cifti tasiyan girdiler alinir; primer3
    diger kombinasyonlara 0 yazar ve bunlarin fiziksel karsiligi yoktur.
    """
    rows = _data_rows(path)
    if len(rows) != 128:
        raise OracleUnavailable(f"{path}: 128 satir bekleniyordu, {len(rows)} bulundu")
    out: dict[str, float] = {}
    n = 0
    for i in BASES:            # 3' bolumu: dangleEnthalpies3[i][k][j]
        for j in BASES:
            for k in BASES:
                label, token = rows[n][0], rows[n][1]
                n += 1
                expected = f"{i}{k}_{j}"
                if label != expected:
                    raise OracleUnavailable(
                        f"{path}: satir {n} etiketi {label}, beklenen {expected}"
                    )
                if (i, j) in WC_PAIRS:
                    out[f"D3:{k}:{i}-{j}"] = _value(token)
    for i in BASES:            # 5' bolumu: dangleEnthalpies5[i][j][k]
        for j in BASES:
            for k in BASES:
                label, token = rows[n][0], rows[n][1]
                n += 1
                expected = f"{i}_{j}{k}"
                if label != expected:
                    raise OracleUnavailable(
                        f"{path}: satir {n} etiketi {label}, beklenen {expected}"
                    )
                if (j, i) in WC_PAIRS:
                    out[f"D5:{k}:{j}-{i}"] = _value(token)
    return out


def _read_loops(path: str) -> dict[str, float]:
    """loops bicimi: uzunluk, INTERNAL, BULGE, HAIRPIN."""
    rows = _data_rows(path)
    if len(rows) != 30:
        raise OracleUnavailable(f"{path}: 30 satir bekleniyordu, {len(rows)} bulundu")
    out: dict[str, float] = {}
    for row in rows:
        size = int(row[0])
        for name, token in zip(("INTERNAL", "BULGE", "HAIRPIN"), row[1:4]):
            out[f"{name}:{size}"] = _value(token)
    return out


def load_thal_oracle(config_dir: str) -> dict[str, dict]:
    """primer3_config tablolari -> params.py sozlesmesinde (kcal, cal/mol K)."""
    stack_dh = _read_quad(os.path.join(config_dir, "stack.dh"))
    stack_ds = _read_quad(os.path.join(config_dir, "stack.ds"))
    mm_dh = _read_quad(os.path.join(config_dir, "stackmm.dh"))
    mm_ds = _read_quad(os.path.join(config_dir, "stackmm.ds"))
    dan_dh = _read_dangle(os.path.join(config_dir, "dangle.dh"))
    dan_ds = _read_dangle(os.path.join(config_dir, "dangle.ds"))
    loop_dh = _read_loops(os.path.join(config_dir, "loops.dh"))
    loop_ds = _read_loops(os.path.join(config_dir, "loops.ds"))

    wc: dict[str, tuple[float, float]] = {}
    for key, dh in stack_dh.items():
        if dh == float("inf"):
            continue
        i, ii, j, jj = key
        if (i, j) not in WC_PAIRS or (ii, jj) not in WC_PAIRS:
            continue
        wc[f"{i}{ii}"] = (dh / 1000.0, stack_ds[key])

    mismatch: dict[str, tuple[float, float]] = {}
    for key, dh in mm_dh.items():
        i, ii, j, jj = key
        left_wc = (i, j) in WC_PAIRS
        right_wc = (ii, jj) in WC_PAIRS
        if left_wc == right_wc:        # ikisi de eslesmis ya da ikisi de degil
            continue
        if dh == float("inf"):
            continue
        mismatch[f"{i}{ii}/{j}{jj}"] = (dh / 1000.0, mm_ds[key])

    dangling = {
        key: (dh / 1000.0, dan_ds[key])
        for key, dh in dan_dh.items()
        if dh != float("inf")
    }

    loops: dict[str, float] = {}
    for key, dh in loop_dh.items():
        ds = loop_ds[key]
        if dh == float("inf") or ds == float("inf"):
            continue
        loops[key] = (dh - T37 * ds) / 1000.0

    return {"wc": wc, "mismatch": mismatch, "dangling": dangling, "loops": loops}


# ---------------------------------------------------------------------------
# oligotm.c ayristirma
# ---------------------------------------------------------------------------
def _c_int_matrix(source: str, name: str) -> list[list[int]]:
    match = re.search(
        r"static\s+const\s+int\s+" + re.escape(name) + r"\s*\[5\]\[5\]\s*=\s*\{(.*?)\};",
        source,
        re.S,
    )
    if not match:
        raise OracleUnavailable(f"oligotm.c: {name} bulunamadi")
    body = re.sub(r"/\*.*?\*/", "", match.group(1), flags=re.S)
    rows = re.findall(r"\{([^{}]*)\}", body)
    if len(rows) != 5:
        raise OracleUnavailable(f"oligotm.c: {name} 5 satir degil")
    return [[int(x) for x in row.replace(" ", "").split(",") if x] for row in rows]


def load_oligotm_oracle(oligotm_c: str) -> dict[str, dict]:
    """oligotm.c santalucia_auto sabitleri -> params.py birimleri."""
    with open(oligotm_c, encoding="utf-8", errors="replace") as handle:
        source = handle.read()
    dh = _c_int_matrix(source, "SantaLucia_1998_dH")
    ds = _c_int_matrix(source, "SantaLucia_1998_dS")
    wc = {}
    for a, row_h in enumerate(dh[:4]):
        for b in range(4):
            key = BASES[a] + BASES[b]
            wc[key] = (-100.0 * row_h[b] / 1000.0, -0.1 * ds[a][b])

    block = re.search(
        r"if\s*\(\s*tm_method\s*==\s*santalucia_auto\s*\)\s*\{(.*?)\}\s*else\s+if",
        source,
        re.S,
    )
    if not block:
        raise OracleUnavailable("oligotm.c: santalucia_auto blogu bulunamadi")
    body = block.group(1)
    sym = re.search(r"sym\s*==\s*1\s*\)\s*\{\s*ds\s*\+=\s*(-?\d+)\s*;", body)
    at = re.search(
        r"case\s*'A'\s*:\s*case\s*'T'\s*:\s*ds\s*\+=\s*(-?\d+)\s*;\s*dh\s*\+=\s*(-?\d+)\s*;",
        body,
    )
    gc = re.search(
        r"case\s*'C'\s*:\s*case\s*'G'\s*:\s*ds\s*\+=\s*(-?\d+)\s*;\s*dh\s*\+=\s*(-?\d+)\s*;",
        body,
    )
    if not (sym and at and gc):
        raise OracleUnavailable("oligotm.c: baslatma/simetri terimleri okunamadi")
    init = {
        "AT": (-100.0 * int(at.group(2)) / 1000.0, -0.1 * int(at.group(1))),
        "GC": (-100.0 * int(gc.group(2)) / 1000.0, -0.1 * int(gc.group(1))),
    }
    symmetry = (0.0, -0.1 * int(sym.group(1)))
    return {"wc": wc, "init": init, "symmetry": symmetry}


# ---------------------------------------------------------------------------
# Karsilastirma
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CrossCheckEntry:
    """Tek bir parametre girdisinin iki kaynak arasindaki farki."""

    table_id: str
    oracle: str
    key: str
    quantity: str          # dH | dS | dG37
    ours: float | None
    theirs: float | None
    tolerance: float
    status: str            # MATCH | SOURCE_NOTE | DIFF | MISSING_OURS | MISSING_ORACLE

    @property
    def delta(self) -> float | None:
        if self.ours is None or self.theirs is None:
            return None
        return self.ours - self.theirs

    def as_row(self) -> dict[str, str]:
        delta = self.delta
        return {
            "table_id": self.table_id,
            "oracle": self.oracle,
            "key": self.key,
            "quantity": self.quantity,
            "ours": "" if self.ours is None else f"{self.ours:.4f}",
            "oracle_value": "" if self.theirs is None else f"{self.theirs:.4f}",
            "delta": "" if delta is None else f"{delta:+.4f}",
            "tolerance": f"{self.tolerance:.4f}",
            "status": self.status,
        }


def _compare(table_id, oracle, key, quantity, ours, theirs, tol) -> CrossCheckEntry:
    if ours is None:
        status = "MISSING_OURS"
    elif theirs is None:
        status = "MISSING_ORACLE"
    elif abs(ours - theirs) <= tol:
        status = "MATCH"
    else:
        status = "DIFF"
    return CrossCheckEntry(table_id, oracle, key, quantity, ours, theirs, tol, status)


def _rotate_mismatch_key(key: str) -> str:
    """XY/ZW yigini 180 derece donunce WZ/YX olur; ayni fiziksel yigin."""
    top, bottom = key.split("/")
    return f"{bottom[1]}{bottom[0]}/{top[1]}{top[0]}"


def crosscheck_wc(oracle_name: str, oracle_wc: dict) -> list[CrossCheckEntry]:
    """10 bagimsiz WC degeri; oracle 16 din. kombinasyonun hepsini tasir."""
    data = WC_SANTALUCIA_1998.data["nn"]
    table_id = WC_SANTALUCIA_1998.provenance.table_id
    out = []
    for first in BASES:
        for second in BASES:
            key = first + second
            revcomp = COMPLEMENT[second] + COMPLEMENT[first]
            mine = data.get(key) or data.get(revcomp)
            theirs = oracle_wc.get(key)
            out.append(_compare(table_id, oracle_name, key, "dH",
                                None if mine is None else mine[0],
                                None if theirs is None else theirs[0],
                                TOL_DH_KCAL))
            out.append(_compare(table_id, oracle_name, key, "dS",
                                None if mine is None else mine[1],
                                None if theirs is None else theirs[1],
                                TOL_DS_CAL))
    return out


def crosscheck_init(oracle_name: str, oracle: dict) -> list[CrossCheckEntry]:
    table_id = WC_SANTALUCIA_1998.provenance.table_id
    out = []
    init = WC_SANTALUCIA_1998.data["init"]
    for key in ("AT", "GC"):
        mine = init.get(key)
        theirs = oracle["init"].get(key)
        out.append(_compare(table_id, oracle_name, f"init:{key}", "dH",
                            None if mine is None else mine[0],
                            None if theirs is None else theirs[0], TOL_DH_KCAL))
        out.append(_compare(table_id, oracle_name, f"init:{key}", "dS",
                            None if mine is None else mine[1],
                            None if theirs is None else theirs[1], TOL_DS_CAL))
    mine_sym = WC_SANTALUCIA_1998.data["symmetry"]
    theirs_sym = oracle["symmetry"]
    out.append(_compare(table_id, oracle_name, "symmetry", "dH",
                        mine_sym[0], theirs_sym[0], TOL_DH_KCAL))
    out.append(_compare(table_id, oracle_name, "symmetry", "dS",
                        mine_sym[1], theirs_sym[1], TOL_DS_CAL))
    return out


def crosscheck_mismatch(oracle_name: str, oracle_mm: dict) -> list[CrossCheckEntry]:
    data = MISMATCH_ALLAWI_PEYRET.data
    table_id = MISMATCH_ALLAWI_PEYRET.provenance.table_id
    out = []
    for key in sorted(data):
        mine = data[key]
        theirs = oracle_mm.get(key) or oracle_mm.get(_rotate_mismatch_key(key))
        out.append(_compare(table_id, oracle_name, key, "dH", mine[0],
                            None if theirs is None else theirs[0], TOL_DH_KCAL))
        out.append(_compare(table_id, oracle_name, key, "dS", mine[1],
                            None if theirs is None else theirs[1], TOL_DS_CAL))
    # Oracle'da olup bizde olmayan fiziksel yiginlar (donme esdegeri dahil).
    covered = set(data) | {_rotate_mismatch_key(k) for k in data}
    for key in sorted(set(oracle_mm) - covered):
        out.append(_compare(table_id, oracle_name, key, "dH", None,
                            oracle_mm[key][0], TOL_DH_KCAL))
    return out


def crosscheck_dangling(oracle_name: str, oracle_dan: dict) -> list[CrossCheckEntry]:
    data = DANGLING_BOMMARITO_2000.data
    table_id = DANGLING_BOMMARITO_2000.provenance.table_id
    out = []
    for key in sorted(set(data) | set(oracle_dan)):
        mine = data.get(key)
        theirs = oracle_dan.get(key)
        out.append(_compare(table_id, oracle_name, key, "dH",
                            None if mine is None else mine[0],
                            None if theirs is None else theirs[0], TOL_DH_KCAL))
        out.append(_compare(table_id, oracle_name, key, "dS",
                            None if mine is None else mine[1],
                            None if theirs is None else theirs[1], TOL_DS_CAL))
    return out


def crosscheck_loops(oracle_name: str, oracle_loops: dict) -> list[CrossCheckEntry]:
    data = DNA_LOOP_SANTALUCIA_HICKS_2004.data
    table_id = DNA_LOOP_SANTALUCIA_HICKS_2004.provenance.table_id
    out = []
    for key in sorted(data):
        mine = data[key][0]
        theirs = oracle_loops.get(key)
        out.append(_compare(table_id, oracle_name, key, "dG37", mine, theirs,
                            TOL_DG37_KCAL))
    return out


def crosscheck_terminal(config_dir: str) -> list[CrossCheckEntry]:
    table = TERMINAL_MISMATCH_SANTALUCIA_PEYRET
    dh = _read_quad(os.path.join(config_dir, "tstack2.dh"))
    ds = _read_quad(os.path.join(config_dir, "tstack2.ds"))
    # Document only these exact source discrepancies, and surface any new drift.
    source_notes = {
        ("AC/TC", "dS"): (0.5, -0.5),
        ("GG/CG", "dH"): (-6.91488, -0.7),
    }
    rows = []
    for key, (ours_h, ours_s) in sorted(table.data.items()):
        i, ii, j, jj = key.replace("/", "")
        coords = i, ii, j, jj
        for qty, ours, oracle, tol in (
            ("dH", ours_h, dh[coords] / 1000.0, TOL_DH_KCAL),
            ("dS", ours_s, ds[coords], TOL_DS_CAL),
        ):
            row = _compare(table.provenance.table_id, ORACLE_TERMINAL, key, qty, ours, oracle, tol)
            known = source_notes.get((key, qty))
            if known and abs(ours-known[0]) < 1e-9 and abs(oracle-known[1]) < 1e-9:
                row = replace(row, status="SOURCE_NOTE")
            rows.append(row)
    return rows


def run_crosscheck(primer3_src: str | None = None) -> list[CrossCheckEntry]:
    """Mevcut tum oracle'lara karsi tam denetim. Kaynak yoksa hata verir."""
    sources = oracle_sources(primer3_src)
    if not sources:
        raise OracleUnavailable(
            "primer3 kaynagi bulunamadi. OLIGOSCOUT_PRIMER3_SRC ortam "
            "degiskenini primer3 'src' klasorune ayarlayin."
        )
    entries: list[CrossCheckEntry] = []
    if "config_dir" in sources:
        thal = load_thal_oracle(sources["config_dir"])
        entries += crosscheck_wc(ORACLE_THAL, thal["wc"])
        entries += crosscheck_mismatch(ORACLE_THAL, thal["mismatch"])
        entries += crosscheck_dangling(ORACLE_THAL, thal["dangling"])
        entries += crosscheck_loops(ORACLE_THAL, thal["loops"])
    if "oligotm_c" in sources:
        tm = load_oligotm_oracle(sources["oligotm_c"])
        entries += crosscheck_wc(ORACLE_OLIGOTM, tm["wc"])
        entries += crosscheck_init(ORACLE_OLIGOTM, tm)
    if "terminal_dir" in sources:
        entries += crosscheck_terminal(sources["terminal_dir"])
    return entries


def summarise(entries: list[CrossCheckEntry]) -> list[dict[str, str]]:
    """Tablo x oracle kirilimi: karsilastirilan, eslesen, sapan, eksik."""
    keys = []
    buckets: dict[tuple[str, str], list[CrossCheckEntry]] = {}
    for entry in entries:
        bucket = (entry.table_id, entry.oracle)
        if bucket not in buckets:
            buckets[bucket] = []
            keys.append(bucket)
        buckets[bucket].append(entry)
    rows = []
    for table_id, oracle in keys:
        group = buckets[(table_id, oracle)]
        compared = [e for e in group if e.status in ("MATCH", "DIFF", "SOURCE_NOTE")]
        deltas = [abs(e.delta) for e in compared if e.delta is not None]
        rows.append({
            "table_id": table_id,
            "oracle": oracle,
            "compared": str(len(compared)),
            "match": str(sum(1 for e in group if e.status == "MATCH")),
            "diff": str(sum(1 for e in group if e.status == "DIFF")),
            "source_notes": str(sum(1 for e in group if e.status == "SOURCE_NOTE")),
            "missing_ours": str(sum(1 for e in group if e.status == "MISSING_OURS")),
            "missing_oracle": str(sum(1 for e in group if e.status == "MISSING_ORACLE")),
            "max_abs_delta": f"{max(deltas):.4f}" if deltas else "",
        })
    return rows
