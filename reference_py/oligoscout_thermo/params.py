"""Termodinamik parametre tablolari ve kunyeleri.

Tasarim kurali (QC spesifikasyonu, bolum 1):
  * Parametreler *deneysel veridir*, yazilim degildir. Birincil makalelerden
    aktarilir; baska yazilimlarin parametre dosyalari kopyalanmaz.
  * Her tablo kunyesini (kaynak, tablo no, birim, kosul) ve icerigin
    SHA-256'sini tasir. Bu checksum her rapor satirina yazilir (REP-02).
  * Her tablonun bir `verification_status` alani vardir. Bir tablo ancak
    birincil PDF'ten ikinci bir kisi tarafindan bagimsiz aktarilip diff
    alindiktan sonra "verified" olur. O ana kadar "transcribed_unverified".
  * Ayrica her tablonun bir `cross_check` alani vardir: ayni birincil
    makalelerden BAGIMSIZ olarak aktarilmis bir yazilimin tablolarina karsi
    yapilan diferansiyel denetimin sonucu. Bu denetim transkripsiyon
    hatasini yakalar ama ikinci insan aktariminin yerine GECMEZ; bu yuzden
    `verification_status`u degistirmez. Bkz. `crosscheck.py`,
    `cli/param_crosscheck.py`.

Durum (bu dosyanin yazildigi an):
  * Watson-Crick NN tablosu  : transcribed_unverified
      SantaLucia 1998 PNAS makalesinin PMC kaydinda yalnizca ozet bulunuyor,
      Tablo 1 govdesi acik erisimde degil. Degerler literaturden aktarildi ve
      bagimsiz bir uygulamaya karsi diferansiyel testle sinaniyor
      (tests/test_thermo.py; oracle yalniz gelistirmede kullanilir, urunde
      yer almaz). Birincil PDF elinize gectiginde ikinci aktarimi yapip
      durumu "verified"e cekin.
  * Uyumsuzluk (mismatch) tablosu : transcribed_unverified
      Allawi & SantaLucia 1997/1998 ve Peyret 1999 tablolari.
  * Sarkan uc tablosu : transcribed_unverified
      Bommarito 2000 Tablo 2 aktarimi OCR ile karsilastirildi.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any

R_CAL = 1.9872  # gaz sabiti, cal / (mol K)
KELVIN = 273.15


# ---------------------------------------------------------------------------
# Kunye
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TableProvenance:
    """Bir parametre tablosunun kunyesi."""

    table_id: str
    citation: str
    table_in_source: str
    units: str
    conditions: str
    verification_status: str  # verified | transcribed_unverified | pending
    notes: str = ""
    cross_check: str = ""
    # cross_check, `verification_status`tan AYRI bir alandir ve onun yerine
    # GECMEZ. Bagimsiz bir ikinci aktarima (ornegin primer3'un tablolari)
    # karsi yapilan diferansiyel denetimin sonucunu tasir:
    #   "<oracle> <eslesen>/<karsilastirilan> max|d|=<deger>" biciminde.
    # Uretimi ve dogrulanmasi: cli/param_crosscheck.py,
    # tests/test_param_crosscheck.py. Bir tablo ancak birincil PDF'ten
    # ikinci bir KISI aktarip diff aldiktan sonra "verified" olur.

    def as_dict(self) -> dict[str, str]:
        return {
            "table_id": self.table_id,
            "citation": self.citation,
            "table_in_source": self.table_in_source,
            "units": self.units,
            "conditions": self.conditions,
            "verification_status": self.verification_status,
            "notes": self.notes,
            "cross_check": self.cross_check,
        }


class ParameterTableUnavailable(RuntimeError):
    """Gerekli parametre tablosu aktarilmamis."""


@dataclass(frozen=True)
class ParameterTable:
    """Veri + kunye + checksum."""

    provenance: TableProvenance
    data: dict[str, Any] = field(default_factory=dict)

    @property
    def sha256(self) -> str:
        """Kanonik JSON uzerinden icerik checksum'i (REP-02)."""
        blob = json.dumps(
            {"provenance": self.provenance.as_dict(), "data": self.data},
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("ascii")
        return hashlib.sha256(blob).hexdigest()

    @property
    def is_available(self) -> bool:
        return bool(self.data) and self.provenance.verification_status != "pending"

    def require(self) -> dict[str, Any]:
        if not self.is_available:
            raise ParameterTableUnavailable(
                f"'{self.provenance.table_id}' tablosu yuklu degil "
                f"(durum: {self.provenance.verification_status}). "
                f"Kaynak: {self.provenance.citation}. "
                "Degerleri birincil makaleden aktarmadan bu hesap yapilamaz."
            )
        return self.data


# ---------------------------------------------------------------------------
# Watson-Crick nearest-neighbor (birlesik set)
# ---------------------------------------------------------------------------
# Anahtar: 5'->3' dinukleotit (ust iplik). Karsilik gelen alt iplik ters
# tamamlayicidir, bu yuzden 10 bagimsiz deger yeterlidir.
# Deger: (dH kcal/mol, dS cal/(mol K)), 1 M NaCl.
_WC_NN: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
}

# Her ucta, o ucun baz ciftine gore uygulanan baslatma terimi.
_WC_INIT: dict[str, tuple[float, float]] = {
    "AT": (2.3, 4.1),   # ucta A veya T
    "GC": (0.1, -2.8),  # ucta G veya C
}

# Kendine tamamlayici dupleksler icin simetri duzeltmesi.
_WC_SYMMETRY: tuple[float, float] = (0.0, -1.4)

WC_SANTALUCIA_1998 = ParameterTable(
    provenance=TableProvenance(
        table_id="wc_nn_santalucia_1998_unified",
        citation="SantaLucia J Jr. PNAS 1998;95:1460-1465. doi:10.1073/pnas.95.4.1460",
        table_in_source="Table 1 (unified oligonucleotide NN set)",
        units="dH kcal/mol, dS cal/(mol K)",
        conditions="1 M NaCl, pH 7",
        verification_status="transcribed_unverified",
        notes=(
            "Birincil PDF acik erisimde degil; PMC kaydinda yalniz ozet var. "
            "Bagimsiz bir uygulamaya karsi diferansiyel test uygulanir."
        ),
        cross_check="primer3_thal 32/32 max|d|=0.0000; "
                    "primer3_oligotm 38/38 max|d|=0.0000",
    ),
    data={
        "nn": {k: list(v) for k, v in _WC_NN.items()},
        "init": {k: list(v) for k, v in _WC_INIT.items()},
        "symmetry": list(_WC_SYMMETRY),
    },
)


# ---------------------------------------------------------------------------
# Ic tekli uyumsuzluklar
# Anahtar: XY/ZW, ust iplik 5'->3', alt iplik 3'->5'.
# Deger: [dH kcal/mol, dS cal/(mol K)].
MISMATCH_ALLAWI_PEYRET = ParameterTable(
    provenance=TableProvenance(
        table_id="mismatch_nn_allawi_peyret",
        citation=(
            "Allawi & SantaLucia, Biochemistry 1997;36:10581 (G.T); "
            "1998;37:2170 (G.A); 1998;37:9435 (A.C); NAR 1998;26:2694 (C.T); "
            "Peyret et al., Biochemistry 1999;38:3468 (A.A, C.C, G.G, T.T)"
        ),
        table_in_source="G.T Tablo 5; G.A, A.C, C.T Tablo 4; Peyret Tablo 2",
        units="dH kcal/mol, dS cal/(mol K)",
        conditions="1 M NaCl, pH 7",
        verification_status="transcribed_unverified",
        notes="48 kanonik komsu yigin. Yalniz tekli ic uyumsuzluk.",
        cross_check="primer3_thal 96/96 max|d|=0.0000",
    ),
    data={
        "AA/AT": [4.7, 12.9],
        "AA/CT": [7.6, 20.2],
        "AA/GT": [3.0, 7.4],
        "AA/TA": [1.2, 1.7],
        "AA/TC": [2.3, 4.6],
        "AA/TG": [-0.6, -2.3],
        "AC/AG": [-2.9, -9.8],
        "AC/CG": [-0.7, -3.8],
        "AC/GG": [0.5, 3.2],
        "AC/TA": [5.3, 14.6],
        "AC/TC": [0.0, -4.4],
        "AC/TT": [0.7, 0.2],
        "AG/AC": [-0.9, -4.2],
        "AG/CC": [0.6, -0.6],
        "AG/GC": [-4.0, -13.2],
        "AG/TA": [-0.7, -2.3],
        "AG/TG": [-3.1, -9.5],
        "AG/TT": [1.0, 0.9],
        "AT/TC": [-1.2, -6.2],
        "AT/TG": [-2.5, -8.3],
        "AT/TT": [-2.7, -10.8],
        "CA/AT": [3.4, 8.0],
        "CA/CT": [6.1, 16.4],
        "CA/GC": [1.9, 3.7],
        "CA/GG": [-0.7, -2.3],
        "CA/TT": [1.0, 0.7],
        "CC/AG": [5.2, 14.2],
        "CC/CG": [3.6, 8.9],
        "CC/GC": [-1.5, -7.2],
        "CC/GT": [-0.8, -4.5],
        "CC/TG": [5.2, 13.5],
        "CG/GG": [-4.9, -15.3],
        "CG/GT": [-4.1, -11.7],
        "CG/TC": [-1.5, -6.1],
        "CT/GG": [-2.8, -8.0],
        "CT/GT": [-5.0, -15.8],
        "GA/AT": [0.7, 0.7],
        "GA/CG": [-0.6, -1.0],
        "GA/GT": [1.6, 3.6],
        "GA/TT": [-1.3, -5.3],
        "GC/CT": [2.3, 5.4],
        "GC/GG": [-6.0, -15.8],
        "GC/TG": [-4.4, -12.3],
        "GG/CT": [3.3, 10.4],
        "GT/CT": [-2.2, -8.4],
        "TA/CT": [1.2, 0.7],
        "TA/GT": [-0.1, -1.7],
        "TA/TT": [0.2, -1.5],
    },
)


# ---------------------------------------------------------------------------
# Sarkan uclar
# Anahtar: D3/D5:<sarkan baz>:<sarkan iplikteki baz>-<karsi baz>.
# Deger: [dH kcal/mol, dS cal/(mol K)], Bommarito Tablo 2.
# ---------------------------------------------------------------------------
DANGLING_BOMMARITO_2000 = ParameterTable(
    provenance=TableProvenance(
        table_id="dangling_ends_bommarito_2000",
        citation=(
            "Bommarito, Peyret & SantaLucia, NAR 2000;28:1929-1934. "
            "doi:10.1093/nar/28.9.1929"
        ),
        table_in_source="Table 2 (5' ve 3' sarkan uc NN parametreleri)",
        units="dH kcal/mol, dS cal/(mol K)",
        conditions="1 M NaCl, pH 7",
        verification_status="transcribed_unverified",
        notes="Tablo 2 gorsellerinden aktarildi ve OCR ile karsilastirildi.",
        cross_check="primer3_thal 64/64 max|d|=0.0000",
    ),
    data={
        "D3:A:A-T": [-0.5, -1.1],
        "D3:A:C-G": [-5.9, -16.5],
        "D3:A:G-C": [-2.1, -3.9],
        "D3:A:T-A": [-0.7, -0.8],
        "D3:C:A-T": [4.7, 14.2],
        "D3:C:C-G": [-2.6, -7.4],
        "D3:C:G-C": [-0.2, -0.1],
        "D3:C:T-A": [4.4, 14.9],
        "D3:G:A-T": [-4.1, -13.1],
        "D3:G:C-G": [-3.2, -10.4],
        "D3:G:G-C": [-3.9, -11.2],
        "D3:G:T-A": [-1.6, -3.6],
        "D3:T:A-T": [-3.8, -12.6],
        "D3:T:C-G": [-5.2, -15.0],
        "D3:T:G-C": [-4.4, -13.1],
        "D3:T:T-A": [2.9, 10.4],
        "D5:A:A-T": [0.2, 2.3],
        "D5:A:C-G": [-6.3, -17.1],
        "D5:A:G-C": [-3.7, -10.0],
        "D5:A:T-A": [-2.9, -7.6],
        "D5:C:A-T": [0.6, 3.3],
        "D5:C:C-G": [-4.4, -12.6],
        "D5:C:G-C": [-4.0, -11.9],
        "D5:C:T-A": [-4.1, -13.0],
        "D5:G:A-T": [-1.1, -1.6],
        "D5:G:C-G": [-5.1, -14.0],
        "D5:G:G-C": [-3.9, -10.9],
        "D5:G:T-A": [-4.2, -15.0],
        "D5:T:A-T": [-6.9, -20.0],
        "D5:T:C-G": [-4.0, -10.9],
        "D5:T:G-C": [-4.9, -13.8],
        "D5:T:T-A": [-0.2, -0.5],
    },
)


# ---------------------------------------------------------------------------
# DNA dongulerinin uzunluga bagli dG37 artimlari
# ---------------------------------------------------------------------------
DNA_LOOP_SANTALUCIA_HICKS_2004 = ParameterTable(
    provenance=TableProvenance(
        table_id="dna_loop_length_santalucia_hicks_2004",
        citation=("SantaLucia & Hicks, Annu Rev Biophys Biomol Struct 2004;33:415-440. "
                  "doi:10.1146/annurev.biophys.32.110601.141800"),
        table_in_source="Table 4 (length terms), p. 427",
        units="dG37 kcal/mol",
        conditions="37 C, 1 M NaCl",
        verification_status="transcribed_unverified",
        notes="Tablo 4 uzunluk terimleri. dH=0 kaynak varsayimidir.",
        cross_check="primer3_thal 47/47 max|d|=0.0030",
    ),
    data={
        "BULGE:1": [4.0],
        "BULGE:10": [4.3],
        "BULGE:12": [4.5],
        "BULGE:14": [4.8],
        "BULGE:16": [5.0],
        "BULGE:18": [5.2],
        "BULGE:2": [2.9],
        "BULGE:20": [5.3],
        "BULGE:25": [5.6],
        "BULGE:3": [3.1],
        "BULGE:30": [5.9],
        "BULGE:4": [3.2],
        "BULGE:5": [3.3],
        "BULGE:6": [3.5],
        "BULGE:7": [3.7],
        "BULGE:8": [3.9],
        "BULGE:9": [4.1],
        "HAIRPIN:10": [4.6],
        "HAIRPIN:12": [5.0],
        "HAIRPIN:14": [5.1],
        "HAIRPIN:16": [5.3],
        "HAIRPIN:18": [5.5],
        "HAIRPIN:20": [5.7],
        "HAIRPIN:25": [6.1],
        "HAIRPIN:3": [3.5],
        "HAIRPIN:30": [6.3],
        "HAIRPIN:4": [3.5],
        "HAIRPIN:5": [3.3],
        "HAIRPIN:6": [4.0],
        "HAIRPIN:7": [4.2],
        "HAIRPIN:8": [4.3],
        "HAIRPIN:9": [4.5],
        "INTERNAL:10": [4.9],
        "INTERNAL:12": [5.2],
        "INTERNAL:14": [5.4],
        "INTERNAL:16": [5.6],
        "INTERNAL:18": [5.8],
        "INTERNAL:20": [5.9],
        "INTERNAL:25": [6.3],
        "INTERNAL:3": [3.2],
        "INTERNAL:30": [6.6],
        "INTERNAL:4": [3.6],
        "INTERNAL:5": [4.0],
        "INTERNAL:6": [4.4],
        "INTERNAL:7": [4.6],
        "INTERNAL:8": [4.8],
        "INTERNAL:9": [4.9],
    },
)


# ---------------------------------------------------------------------------
# Terminal mismatch NN increments transcribed from the primary PDF, pp. 52–54.
# First base pair is Watson-Crick, second base pair is the terminal mismatch.
TERMINAL_MISMATCH_SANTALUCIA_PEYRET = ParameterTable(
    provenance=TableProvenance(
        table_id="terminal_mismatch_santalucia_peyret_2001",
        citation="SantaLucia & Peyret (2001), WO2001094611A2",
        table_in_source="Tables 2–3, printed pages 52–54",
        units="dH kcal/mol, dS cal/(mol K)",
        conditions="1 M NaCl",
        verification_status="transcribed_unverified",
        notes=("Terminal increment on WC core. GG/CG: printed dH=-0.7 is inconsistent "
               "with dS=-19.2 and dG37=-0.96. Model dH is reconstructed from dG37 and dS."),
        cross_check="primer3_tstack2 94/96 max|d|=6.2149",
    ),
    data={
        "AA/TA": [-3.1, -7.8], "TA/AA": [-2.5, -6.3], "CA/GA": [-4.3, -10.7], "GA/CA": [-8.0, -22.5],
        "AC/TC": [-0.1, 0.5], "TC/AC": [-0.7, -1.3], "CC/GC": [-2.1, -5.1], "GC/CC": [-3.9, -10.6],
        "AG/TG": [-1.1, -2.1], "TG/AG": [-1.1, -2.7], "CG/GG": [-3.8, -9.5],
        "GG/CG": [-0.96 + 310.15 * -19.2 / 1000, -19.2],
        "AT/TT": [-2.4, -6.5], "TT/AT": [-3.2, -8.9], "CT/GT": [-6.1, -16.9], "GT/CT": [-7.4, -21.2],
        "AA/TC": [-1.6, -4.0], "AC/TA": [-1.8, -3.8], "CA/GC": [-2.6, -5.9], "CC/GA": [-2.7, -6.0],
        "GA/CC": [-5.0, -13.8], "GC/CA": [-3.2, -7.1], "TA/AC": [-2.3, -5.9], "TC/AA": [-2.7, -7.0],
        "AC/TT": [-0.9, -1.7], "AT/TC": [-2.3, -6.3], "CC/GT": [-3.2, -8.0], "CT/GC": [-3.9, -10.6],
        "GC/CT": [-4.9, -13.5], "GT/CC": [-3.0, -7.8], "TC/AT": [-2.5, -6.3], "TT/AC": [-0.7, -1.2],
        "AA/TG": [-1.9, -4.4], "AG/TA": [-2.5, -5.9], "CA/GG": [-3.9, -9.6], "CG/GA": [-6.0, -15.5],
        "GA/CG": [-4.3, -11.1], "GG/CA": [-4.6, -11.4], "TA/AG": [-2.0, -4.7], "TG/AA": [-2.4, -5.8],
        "AG/TT": [-3.2, -8.7], "AT/TG": [-3.5, -9.4], "CG/GT": [-3.8, -9.0], "CT/GG": [-6.6, -18.7],
        "GG/CT": [-5.7, -15.9], "GT/CG": [-5.9, -16.1], "TG/AT": [-3.9, -10.5], "TT/AG": [-3.6, -9.8],
    },
)


# Kayit defteri
# ---------------------------------------------------------------------------
# DNA/RNA hibrit NN (RT-05) — AKTARIM BEKLIYOR
#
# RT primerinin RNA kalibina baglanmasi bir DNA/RNA hibrit dupleksidir ve
# saf DNA NN tablosuyla hesaplanamaz. Sugimoto 1995 16 dinukleotit hibrit
# parametresi verir. Bu tablo BILEREK BOS birakilmistir: degerleri hafizadan
# yazmak yeni bir transkripsiyon hatasi riskidir. Birincil PDF elinizdedir
# (kaynak listesi). Her anahtari [dH_kcal/mol, dS_cal/(mol K)] olarak doldurun,
# `verification_status`u "transcribed_unverified" yapin; tablo bos oldugu
# surece `require()` hata verir ve RT-05 "pending" raporlar.
#
# Anahtar: DNA ipliginin 5'->3' dinukleotiti (RNA ipligi tersi tamamlayici).
# Doldurulacak 16 anahtar asagida (hepsi bos).
HYBRID_SUGIMOTO_1995 = ParameterTable(
    provenance=TableProvenance(
        table_id="dna_rna_hybrid_nn_sugimoto_1995",
        citation=("Sugimoto et al., Biochemistry 1995;34:11211-11216. "
                  "doi:10.1021/bi00035a029"),
        table_in_source="Table (16 DNA/RNA hibrit NN parametresi)",
        units="dH kcal/mol, dS cal/(mol K)",
        conditions="1 M NaCl, pH 7",
        verification_status="transcribed_unverified",
        notes=("DNA ipligi 5'->3' dinukleotit anahtari (RNA ters tamamlayici). "
               "Degerler Biopython Bio.SeqUtils.MeltingTemp.R_DNA_NN1 "
               "tablosundan aktarildi; bu tablo Sugimoto et al. 1995'i "
               "kaynak gosterir. Pay-anahtar eslemesi (ic iplik 5'->3') "
               "Biopython'un tam-anahtar NN toplamiyla birebir dogrulandi. "
               "Baslatma terimi init=(dH 1.9, dS -3.9) ayrica uygulanir "
               "(HYBRID_SUGIMOTO_1995_INIT). IKINCI insan aktarimiyla birincil "
               "PDF'e karsi dogrulama HALA gerekir (verified degil)."),
    ),
    data={
        "AA": [-7.8, -21.9],
        "AC": [-5.9, -12.3],
        "AG": [-9.1, -23.5],
        "AT": [-8.3, -23.9],
        "CA": [-9.0, -26.1],
        "CC": [-9.3, -23.2],
        "CG": [-16.3, -47.1],
        "CT": [-7.0, -19.7],
        "GA": [-5.5, -13.5],
        "GC": [-8.0, -17.1],
        "GG": [-12.8, -31.9],
        "GT": [-7.8, -21.6],
        "TA": [-7.8, -23.2],
        "TC": [-8.6, -22.9],
        "TG": [-10.4, -28.4],
        "TT": [-11.5, -36.4],
    },
)

# Sugimoto 1995 DNA/RNA hibrit baslatma terimi (dH kcal/mol, dS cal/(mol K)).
# Toplam dupleks dH/dS'e BIR KEZ eklenir. Kaynak: Biopython R_DNA_NN1["init"].
HYBRID_SUGIMOTO_1995_INIT = (1.9, -3.9)

# Transkripsiyon sablonu: aktarim yapan kisi bu 16 anahtari doldurur.
HYBRID_NN_KEYS = ("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                  "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")


# ---------------------------------------------------------------------------
REGISTRY: dict[str, ParameterTable] = {
    t.provenance.table_id: t
    for t in (
        WC_SANTALUCIA_1998,
        MISMATCH_ALLAWI_PEYRET,
        DANGLING_BOMMARITO_2000,
        DNA_LOOP_SANTALUCIA_HICKS_2004,
        TERMINAL_MISMATCH_SANTALUCIA_PEYRET,
        HYBRID_SUGIMOTO_1995,
    )
}


# ---------------------------------------------------------------------------
# Beyan edilmis model sabitleri
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ModelConstant:
    """Literatur tablosu OLMAYAN, acikca beyan edilmis model sabiti.

    Parametre tablolarindan ayri tutulur: bunlar bir makaleden aktarilmis
    olcum degil, modelleme kararidir. Bu yuzden `verification_status`
    tasimazlar; yerine `rationale` ve `status` tasirlar. Yine de
    `registry_digest()`e girerler, cunku degeri degistirmek uretilen
    sayilari degistirir ve kayitlardaki `param_set_sha256` bunu
    gostermelidir (REP-02).
    """

    const_id: str
    value: float
    units: str
    applies_to: str
    rationale: str
    status: str = "project_modelling_constant"

    def as_dict(self) -> dict[str, str]:
        return {
            "const_id": self.const_id,
            "value": repr(self.value),
            "units": self.units,
            "applies_to": self.applies_to,
            "rationale": self.rationale,
            "status": self.status,
        }


AT_BULGE_CLOSURE = ModelConstant(
    const_id="at_bulge_closure_dg37",
    value=0.5,
    units="dG37 kcal/mol, kapanan A/T cifti basina",
    applies_to="align.py tek bazlik ic bulge (DIM-01..DIM-04, DIM-06, DIM-08)",
    rationale=(
        "SantaLucia & Hicks 2004 A-T kapanis cezasini SAC TOKASI dongusu "
        "icin tanimlar (denklem 8; fold.py'de oyle kullanilir). Bulge "
        "kapanisina ayni buyuklukte bir ceza uygulamak LITERATURDEN "
        "AKTARILMIS bir deger degil, bu projenin modelleme kararidir. "
        "Islak-lab kalibrasyonu beklemektedir (REP-03)."
    ),
)

MODEL_CONSTANTS: dict[str, ModelConstant] = {
    c.const_id: c for c in (AT_BULGE_CLOSURE,)
}


def model_constants_summary() -> list[dict[str, str]]:
    return [MODEL_CONSTANTS[k].as_dict() for k in sorted(MODEL_CONSTANTS)]


def registry_digest() -> str:
    """Tablolarin ve model sabitlerinin birlesik checksum'i (REP-02)."""
    joined = "|".join(f"{k}:{REGISTRY[k].sha256}" for k in sorted(REGISTRY))
    joined += "||" + "|".join(
        f"{k}:{MODEL_CONSTANTS[k].value!r}" for k in sorted(MODEL_CONSTANTS)
    )
    return hashlib.sha256(joined.encode("ascii")).hexdigest()


def registry_summary() -> list[dict[str, str]]:
    rows = []
    for key in sorted(REGISTRY):
        t = REGISTRY[key]
        row = t.provenance.as_dict()
        row["sha256"] = t.sha256
        nested_entries = sum(len(v) for v in t.data.values() if isinstance(v, dict))
        row["entries"] = str(nested_entries or len(t.data))
        rows.append(row)
    return rows
