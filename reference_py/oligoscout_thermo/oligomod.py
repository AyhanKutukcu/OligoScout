"""Modifiye oligo notasyonunu ayristirir (LNA, 5'/3' isaretleri).

Girdi ornekleri:
    AGCATCTCCAGCAGCAGG+T+CA+T+AGAGGG   (+ sonraki baz LNA)
    /56-FAM/CCC+CT+C+T+C+TGAC/3IABkFQ/ (uc isaretleri /.../ arasinda)

Ayristirici DNA iskeletini (A/C/G/T), LNA konumlarini (0-tabanli) ve taninan
uc modifikasyonlarini dondurur. Termodinamik hesap DNA iskeleti uzerinden
yapilir; LNA kararlilik katkisi AYRI bir parametre tablosu gerektirir
(McTigue 2004) ve bu motorda MODELLENMEZ — bu durum acikca raporlanir.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

# Yaygin 5'/3' etiketleri (bilgi; termodinamik etkisi ayrica modellenmez).
_END_TOKENS = {
    "FAM", "56-FAM", "6-FAM", "HEX", "VIC", "TAMRA", "ROX", "CY5", "CY3",
    "IABKFQ", "3IABKFQ", "IABRQSP", "BHQ1", "BHQ2", "BHQ3", "MGB", "MGBNFQ",
    "3PHOS", "PHOS", "5PHOS", "BIOSG", "BIOTIN",
}


@dataclass(frozen=True)
class ParsedOligo:
    dna: str                       # yalniz A/C/G/T, 5'->3'
    lna_positions: frozenset[int]  # LNA bazlarin 0-tabanli konumu
    end_mods: tuple[str, ...]      # taninan uc etiketleri
    warnings: tuple[str, ...] = ()

    @property
    def lna_count(self) -> int:
        return len(self.lna_positions)

    def as_dict(self) -> dict:
        return {
            "dna": self.dna,
            "lna_positions": sorted(self.lna_positions),
            "lna_count": self.lna_count,
            "end_mods": list(self.end_mods),
            "warnings": list(self.warnings),
        }


def parse(raw: str) -> ParsedOligo:
    """Modifiye oligo dizisini ayristir. `raw` kullanici notasyonudur."""
    s = raw.strip()
    end_mods: list[str] = []
    warnings: list[str] = []

    # /.../ token'lari uc modifikasyonlaridir; cikar ve kaydet.
    for token in re.findall(r"/([^/]+)/", s):
        key = token.upper().replace(" ", "")
        if key in _END_TOKENS:
            end_mods.append(token)
        else:
            warnings.append(f"taninmayan uc modifikasyonu: {token}")
            end_mods.append(token)
    s = re.sub(r"/[^/]+/", "", s)

    dna_chars: list[str] = []
    lna: set[int] = set()
    is_lna_next = False
    for ch in s:
        if ch == "+":
            is_lna_next = True
            continue
        up = ch.upper()
        if up in "ACGT":
            if is_lna_next:
                lna.add(len(dna_chars))
            dna_chars.append(up)
            is_lna_next = False
        elif up in " -_":
            is_lna_next = False        # ayirici; yok say
        else:
            warnings.append(f"tanimsiz karakter yok sayildi: {ch}")
            is_lna_next = False
    if is_lna_next:
        warnings.append("dizi sonunda bos '+' (LNA bazi yok)")

    return ParsedOligo(dna="".join(dna_chars), lna_positions=frozenset(lna),
                       end_mods=tuple(end_mods), warnings=tuple(warnings))
