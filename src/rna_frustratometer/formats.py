from __future__ import annotations

from enum import StrEnum


class SupportedFormats(StrEnum):
    PDB = "PDB"
    CIF = "CIF"
    MMCIF = "MMCIF"

    def normalize_format(self) -> SupportedFormats:
        if self == SupportedFormats.CIF:
            return SupportedFormats.MMCIF
        return self
