"""Custom exceptions raised by flexres."""

from __future__ import annotations


class FlexresError(Exception):
    """Base class for expected flexres failures."""


class ReferenceStructureError(FlexresError):
    """The reference structure failed validation."""


class ChainNotFoundError(FlexresError):
    """A requested chain was not found."""


class LigandNotFoundError(FlexresError):
    """A requested ligand was not found."""


class AmbiguousLigandError(FlexresError):
    """A ligand specification matched multiple residues."""


class UniProtMappingError(FlexresError):
    """A UniProt mapping could not be found."""


class AmbiguousUniProtMappingError(FlexresError):
    """A chain maps to more than one UniProt accession."""


class StructureDownloadError(FlexresError):
    """A structure download failed."""


class StructureParseError(FlexresError):
    """A structure could not be parsed."""


class AlignmentError(FlexresError):
    """A structural alignment failed."""
